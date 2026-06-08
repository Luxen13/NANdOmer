import csv
import re
from pathlib import Path

import pysam
from tqdm import tqdm

from . import parer
from .utils import fasta_parser
from .utils import load_references as _load_references
from .utils import parse_cigar
from .threshold import determine_threshold
from .defaults import (
    CALLER_ANCHOR_LEN,
    CALLER_MAX_ERR_RATIO,
    CALLER_MIN_ALIGN_PCT_IDENT_END_MIN,
    CALLER_MIN_GAP_OVERLAP,
    CALLER_MIN_ANCHOR_PCT_IDENT,
    IUPAC_DICT,
)



tsv_fields = [
    "read_id",
    "ref_name",
    "alignment_score",
    "rel_error_ratio",
    "pct_ident",
    "align_pct_ident_end_min",
    "align_pct_ident_5",
    "align_pct_ident_3",
    "anchor_pct_ident",
    "gap_pct_ratio",
    "gap_sequence",
    "gap_seq_match",
    "len_gap_overlap_5",
    "len_gap_overlap_3",
    "cigar_string",
    "ref_begin",
    "query_begin",
]

_FAILURE_FIELDS = [
    "read_id",
    "ref_name",
    "alignment_score",
    "align_pct_ident_5",
    "align_pct_ident_3",
    "align_end_ident_ratio",
    "min_align_end_ident_ratio",
    "pct_ident",
    "min_pct_ident",
    "anchor_ident_ratio",
    "min_anchor_ident_ratio",
    "failure_reason",
]



def _count_bam_records(source_bam):
    with pysam.AlignmentFile(str(source_bam), "rb", check_sq=False) as bam:
        return sum(1 for _ in bam.fetch(until_eof=True))


def _failure_row(**overrides):
    row = {field: "N.a" for field in _FAILURE_FIELDS}
    row.update(overrides)
    return row


def _tsv_round_row(row, ndigits=3):
    rounded = {}
    for key, value in row.items():
        if isinstance(value, float):
            rounded[key] = round(value, ndigits)
        else:
            rounded[key] = value
    return rounded



def _parse_cigar_validated(cigar_string):
    ops = parse_cigar(cigar_string)
    if any(op == "M" for _, op in ops):
        raise ValueError("Unexpected M op in CIGAR")
    return ops


def _safe_filename(name):
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._")
    return safe or "motif"

def find_gap_region(reference_seq):
    """Find the span of the variable/ambiguous region using all IUPAC codes."""
    match_positions = list(re.finditer(r"[^ACGT]+", reference_seq.upper()))
    if not match_positions:
        return None, None, 0
    gap_start = match_positions[0].start()
    gap_end = match_positions[-1].end()
    gap_length = gap_end - gap_start
    return gap_start, gap_end, gap_length

def _walk_cigar(
    ops,
    query_begin,
    ref_begin,
    gap_ref_start,
    gap_ref_end,
    left_anchor_start,
    left_anchor_end,
    right_anchor_start,
    right_anchor_end,
):
    ref_pos = ref_begin
    query_pos = query_begin

    identity_vec = []
    query_gap_start = None
    query_gap_end = None
    left_anchor_vals = []
    right_anchor_vals = []
    covered_n = 0
    gap_start_col = None
    gap_end_col = None
    anchor_specs = (
        (left_anchor_start, left_anchor_end, left_anchor_vals),
        (right_anchor_start, right_anchor_end, right_anchor_vals),
        )  
      #walker logic, the walker can only progress along the ops and enxtend the identity vec if the the operation is a base consuming operation.
    for length, op in ops:
        col_start = len(identity_vec)
        op_ref_end = ref_pos + length if op in {"=", "X", "D", "N"} else ref_pos

        if gap_ref_start is not None and op in {"=", "X", "D", "N"}:
            if gap_start_col is None and ref_pos <= gap_ref_start < op_ref_end:
                gap_start_col = col_start + (gap_ref_start - ref_pos)

#changed from  [if gap_end_col is None and ref_pos < gap_ref_end <= op_ref_end] the is none disturbes the insertion handler later and doesnt update the logic.
            if ref_pos < gap_ref_end <= op_ref_end:
                gap_end_col = col_start + (gap_ref_end - ref_pos)

        if op in {"=", "X"}:
            value = 1 if op == "=" else 0
            identity_vec.extend([value] * length)

            if gap_ref_start is not None:
                overlap_start = max(ref_pos, gap_ref_start)
                overlap_end = min(op_ref_end, gap_ref_end)
                if overlap_start < overlap_end:
                    if query_gap_start is None:
                        query_gap_start = query_pos + (overlap_start - ref_pos)
                    query_gap_end = query_pos + (overlap_end - ref_pos)
                    covered_n += overlap_end - overlap_start

            for anchor_start, anchor_end, anchor_values in anchor_specs:
                if anchor_start is None or anchor_end is None or anchor_start >= anchor_end:
                    continue
                overlap_start = max(ref_pos, anchor_start)
                overlap_end = min(op_ref_end, anchor_end)
                overlap = max(0, overlap_end - overlap_start)
                if overlap:
                    anchor_values.extend([value] * overlap)
        #insertion handler. Insertions in gap regions are especially annoying to deal with. They don't consume reference bases, so they don't have a clear position relative to the gap region, but they can still be informative if they occur in the right place in the query sequence.
        #  The current approach is to assign them a position based on where they occur in the reference alignment and extend the gap region in the query accordingly.
        elif op == "I":
            identity_vec.extend([0] * length)
            if gap_ref_start is not None and gap_ref_start <= ref_pos < gap_ref_end:
                if query_gap_start is None:
                    query_gap_start = query_pos
                query_gap_end = query_pos + length
                if gap_start_col is None:
                    gap_start_col = col_start
                gap_end_col = len(identity_vec)
            for anchor_start, anchor_end, anchor_values in anchor_specs:
                if anchor_start is None or anchor_end is None or anchor_start >= anchor_end:
                    continue
                if anchor_start <= ref_pos < anchor_end:
                    anchor_values.extend([0] * length)
# N operations should not be here because parasail alignment doesnt produce them but safety first.
        elif op in {"D", "N"}:
            identity_vec.extend([0] * length)

        if op in {"D", "N"}:
            if gap_ref_start is not None:
                overlap_start = max(ref_pos, gap_ref_start)
                overlap_end = min(op_ref_end, gap_ref_end)
                covered_n += max(0, overlap_end - overlap_start)

            for anchor_start, anchor_end, anchor_values in anchor_specs:
                if anchor_start is None or anchor_end is None or anchor_start >= anchor_end:
                    continue
                overlap_start = max(ref_pos, anchor_start)
                overlap_end = min(op_ref_end, anchor_end)
                overlap = max(0, overlap_end - overlap_start)
                if overlap:
                    anchor_values.extend([0] * overlap)

        if op in {"=", "X", "I", "S"}:
            query_pos += length
        if op in {"=", "X", "D", "N"}:
            ref_pos += length

    if gap_ref_start is not None:
        if gap_start_col is None and ref_pos >= gap_ref_start:
            gap_start_col = len(identity_vec)
        if gap_end_col is None and ref_pos >= gap_ref_end:
            gap_end_col = len(identity_vec)

    return {
        "identity_vec": identity_vec,
        "query_gap_start": query_gap_start,
        "query_gap_end": query_gap_end,
        "left_anchor_vals": left_anchor_vals,
        "right_anchor_vals": right_anchor_vals,
        "covered_n": covered_n,
        "gap_start_col": gap_start_col,
        "gap_end_col": gap_end_col,
    }


def _pct_identity(values):
    if not values:
        return 0.0
    return (sum(values) / len(values)) * 100.0


def _compute_metrics(
    read_id,
    ref_name,
    query,
    alignment,
    gap_region,
    reference_seq,
    anchor_size=CALLER_ANCHOR_LEN,
    min_gap_overlap= CALLER_MIN_GAP_OVERLAP,
):
    cigar_string = alignment["cigar_string"]
    ref_begin = alignment["ref_begin"]
    reference_len = len(reference_seq)

    gap_ref_start, gap_ref_end, total_n_count = gap_region
    if gap_ref_start is None or gap_ref_end is None:
        left_anchor_start = None
        left_anchor_end = None
        right_anchor_start = None
        right_anchor_end = None
    else:
        left_anchor_start = max(0, gap_ref_start - anchor_size)
        left_anchor_end = gap_ref_start
        right_anchor_start = gap_ref_end
        right_anchor_end = min(reference_len, gap_ref_end + anchor_size)

    ops = _parse_cigar_validated(cigar_string)
    ref_end = ref_begin + sum(length for length, op in ops if op in {"=", "X", "D", "N"})

    len_gap_overlap_5 = 0
    len_gap_overlap_3 = 0
    if gap_ref_start is not None and gap_ref_end is not None:
        len_gap_overlap_5 = max(0, gap_ref_start - ref_begin)
        len_gap_overlap_3 = max(0, ref_end - gap_ref_end)

    walked = _walk_cigar(
        ops,
        alignment["query_begin"],
        ref_begin,
        gap_ref_start,
        gap_ref_end,
        left_anchor_start,
        left_anchor_end,
        right_anchor_start,
        right_anchor_end,
    )

    identity_vec = walked["identity_vec"]

    total_columns = len(identity_vec)
    error_columns = total_columns - sum(identity_vec)
    abs_error_ratio = round((error_columns / total_columns) if total_columns else 0.0, 6)

    q_start = walked["query_gap_start"]
    q_end = walked["query_gap_end"]
    if q_start is None or q_end is None or q_end <= q_start:
        gap_sequence = ""
    else:
        q_start = max(0, min(len(query), q_start))
        q_end = max(0, min(len(query), q_end))
        gap_sequence = query[q_start:q_end]

    gap_pct_ratio = (walked["covered_n"] / total_n_count) if total_n_count else 0.0

    gap_start_col = walked.get("gap_start_col")
    gap_end_col = walked.get("gap_end_col")

    if gap_start_col is None or gap_end_col is None:
        identity_before_gap = identity_vec
        identity_after_gap = identity_vec
        anchor_pct_ident = 0.0
    else:
        identity_before_gap = identity_vec[:gap_start_col]
        identity_after_gap = identity_vec[gap_end_col:]
        anchor_pct_ident = _pct_identity(walked["left_anchor_vals"] + walked["right_anchor_vals"])

    """
    OLD LOGIC. SLICES IDENTITY VEC BASED ON GAP COLUMS, replaced with else above 
    else:
        identity_before_gap = identity_vec[:gap_start_col]
        identity_after_gap = identity_vec[gap_end_col:]
        left_edge = identity_vec[max(0, gap_start_col - anchor_size):gap_start_col]
        right_edge = identity_vec[gap_end_col:gap_end_col + anchor_size]
        anchor_pct_ident = _pct_identity(left_edge + right_edge)
     """
    
    align_pct_ident_5 = _pct_identity(identity_before_gap)
    align_pct_ident_3 = _pct_identity(identity_after_gap)
    align_pct_ident_end_min = min(align_pct_ident_5, align_pct_ident_3)
   
    if (
        gap_ref_start is not None
        and (len_gap_overlap_5 < min_gap_overlap[0] or len_gap_overlap_3 < min_gap_overlap[1])
    ):
        return {
            "read_id": read_id,
            "ref_name": ref_name,
            "alignment_score": int(alignment["score"]),
            "_abs_error_ratio": abs_error_ratio,
            "rel_error_ratio": 1.0,  # placeholder, updated later based on best alignment
            "pct_ident": round(_pct_identity(identity_vec), 6),
            "align_pct_ident_end_min": 0.0,
            "align_pct_ident_5": round(align_pct_ident_5, 6),
            "align_pct_ident_3": round(align_pct_ident_3, 6),
            "anchor_pct_ident": 0.0,
            "gap_pct_ratio": 0.0,
            "gap_sequence": "",
            "gap_seq_match": False,
            "len_gap_overlap_5": len_gap_overlap_5,
            "len_gap_overlap_3": len_gap_overlap_3,
            "cigar_string": alignment["cigar_string"],
            "ref_begin": alignment["ref_begin"],
            "query_begin": alignment["query_begin"],
        }

    if gap_ref_start is not None and gap_ref_end is not None:
        ref_gap_seq = reference_seq[gap_ref_start:gap_ref_end]
    else:
        ref_gap_seq = ""
    #turn True if gap sequence and and gap are same length and fulfill iupac dict     
    gap_seq_match = (
        len(gap_sequence) == len(ref_gap_seq)
        and all(
            qb.upper() in IUPAC_DICT.get(rb.upper(), set())
            for qb, rb in zip(gap_sequence, ref_gap_seq)
        )
    )

    row = {
        "read_id": read_id,
        "ref_name": ref_name,
        "alignment_score": int(alignment["score"]),
        "_abs_error_ratio": abs_error_ratio,
        "rel_error_ratio": 1.0,  # placeholder, updated later based on best alignment
        "pct_ident": round(_pct_identity(identity_vec), 6),
        "align_pct_ident_end_min": round(align_pct_ident_end_min, 6),
        "align_pct_ident_5": round(align_pct_ident_5, 6),
        "align_pct_ident_3": round(align_pct_ident_3, 6),
        "anchor_pct_ident": round(anchor_pct_ident, 6),
        "gap_pct_ratio": round(gap_pct_ratio, 6),
        "gap_sequence": gap_sequence,
        "gap_seq_match": gap_seq_match,
        "len_gap_overlap_5": len_gap_overlap_5,
        "len_gap_overlap_3": len_gap_overlap_3,
        "cigar_string": alignment["cigar_string"],
        "ref_begin": alignment["ref_begin"],
        "query_begin": alignment["query_begin"],
    }
    return row


# TODO: this is a double wrapped load ref. check if  not just _load_references is enough.
def load_references(reference_fasta, include_revcomp=False):
    references = _load_references(
        reference_fasta,
        parser=fasta_parser,
        reject_duplicate_names=True,
        reject_duplicate_sequences=True,
        include_revcomp=include_revcomp,
    )
    return references

def call_randomers(
    source_bam,
    reference_fasta,
    output_dir,
    config=None,
    matrix=None,
    mode="local",
    threshold = "moderate",
    max_err_ratio=CALLER_MAX_ERR_RATIO,
    read_len_range=(None, None),
    include_revcomp=False,
    show_progress=False,
):
    if include_revcomp:
        references = load_references(reference_fasta, include_revcomp=True)
    else:
        references = load_references(reference_fasta)
    if not references:
        raise ValueError("No reference sequences were loaded from the provided FASTA file.")

    if config is None:
        config = parer.align_config()
    if matrix is None:
        matrix = parer.build_matrix(config, parer.full_alphabet)
    score_threshold = config["score_threshold"]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    gap_regions = {ref_name: find_gap_region(ref_seq) for ref_name, ref_seq in references.items()}

    writers = {}
    file_handles = {}
    output_paths = {}
    failures_path = output_dir / "failures.tsv"
    failures_fh = failures_path.open("w", newline="")
    failures_writer = csv.DictWriter(failures_fh, fieldnames=_FAILURE_FIELDS, delimiter="\t")
    failures_writer.writeheader()
    pct_ident_values = []

    output_paths["failures"] = str(failures_path)
    file_handles["failures"] = failures_fh
    try:
        for ref_name in references:
            out_path = output_dir / f"{_safe_filename(ref_name)}_call.tsv"
            out_fh = out_path.open("w", newline="")
            writer = csv.DictWriter(out_fh, fieldnames=tsv_fields, delimiter="\t")
            writer.writeheader()

            writers[ref_name] = writer
            file_handles[ref_name] = out_fh
            output_paths[ref_name] = str(out_path)

        total_records = _count_bam_records(source_bam) if show_progress else None

        with pysam.AlignmentFile(str(source_bam), "rb", check_sq=False) as bam:
            record_iter = bam.fetch(until_eof=True)
            record_iter = tqdm(
                record_iter,
                total=total_records,
                desc="Align + TSV creation",
                unit="read",
                disable=not show_progress,
            )
            for record in record_iter:
                if not record.query_sequence:
                    failures_writer.writerow(
                        _tsv_round_row(
                            _failure_row(
                            read_id=record.query_name,
                            ref_name="*",
                            failure_reason="no_sequence",
                            )
                        )
                    )
                    continue

                query = record.query_sequence.upper()
                read_id = record.query_name
                seq_len = len(query)

                if read_len_range[0] is not None and seq_len < read_len_range[0]:
                    failures_writer.writerow(
                        _tsv_round_row(
                            _failure_row(
                            read_id=read_id,
                            ref_name="*",
                            failure_reason="too_short",
                            )
                        )
                    )
                    continue

                if read_len_range[1] is not None and seq_len > read_len_range[1]:
                    failures_writer.writerow(
                        _tsv_round_row(
                            _failure_row(
                            read_id=read_id,
                            ref_name="*",
                            failure_reason="too_long",
                            )
                        )
                    )
                    continue

                read_alignments = []

                for ref_name, ref_seq in references.items():
                    alignment = parer.align(query, ref_seq, config, matrix, mode=mode)

                    if alignment.get("saturated"):
                        failures_writer.writerow(
                            _tsv_round_row(
                                _failure_row(
                                read_id=read_id,
                                ref_name=ref_name,
                                failure_reason="saturated",
                                )
                            )
                        )
                        continue

                    cigar_string = alignment.get("cigar_string")
                    if not cigar_string or cigar_string == "*":
                        failures_writer.writerow(
                            _tsv_round_row(
                                _failure_row(
                                read_id=read_id,
                                ref_name=ref_name,
                                alignment_score=alignment.get("score", "N.a"),
                                failure_reason="no_cigar",
                                )
                            )
                        )
                        continue

                    if alignment["score"] < score_threshold:
                        failures_writer.writerow(
                            _tsv_round_row(
                                _failure_row(
                                read_id=read_id,
                                ref_name=ref_name,
                                alignment_score=alignment["score"],
                                failure_reason="below_threshold",
                                )
                            )
                        )
                        continue

                    row = _compute_metrics(
                        read_id,
                        ref_name,
                        query,
                        alignment,
                        gap_regions[ref_name],
                        ref_seq,
                    )
                    
                    #empty gap sequence check
                    if row["gap_sequence"] == "" and gap_regions[ref_name][0] is not None:
                        failures_writer.writerow(
                            _tsv_round_row(
                                _failure_row(
                                read_id=read_id,
                                ref_name=ref_name,
                                alignment_score=row["alignment_score"],
                                failure_reason="empty_gap_sequence",
                                )
                            )
                        )
                        continue
                    # Metrics are stored as percentages (0-100). End/anchor thresholds are configured as ratios.
                    align_end_ident_ratio = row["align_pct_ident_end_min"] / 100.0
                    anchor_ident_ratio = row["anchor_pct_ident"] / 100.0
                    failure_reason = None
                    if align_end_ident_ratio < CALLER_MIN_ALIGN_PCT_IDENT_END_MIN:
                        failure_reason = "below_align_end_identity_threshold"
                    elif anchor_ident_ratio < CALLER_MIN_ANCHOR_PCT_IDENT:
                        failure_reason = "below_anchor_identity_threshold"

                    if failure_reason is not None:
                        failures_writer.writerow(
                            _tsv_round_row(
                                _failure_row(
                                read_id=read_id,
                                ref_name=ref_name,
                                alignment_score=row["alignment_score"],
                                align_pct_ident_5=row["align_pct_ident_5"],
                                align_pct_ident_3=row["align_pct_ident_3"],
                                align_end_ident_ratio=round(align_end_ident_ratio, 6),
                                min_align_end_ident_ratio=CALLER_MIN_ALIGN_PCT_IDENT_END_MIN,
                                pct_ident=row["pct_ident"],
                                anchor_ident_ratio=round(anchor_ident_ratio, 6),
                                min_anchor_ident_ratio=CALLER_MIN_ANCHOR_PCT_IDENT,
                                failure_reason=failure_reason,
                                )
                            )
                        )
                        continue

                    read_alignments.append((row["_abs_error_ratio"], ref_name, row))

                if not read_alignments:
                    continue
                
                #best hit selection with sort first to find best hit with min error at the start of the list.  
                read_alignments.sort(key=lambda x: x[0])
                min_err = read_alignments[0][0]

                for abs_error_ratio, ref_name, row in read_alignments:
                    if min_err == 0:
                        if abs_error_ratio > 0:
                            break
                        row["rel_error_ratio"] = 1.0
                    else:
                        ratio = abs_error_ratio / min_err
                        if ratio > max_err_ratio:
                            break
                        row["rel_error_ratio"] = round(ratio, 3)

                    row.pop("_abs_error_ratio", None)
                    writers[ref_name].writerow(_tsv_round_row(row))
                    #pct ident collection for threshold determination.
                    pct_ident_values.append(row["pct_ident"])

        for ref_name in references:
            file_handles[ref_name].close()

        if pct_ident_values:
            threshold = determine_threshold(pct_ident_values, threshold=threshold)
        else:
            threshold = 0.2

        try:
            for ref_name in references:
                tsv_path = output_paths[ref_name]
                rows_pass = []
                rows_fail = []

                with open(tsv_path, "r") as f:
                    reader = csv.DictReader(f, delimiter="\t")
                    for row in reader:
                        if float(row["pct_ident"]) < threshold:
                            rows_fail.append(row)
                        else:
                            rows_pass.append(row)

                with open(tsv_path, "w", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=tsv_fields, delimiter="\t")
                    writer.writeheader()
                    writer.writerows(rows_pass)

                for row in rows_fail:
                    failures_writer.writerow(
                        _failure_row(
                            read_id=row["read_id"],
                            ref_name=row["ref_name"],
                            alignment_score=row["alignment_score"],
                            pct_ident=row["pct_ident"],
                            min_pct_ident=threshold,
                            failure_reason="below_dynamic_pct_ident_threshold",
                        )
                    )
        finally:
            failures_fh.close()

        output_paths["threshold"] = round(threshold, 6)
        return output_paths
    finally:
        for fh in file_handles.values():
            if not fh.closed:
                fh.close()
