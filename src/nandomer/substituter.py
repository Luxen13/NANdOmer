import atexit
import csv
from pathlib import Path

import pysam
from tqdm import tqdm

from .utils import load_references
from .aligner import (
    build_reference_header,
    _build_mapped_record,
    _build_unmapped_record,
    _count_bam_records,
)
from .caller import (
    call_randomers,
    find_gap_region,
    _parse_cigar_validated,
)


def _build_modified_reference(reference_seq, gap_start, gap_end, gap_sequence):
    return reference_seq[:gap_start] + gap_sequence + reference_seq[gap_end:]


def _compute_md_tag(cigar_string, modified_ref, ref_begin):
    ops = _parse_cigar_validated(cigar_string)
    ref_pos = ref_begin
    md_parts = []
    match_count = 0

    for length, op in ops:
        if op == "=":
            match_count += length
            ref_pos += length
        elif op == "X":
            for i in range(length):
                md_parts.append(str(match_count))
                match_count = 0
                md_parts.append(modified_ref[ref_pos].upper())
                ref_pos += 1
        elif op == "D":
            md_parts.append(str(match_count))
            match_count = 0
            deleted_bases = modified_ref[ref_pos : ref_pos + length].upper()
            md_parts.append("^" + deleted_bases)
            ref_pos += length
        elif op in ("I", "S"):
            pass
        elif op == "N":
            ref_pos += length

    md_parts.append(str(match_count))

    return "".join(md_parts)


def _collect_passing_reads(caller_output_paths, references):
	"""Read caller TSVs and keep the best passing row per read_id."""
	passing_reads = {}
	for ref_name in references:
		if ref_name not in caller_output_paths:
			continue
		tsv_path = caller_output_paths[ref_name]
		with open(tsv_path, "r") as f:
			reader = csv.DictReader(f, delimiter="\t")
			for row in reader:
				read_id = row["read_id"]
				try:
					rel_error_ratio = float(row.get("rel_error_ratio", "inf"))
				except (TypeError, ValueError):
					rel_error_ratio = float("inf")

				candidate = {
					"ref_name": ref_name,
					"gap_sequence": row["gap_sequence"],
					"gap_seq_match": row["gap_seq_match"],
					"alignment_score": int(row["alignment_score"]),
					"cigar_string": row["cigar_string"],
					"ref_begin": int(row["ref_begin"]),
					"query_begin": int(row["query_begin"]),
					"rel_error_ratio": rel_error_ratio,
				}
				if read_id not in passing_reads:
					passing_reads[read_id] = candidate
					continue

				current = passing_reads[read_id]
				if (
					candidate["rel_error_ratio"] < current["rel_error_ratio"]
					or (
						candidate["rel_error_ratio"] == current["rel_error_ratio"]
						and candidate["alignment_score"] > current["alignment_score"]
					)
				):
					passing_reads[read_id] = candidate
	return passing_reads

def _rewrite_cigar_after_substitution(cigar_string, query, modified_ref, query_begin, ref_begin):
    """Rewrite = and X ops to reflect the modified (N-substituted) reference.
 
    After substitution, positions where parasail reported X in N positions (because N != base)
    are now matches. Walk the CIGAR position-by-position and correct each
    =/X op based on the actual query vs modified reference comparison.
    """
    ops = _parse_cigar_validated(cigar_string)
    new_ops = []
    q_pos = query_begin
    r_pos = ref_begin
 
    for length, op in ops:
        if op in ("=", "X"):
            for i in range(length):
                if query[q_pos].upper() == modified_ref[r_pos].upper():
                    new_ops.append("=")
                else:
                    new_ops.append("X")
                q_pos += 1
                r_pos += 1
        else:
            new_ops.append((length, op))
            if op in ("I", "S"):
                q_pos += length
            if op in ("D", "N"):
                r_pos += length
 
    # Merge consecutive identical ops
    merged = []
    for item in new_ops:
        if isinstance(item, str):
            op = item
            if merged and merged[-1][1] == op:
                merged[-1] = (merged[-1][0] + 1, op)
            else:
                merged.append((1, op))
        else:
            length, op = item
            if merged and merged[-1][1] == op:
                merged[-1] = (merged[-1][0] + length, op)
            else:
                merged.append((length, op))
 
    return "".join(f"{length}{op}" for length, op in merged)
 

def create_bam_with_md(
	source_bam,
	reference_fasta,
	output_dir,
	config=None,
	matrix=None,
	mode="local",
	caller_kwargs=None,
	show_progress=False,
):
	references = load_references(reference_fasta)
	if not references:
		raise ValueError("No reference sequences were loaded from the provided FASTA file.")

	if config is None or matrix is None:
		from . import parer
		if config is None:
			config = parer.align_config()
		if matrix is None:
			matrix = parer.build_matrix(config, parer.full_alphabet)

	header = build_reference_header(references)
      
	output_dir = Path(output_dir)
	output_dir.mkdir(parents=True, exist_ok=True)
	source_stem = Path(source_bam).stem
	main_bam_path = output_dir / f"{source_stem}.substitute.bam"
	#failure_bam_path = output_dir / f"{source_stem}.failures.bam"

	caller_out_dir = output_dir / "caller_results"
	caller_kw = caller_kwargs or {}
	caller_output_paths = call_randomers(
		source_bam,
		reference_fasta,
		str(caller_out_dir),
		config=config,
		matrix=matrix,
		mode=mode,
		show_progress=show_progress,
		**caller_kw,
	)

	passing_reads = _collect_passing_reads(caller_output_paths, references)

	gap_regions = {
		ref_name: find_gap_region(ref_seq)
		for ref_name, ref_seq in references.items()
	}

	total_records = _count_bam_records(source_bam) if show_progress else None

	with pysam.AlignmentFile(str(source_bam), "rb", check_sq=False) as source:
		bam_main = pysam.AlignmentFile(str(main_bam_path), "wb", header=header)
		#bam_fail = pysam.AlignmentFile(str(failure_bam_path), "wb", header=header)
		atexit.register(bam_main.close)
		#atexit.register(bam_fail.close)
		try:
			record_iter = source.fetch(until_eof=True)
			record_iter = tqdm(
				record_iter,
				total=total_records,
				desc="Substitute + BAM creation",
				unit="read",
				disable=not show_progress,
			)
			for source_record in record_iter:
				if not source_record.query_sequence:
					continue

				read_id = source_record.query_name

				if read_id not in passing_reads:
					#bam_fail.write(_build_unmapped_record(source_record, header))
					continue

				read_data = passing_reads[read_id]
				
				ref_name = read_data["ref_name"]
				ref_seq = references[ref_name]
				gap_sequence = read_data["gap_sequence"]

				alignment = {
					"ref_name": ref_name,
					"score": read_data["alignment_score"],
					"cigar_string": read_data["cigar_string"],
					"ref_begin": read_data["ref_begin"],
					"query_begin": read_data["query_begin"],
				}

				if read_data["gap_sequence"] == "" or read_data.get("gap_seq_match") == "False":
					#bam_fail.write(_build_mapped_record(source_record, header, alignment))
					continue
				
				gap_start, gap_end, _ = gap_regions[ref_name]
				if gap_start is not None and gap_end is not None and gap_sequence:
					modified_ref = _build_modified_reference(
						ref_seq, gap_start, gap_end, gap_sequence
					)
				else:
					modified_ref = ref_seq

 				#Rewrite CIGAR so N-substituted positions become = instead of X
				alignment["cigar_string"] = _rewrite_cigar_after_substitution(
					alignment["cigar_string"],
					source_record.query_sequence,
					modified_ref,
					alignment["query_begin"],
					alignment["ref_begin"],
				)
				md_tag = _compute_md_tag(
					alignment["cigar_string"],
					modified_ref,
					alignment["ref_begin"],
				)

				record = _build_mapped_record(source_record, header, alignment)
				record.set_tag("MD", md_tag)
				bam_main.write(record)
		finally:
			bam_main.close()
			#bam_fail.close()
			atexit.unregister(bam_main.close)
			#atexit.unregister(bam_fail.close)

	return str(main_bam_path)
