import atexit
from collections.abc import Mapping
from pathlib import Path

import pysam
from tqdm import tqdm

from . import parer
from .utils import load_references
from .utils import parse_cigar



_QUERY_CONSUMING_OPS = {"M", "I", "S", "=", "X"}
_TAGS_TO_SKIP = {"AS", "NM", "MD", "SA", "XA"}

#helper for tqdm
def _count_bam_records(source_bam):
	with pysam.AlignmentFile(str(source_bam), "rb", check_sq=False) as source:
		return sum(1 for _ in source.fetch(until_eof=True))


def bases_consumed_by_cigar(cigar_string, consuming_ops=_QUERY_CONSUMING_OPS):
	return sum(length for length, op in parse_cigar(cigar_string) if op in consuming_ops)


def soft_clipped_cigar(cigar_string, query_begin, query_len):
	if not cigar_string or cigar_string == "*":
		return "*"

	query_end = query_begin + bases_consumed_by_cigar(cigar_string)
	leading_clip = max(0, query_begin)
	trailing_clip = max(0, query_len - query_end)

	cigar_parts = []
	if leading_clip > 0:
		cigar_parts.append(f"{leading_clip}S")
	cigar_parts.append(cigar_string)
	if trailing_clip > 0:
		cigar_parts.append(f"{trailing_clip}S")
	return "".join(cigar_parts)


def build_reference_header(references):
	return pysam.AlignmentHeader.from_dict(
		{
			"HD": {"VN": "1.6", "SO": "unsorted"},
			"SQ": [{"SN": ref_name, "LN": len(ref_seq)} for ref_name, ref_seq in references.items()],
		}
	)

def _build_flag(unmapped=False):
	flag = 0
	if unmapped:
		flag |= 0x4
	return flag


def _calculate_nm_from_cigar(cigar_string):
	nm = 0
	for length, op in parse_cigar(cigar_string):
		if op in {"X", "I", "D"}:
			nm += length
	return nm


def _build_passthrough_tags(source_record):
	tags = []
	for tag, value, type_code in source_record.get_tags(with_value_type=True):
		if tag in _TAGS_TO_SKIP:
			continue
		# pysam infers array ('B') tags correctly from the value object, but
		# passing value_type='B' to set_tags() raises "invalid value type 'B'".
		if type_code == "B":
			tags.append((tag, value))
		else:
			tags.append((tag, value, type_code))
	return tags


def _build_mapped_tags(source_record, alignment, ref_len):
	tags = _build_passthrough_tags(source_record)
	tags.append(("AS", int(alignment["score"])))
	tags.append(("NM", _calculate_nm_from_cigar(alignment["cigar_string"])))
	tags.append(("MN", ref_len))
	return tags


def choose_best_alignment(query, references, config, matrix, mode="local"):
	best = None
	for ref_name, ref_seq in references.items():
		alignment = parer.align(query, ref_seq, config, matrix, mode=mode)
		if best is None or alignment["score"] > best["score"]:
			best = {
				"ref_name": ref_name,
				"score": alignment["score"],
				"cigar_string": alignment["cigar_string"],
				"ref_begin": alignment["ref_begin"],
				"query_begin": alignment["query_begin"],
			}
	return best


def _build_mapped_record(source_record, header, alignment):
	query = source_record.query_sequence.upper()
	record = pysam.AlignedSegment(header)
	record.query_name = source_record.query_name
	record.query_sequence = query
	record.flag = _build_flag(unmapped=False)
	record.reference_id = header.get_tid(alignment["ref_name"])
	record.reference_start = alignment["ref_begin"]
	record.mapping_quality = 255
	record.cigarstring = soft_clipped_cigar(alignment["cigar_string"], alignment["query_begin"], len(query))
	if source_record.query_qualities is not None:
		record.query_qualities = source_record.query_qualities

	record.set_tags(
		_build_mapped_tags(
			source_record,
			alignment,
			header.get_reference_length(alignment["ref_name"]),
		)
	)
	return record


def _build_unmapped_record(source_record, header):
	record = pysam.AlignedSegment(header)
	record.query_name = source_record.query_name
	record.query_sequence = (source_record.query_sequence or "").upper()
	record.flag = _build_flag(unmapped=True)
	record.reference_id = -1
	record.reference_start = 0
	record.mapping_quality = 0
	record.cigarstring = "*"
	if source_record.query_qualities is not None:
		record.query_qualities = source_record.query_qualities

	record.set_tags(_build_passthrough_tags(source_record))
	return record


def create_bam_from_alignments(
	source_bam,
	reference_fasta,
	output_bam,
	config=None,
	matrix=None,
	mode="local",
	show_progress=False,
):
	references = load_references(reference_fasta)
	if not references:
		raise ValueError("No reference sequences were loaded from the provided FASTA file.")

	if config is None:
		config = parer.align_config()
	if matrix is None:
		matrix = parer.build_matrix(config, parer.full_alphabet)

	header = build_reference_header(references)
	streaming_output = output_bam is None
	if streaming_output:
		output_mode = "w"
		final_bam_path = "-"
	else:
		output_bam = Path(output_bam)
		output_mode = "wb"
		final_bam_path = str(output_bam)
	score_threshold = config["score_threshold"]

	total_records = _count_bam_records(source_bam) if show_progress else None

	with pysam.AlignmentFile(str(source_bam), "rb", check_sq=False) as source:
		bam_out = pysam.AlignmentFile(final_bam_path, output_mode, header=header)
		atexit.register(bam_out.close)
		try:
			record_iter = source.fetch(until_eof=True)
			record_iter = tqdm(
				record_iter,
				total=total_records,
				desc="Align + BAM creation",
				unit="read",
				disable=not show_progress,
			)
			for source_record in record_iter:
				if not source_record.query_sequence:
					continue

				query = source_record.query_sequence.upper()
				alignment = choose_best_alignment(query, references, config, matrix, mode=mode)
				if alignment is None or alignment["score"] < score_threshold:
					bam_out.write(_build_unmapped_record(source_record, header))
					continue

				bam_out.write(_build_mapped_record(source_record, header, alignment))
		finally:
			bam_out.close()
			atexit.unregister(bam_out.close)

	return str(final_bam_path)



