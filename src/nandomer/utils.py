import os
import re
import pysam

_CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=XB])")

def fasta_parser(file_path):
    """
    Parse FASTA/FASTQ files using pysam.

    Args:
        file_path (str): Path to the FASTA/FASTQ file (supports .gz)

    Yields:
        tuple: (name, sequence) pairs
    """
    with pysam.FastxFile(file_path) as fh:
        for entry in fh:
            yield entry.name, entry.sequence.upper()


def bam_parser(bam_file):
    """
    Parse BAM files and yield read sequences.

    Args:
        bam_file (str): Path to the BAM file

    Yields:
        tuple: (read_name, sequence) pairs
    """
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
        for read in bam.fetch(until_eof=True):
            if read.query_sequence:
                yield read.query_name, read.query_sequence.upper()


def read_sequences(file_path):
    """
    Unified reader that dispatches based on file extension.

    Args:
        file_path (str): Path to FASTA, FASTQ, or BAM file

    Yields:
        tuple: (name, sequence) pairs
    """
    if file_path.lower().endswith(".bam"):
        yield from bam_parser(file_path)
    else:
        yield from fasta_parser(file_path)

# added a revComp function so load_references always returns references*2 one forward and one revComped. 

def rev_comp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(complement)[::-1]

def load_references(
    reference_fasta,
    parser=None,
    reject_duplicate_names=False,
    reject_duplicate_sequences=False,
    include_revcomp=False,
):
    """Load reference sequences from FASTA/FASTQ into a name -> sequence mapping."""
    parser = fasta_parser if parser is None else parser

    references = {}
    for name, sequence in parser(str(reference_fasta)):
        if reject_duplicate_names and name in references:
            raise ValueError(f"Duplicate reference name '{name}' in {reference_fasta}")

        if reject_duplicate_sequences:
            for existing_name, existing_seq in references.items():
                if existing_seq == sequence:
                    raise ValueError(
                        f"Duplicate reference sequence: '{name}' == '{existing_name}'"
                    )

        references[name] = sequence

    if include_revcomp:
        for name, sequence in list(references.items()):
            references[f"{name}_revcomp"] = rev_comp(sequence)

    return references



def parse_cigar(cigar_string):
    if not cigar_string or cigar_string == "*":
        return []
    return [(int(length), op) for length, op in _CIGAR_PATTERN.findall(cigar_string)]