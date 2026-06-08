import shutil
import subprocess
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


def _require_tool(name):
    if shutil.which(name) is None:
        raise RuntimeError(f"{name} not found in PATH.")


def _run(args):
    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{' '.join(args)}\n{result.stderr.strip()}")


def sort_and_index_bam(bam_path, threads=4):
    """Sort and index a BAM file without recalculating MD tags."""
    _require_tool("samtools")

    bam = Path(bam_path)
    sorted_bam = bam.with_name(f"{bam.stem}_sorted.bam")
    threads = str(max(1, int(threads)))

    logger.info("Sorting...")
    _run(["samtools", "sort", "-@", threads, "-o", str(sorted_bam), str(bam)])

    logger.info("Indexing...")
    _run(["samtools", "index", "-@", threads, str(sorted_bam)])

    return str(sorted_bam), f"{sorted_bam}.bai"


def postprocess_bam(bam_path, reference_fasta, threads=4):
    _require_tool("samtools")

    bam = Path(bam_path)
    ref = Path(reference_fasta)
    final = bam.with_name(f"{bam.stem}_final.bam")
    threads = str(max(1, int(threads)))

    logger.info("Sorting and adding MD tags...")
    sort_proc = subprocess.Popen(
        ["samtools", "sort", "-@", threads, str(bam)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    calmd_proc = subprocess.Popen(
        ["samtools", "calmd", "-b", "-", str(ref)],
        stdin=sort_proc.stdout,
        stdout=final.open("wb"),
        stderr=subprocess.PIPE,
    )
    sort_proc.stdout.close()

    _, calmd_err = calmd_proc.communicate()
    _, sort_err = sort_proc.communicate()

    if sort_proc.returncode != 0:
        raise RuntimeError(f"samtools sort failed:\n{sort_err.decode().strip()}")
    if calmd_proc.returncode != 0:
        raise RuntimeError(f"samtools calmd failed:\n{calmd_err.decode().strip()}")

    logger.info("Indexing...")
    _run(["samtools", "index", "-@", threads, str(final)])

    return str(final), f"{final}.bai"