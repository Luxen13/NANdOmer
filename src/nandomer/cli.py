import argparse
import re
import sys
from pathlib import Path

from .aligner import create_bam_from_alignments
from .caller import call_randomers, load_references
from .doctor import run_doctor
from .postprocess import postprocess_bam, sort_and_index_bam
from .substituter import create_bam_with_md

_SAFE_NAME_PATTERN = re.compile(r"[^A-Za-z0-9._-]+")


def _safe_ref_name(name):
    safe = _SAFE_NAME_PATTERN.sub("_", name).strip("._")
    return safe or "motif"


def _ensure_file_exists(path_value, label):
    path = Path(path_value)
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def _ensure_output_dir(path_value):
    output_dir = Path(path_value)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def _validate_read_len_range(min_len, max_len):
    if min_len is not None and min_len < 0:
        raise ValueError("--read-len-min must be >= 0")
    if max_len is not None and max_len < 0:
        raise ValueError("--read-len-max must be >= 0")
    if min_len is not None and max_len is not None and min_len > max_len:
        raise ValueError("--read-len-min cannot be greater than --read-len-max")


def _aligned_bam_path(query_bam, output_dir):
    return output_dir / f"{query_bam.stem}.aligned.bam"


def _expected_call_outputs(reference_fasta, output_dir, include_revcomp=False):
    references = load_references(reference_fasta, include_revcomp=include_revcomp)
    expected = [output_dir / "failures.tsv"]
    for ref_name in references:
        expected.append(output_dir / f"{_safe_ref_name(ref_name)}_call.tsv")
    return expected


def _handle_collision(paths, force):
    existing = [path for path in paths if path.exists()]
    if not existing:
        return

    if not force:
        joined = "\n".join(str(path) for path in existing)
        raise FileExistsError(
            "Output file(s) already exist. Re-run with --force to overwrite:\n"
            f"{joined}"
        )

    for path in existing:
        if path.is_file():
            path.unlink()

#TODO: not used anymore as align is legacy command, decide what to do 
def _build_align_parser(prog):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Align reads from a BAM to randomer references and write aligned BAM output.",
    )
    parser.add_argument("-q", "--query", required=True, help="Path to input BAM file")
    parser.add_argument("-r", "--reference", required=True, help="Path to reference FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output directory path")
    parser.add_argument(
        "--mode",
        choices=["local", "semi_global", "global"],
        default="local",
        help="Alignment mode (default: local)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output BAM file",
    )
    parser.add_argument(
        "-postp",
        "--postp",
        action="store_true",
        help="Sort and index aligned BAM using samtools postprocessing",
    )
    return parser


def _build_call_parser(prog):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Call randomer metrics from an aligned BAM and write TSV outputs.",
    )
    parser.add_argument("-q", "--query", required=True, help="Path to aligned BAM file")
    parser.add_argument("-r", "--reference", required=True, help="Path to reference FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output directory path")
    parser.add_argument(
        "--mode",
        choices=["local", "semi_global", "global"],
        default="local",
        help="Alignment mode for caller scoring (default: local)",
    )
    parser.add_argument(
        "-t", "--threshold",
        choices=["moderate", "high", "superhigh"],
        default="moderate",
        help="Dynamic pct_ident threshold for scoring (default: moderate)",
    )
    parser.add_argument(
        "--max-err-ratio",
        type=float,
        default=1.5,
        help="Maximum relative error ratio to keep rows (default: 1.5)",
    )
    parser.add_argument("--read-len-min", type=int, default=None, help="Minimum read length filter")
    parser.add_argument("--read-len-max", type=int, default=None, help="Maximum read length filter")
    parser.add_argument(
        "--revComp",
        action="store_true",
        help="Include reverse-complemented references in calling",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing caller TSV files",
    )
    return parser


def _build_full_parser(prog):
    raise NotImplementedError("The full-run command has been removed; run align and call separately.")

#TODO: legacy as postprocess is not callable anymore, decide what to do with this
def _build_postprocess_parser(prog):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Postprocess an aligned BAM with samtools sort/calmd/index.",
    )
    parser.add_argument("-q", "--query", required=True, help="Path to aligned BAM file")
    parser.add_argument("-r", "--reference", required=True, help="Path to reference FASTA")
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Threads for samtools sort/index (default: 4)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing postprocessed BAM/BAM index outputs",
    )
    return parser


def _run_align(argv):
    args = _build_align_parser("nandomer_align").parse_args(argv)
    print("start process (nandomer_align)")
    query_bam = _ensure_file_exists(args.query, "Input BAM")
    reference_fasta = _ensure_file_exists(args.reference, "Reference FASTA")
    output_dir = _ensure_output_dir(args.output)

    aligned_bam = _aligned_bam_path(query_bam, output_dir)
    output_targets = [aligned_bam]
    if args.postp:
        final_bam = aligned_bam.with_name(f"{aligned_bam.stem}_final.bam")
        output_targets.extend([final_bam, Path(f"{final_bam}.bai")])
    _handle_collision(output_targets, args.force)

    result_bam = create_bam_from_alignments(
        source_bam=str(query_bam),
        reference_fasta=str(reference_fasta),
        output_bam=str(aligned_bam),
        mode=args.mode,
        show_progress=True,
    )
    if args.postp:
        result_bam, result_bai = postprocess_bam(result_bam, str(reference_fasta))
    print(f"Aligned BAM: {result_bam}")
    if args.postp:
        print(f"BAM index: {result_bai}")


def _run_call(argv):
    args = _build_call_parser("nandomer_call").parse_args(argv)
    print("start process (nandomer_call)")
    query_bam = _ensure_file_exists(args.query, "Input BAM")
    reference_fasta = _ensure_file_exists(args.reference, "Reference FASTA")
    _validate_read_len_range(args.read_len_min, args.read_len_max)
    output_dir = _ensure_output_dir(args.output)

    expected_outputs = _expected_call_outputs(
        reference_fasta,
        output_dir,
        include_revcomp=args.revComp,
    )
    _handle_collision(expected_outputs, args.force)

    output_paths = call_randomers(
        source_bam=str(query_bam),
        reference_fasta=str(reference_fasta),
        output_dir=str(output_dir),
        mode=args.mode,
        threshold=args.threshold,
        max_err_ratio=args.max_err_ratio,
        read_len_range=(args.read_len_min, args.read_len_max),
        include_revcomp=args.revComp,
        show_progress=True,
    )

    print("Generated TSV outputs:")
    for ref_name in sorted(k for k in output_paths if k != "failures"):
        print(f"  {ref_name}: {output_paths[ref_name]}")
    print(f"  failures: {output_paths['failures']}")

#TODO: legacy as full-run is not callable anymore, decide what to do with this
def _run_full(argv):
    raise NotImplementedError("The nandomer_full runner is no longer available.")


def _build_substitute_parser(prog):
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Align reads, call gap sequences, and write BAM with MD tags against substituted reference.",
    )
    parser.add_argument("-q", "--query", required=True, help="Path to input BAM file")
    parser.add_argument("-r", "--reference", required=True, help="Path to reference FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output directory path")
    parser.add_argument(
        "--mode",
        choices=["local", "semi_global", "global"],
        default="local",
        help="Alignment mode (default: local)",
    )
    parser.add_argument(
        "-t", "--threshold",
        choices=["moderate", "high", "superhigh"],
        default="moderate",
        help="Dynamic pct_ident threshold for scoring (default: moderate)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output files",
    )
    parser.add_argument(
        "--sort",
        action="store_true",
        help="Sort and index output BAMs (preserves custom MD tags)",
    )
    return parser


def _run_substitute(argv):
    args = _build_substitute_parser("nandomer_substitute").parse_args(argv)
    print("start process (nandomer_substitute)")
    query_bam = _ensure_file_exists(args.query, "Input BAM")
    reference_fasta = _ensure_file_exists(args.reference, "Reference FASTA")
    output_dir = _ensure_output_dir(args.output)

    source_stem = query_bam.stem
    main_bam = output_dir / f"{source_stem}.substitute.bam"
    failure_bam = output_dir / f"{source_stem}.failures.bam"
    targets = [main_bam, failure_bam]
    if args.sort:
        targets.extend([
            main_bam.with_name(f"{source_stem}.substitute_sorted.bam"),
            Path(f"{main_bam.with_name(f'{source_stem}.substitute_sorted.bam')}.bai"),
            failure_bam.with_name(f"{source_stem}.failures_sorted.bam"),
            Path(f"{failure_bam.with_name(f'{source_stem}.failures_sorted.bam')}.bai"),
        ])
    _handle_collision(targets, args.force)

    result_main = create_bam_with_md(
        source_bam=str(query_bam),
        reference_fasta=str(reference_fasta),
        output_dir=str(output_dir),
        mode=args.mode,
        caller_kwargs={"threshold": args.threshold},
        show_progress=True,
    )

    if args.sort:
        result_main, main_bai = sort_and_index_bam(result_main)

    print(f"Substituted BAM output: {result_main}")
    if args.sort:
        print(f"Main BAM index output: {main_bai}")


#TODO: legacy as postprocess is not callable anymore, decide what to do with this
def _run_postprocess(argv):
    args = _build_postprocess_parser("nandomer_postprocess").parse_args(argv)
    print("start process (nandomer_postprocess)")
    query_bam = _ensure_file_exists(args.query, "Input BAM")
    reference_fasta = _ensure_file_exists(args.reference, "Reference FASTA")

    final_bam = query_bam.with_name(f"{query_bam.stem}_final.bam")
    _handle_collision([final_bam, Path(f"{final_bam}.bai")], args.force)

    result_bam, result_bai = postprocess_bam(
        str(query_bam),
        str(reference_fasta),
        threads=args.threads,
    )
    print(f"Postprocessed BAM: {result_bam}")
    print(f"BAM index: {result_bai}")


def _run_doctor(_argv):
    output, ok = run_doctor()
    print(output)
    if not ok:
        raise RuntimeError(
            "Environment check failed. Install missing dependencies to use all features."
        )


def _helper_text():
    return """nandomer command overview

Commands:
    nandomer_call        Call randomer metrics from BAM and write TSV outputs.
    nandomer_substitute  Process reads, call gap sequences, write BAM with substituted MD tags.
    nandomer_doctor      Check environment/tool availability.

Shared required flags:
    -q, --query       Path to BAM file
    -r, --reference   Path to reference FASTA
    -o, --output      Path to output directory
    -t, --threshold   Threshold for caller scoring (moderate/high/superhigh, default: moderate)

Examples:
        nandomer_call -q your_reads.bam -r your_reference.fa -o output_call
        nandomer_substitute -q your_reads.bam -r your_reference.fa -o output_sub
        nandomer_substitute -q your_reads.bam -r your_reference.fa -o output_sub --sort
        nandomer_doctor

Output naming:
    nandomer_call  -> <safe_ref_name>_call.tsv files + failures.tsv
    nandomer_substitute -> <query_stem>.substitute.bam + <query_stem>.failures.bam (--sort: _sorted.bam + .bai)

Note:
    --sort requires samtools in PATH.

Use --help on any command for full options:
    nandomer_call --help
    nandomer_substitute --help
    nandomer_doctor
"""


def _run_and_exit(run_func, argv=None):
    try:
        run_func(sys.argv[1:] if argv is None else argv)
        return 0
    except Exception as exc:  # noqa: BLE001
        print(f"Error: {exc}", file=sys.stderr)
        return 1


def main_help():
    print(_helper_text())
    return 0




def main_call():
    raise SystemExit(_run_and_exit(_run_call))


def main_full():
    raise SystemExit(_run_and_exit(_run_full))


def main_substitute():
    raise SystemExit(_run_and_exit(_run_substitute))


def main_postprocess():
    raise SystemExit(_run_and_exit(_run_postprocess))


def main_doctor():
    raise SystemExit(_run_and_exit(_run_doctor))
