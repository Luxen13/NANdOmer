# NANdOmer

Nandomere is a command-line tool for preparing direct RNA nanopore sequencing reads from randomer-containing constructs for modified-base model training.

## What is Nandomere?

Nandomere aligns reads against randomer-containing construct References, extracts the actual bases each read carried at the N positions, and writes a BAM file where the per-read reference has the called bases substituted in place.

The output is designed for downstream modified-base model training (e.g. Remora), where each read needs to be paired with a reference containing real bases, not ambiguity codes.

## How does it work?

Nandomere requires:

- Nanopore direct RNA sequencing reads in a `.bam` file
- One or more randomer-containing construct references in a `.fa` file

In a first step, Nandomere performs local alignment per read against each reference using a Smith–Waterman implementation (parasail) with a custom IUPAC-aware scoring matrix. The alignments are then scanned to extract the sequence at the randomised positions and validated through quality metrics, including percent identity, anchor identity, and a length-based gap sequence check.

In a second step, passing reads are written to an output BAM. Original tags are carried over, and the MD tag is recomputed against a modified reference where the N positions have been replaced with the called gap sequence. Additional tags (AS, NM, MN) are set from the alignment.

The result is a BAM file where each read is paired with a reference that reflects its actual bases at the randomer positions.

## Installation

First you need to clone this github repository to your desired destination.

```bash
git clone https://github.com/Luxen13/NANdOmer.git
cd NANdOmer
```

Using an isolated environment is recommended to avoid dependency conflicts.

#### Conda or mamba environment (nandomer)

```bash
conda create -n nandomer python=3.11 -y
conda activate nandomer
pip install .
```
Note: `pip install` does not install external tools like `samtools`. The `-sort` option requires `samtools` in `PATH`. See [Postprocess Requirements](#postprocess-requirements).


#### Pixi environment 

This repository also includes a `pixi.toml` for installation via Pixi  https://github.com/prefix-dev/pixi.
The installation via Pixi does not require the additional installation of samtools.

```bash
pixi install
```

To activate the pixi environment in the directory:

```bash
pixi shell
```
### Post Installation

After the installation, the command
```bash
nandomer_doctor
```
can be used to verify the setup.

## CLI Commands

- `nandomer`: helper/overview command with command descriptions and examples
- `nandomer_substitute`: alignment, randomer and gap sequence calling -> writes randomer_calling TSV and BAM with substituted MD tags per Reference
- `nandomer_call`: substep of `nandomere_substitute`. Alignment and randomere_calling -> writes randomer_calling TSV per Reference
- `nandomer_doctor`: checks environment and tool availability

## Threshold 
The -t flag controls how aggressively reads are filtered by percent identity. Three presets are available:

moderate (default): median − MAD — permissive, retains most reads
high: median — filters the lower half of the distribution
superhigh: median + MAD — aggressive, keeps only the highest-quality alignments

```bash
nandomer_substitute -q your_reads.bam -r your_reference.fa -o output_sub -t superhigh
```
## Examples

```bash
nandomer_call -q output_align/your_reads.aligned.bam -r your_reference.fa -o output_call
nandomer_substitute -q your_reads.bam -r your_reference.fa -o output_sub 
nandomer_doctor
```

## Output Naming Policy

- `nandomer_align` writes `<query_stem>.aligned.bam`
- `nandomer_align -postp` also writes postprocessed outputs: `<query_stem>.aligned_final.bam` and `<query_stem>.aligned_final.bam.bai`
- `nandomer_call` writes `<safe_ref_name>_call.tsv` and `failures.tsv`
- `nandomer_substitute` writes `<query_stem>.substitute.bam` and `<query_stem>.failures.bam` in the output directory. With `--sort` it also writes `_sorted.bam` and corresponding `.bai` files.

By default, commands fail if target output files already exist. Use `--force` to overwrite outputs.

## Postprocess Requirements

`-sort` runs samtools sorting after Bam creation. Ensure `samtools` is installed and available in `PATH`.

Common install option via conda:
- conda/mamba: `conda install -c bioconda samtools`

It is recommended to install samtools in the virtual environment of NANdOmer so the postprocess works.
Global installation instances of samtools can not be recognized.

Verify setup with:

```bash
nandomer_doctor
```
