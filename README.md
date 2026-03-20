# nanocargo

## Overview

nanocargo is a set of shell and Python scripts for analyzing Oxford Nanopore
long-read amplicon sequencing data containing known repeated tag sequences at
both ends. It is particularly suited for mixed-length LR-PCR products.

The core module is `nanocargo_repseq`. Optional modules for preprocessing, 
tag-read filtering, and polishing are also provided but can be replaced with 
user-defined steps.

---

## Citation

If you use nanocargo in your research, please cite this repository:

Wang J. (2026). nanocargo. GitHub. https://github.com/wang-jiaying/nanocargo

---

## Installation

### 1. Dependencies

All dependencies can be installed using the provided `nanocargo_env.yml`.

### 2. Clone the repository

```bash
git clone https://github.com/wang-jiaying/nanocargo.git
cd nanocargo
```

### 3. Create the environment

Using conda:
```bash
conda env create -f nanocargo_env.yml
```

Using micromamba (recommended):
```bash
micromamba env create -f nanocargo_env.yml
```

### 4. Activate the environment

Using conda:
```bash
conda activate nanocargo
export PATH=/path/to/nanocargo/scripts:$PATH
```

Using micromamba (recommended):
```bash
micromamba activate nanocargo
export PATH=/path/to/nanocargo/scripts:$PATH
```

To deactivate, run `conda deactivate` or `micromamba deactivate`.

### 5. Locate files

You can find the nanocargo directory by running:
```bash
dirname $(dirname $(which nanocargo_pipeline))
```

| File | Location |
|---|---|
| Advanced configuration | `nanocargo/scripts/config.yaml` |
| Annotated configuration template | `nanocargo/scripts/config.example.yaml` |
| Example test data | `nanocargo/test_data/example.fastq.gz` |

---

## Quick Start

The following example uses `nanocargo_pipeline` to run all steps with default
parameters on the provided test data.


```bash
cd /path/to/nanocargo/test_data

nanocargo_pipeline \
    -i example.fastq.gz \
    -o output_test \
    -s starttag_IS1071.fasta \
    -e endtag_IS1071.fasta
```

The final output of the nanocargo pipeline is:
```
output_test/polish/consensus_racon.fa
```

Intermediate results are stored in the following directories:

| Directory | Contents |
|---|---|
| `preprocess/` | Read trimming and quality filtering results |
| `taggedseq/` | Reads containing tag sequences |
| `repseq/` | Representative sequence generation |
| `polish/` | Consensus polishing |

> **Note:** For real datasets, run each module individually with tuned parameters — see [Usage](#usage) 
> and [Configuration](#configuration).

---

## Usage

### `nanocargo_preprocess`

Trims raw Nanopore reads and applies quality and length filtering. 

```
nanocargo_preprocess -i <reads.fastq.gz> -o <output_dir> [options]
nanocargo_preprocess -d <directory> -o <output_dir> [options]
```

> **Note:** Use either `-i` or `-d`, but not both.

| Flag | Default | Description |
|---|---|---|
| `-h` | — | Show help message |
| `-i` | — | Input FASTQ file ( `.fastq.gz`) |
| `-d` | — | Input directory containing multiple FASTQ files |
| `-o` | — | Output directory |
| `-t INT` | 16 | Number of threads |
| `-m INT` | 2000 | Minimum read length to retain (bp) |
| `-q INT` | 90 | Minimum mean read quality score to retain (scaled 0–100 as defined by Filtlong; equivalent to `--min_mean_q`) |

---

### `nanocargo_tag_filter`

Filters reads based on the presence of tag sequences at both ends of the amplicon.
Tag sequences for the start and end of the amplicon should be provided according
to the experimental design. See `test_data/starttag_IS1071.fasta` and
`test_data/endtag_IS1071.fasta` for IS1071 composite transposon examples.

```
nanocargo_tag_filter -i <input_reads> -o <outdir> -s <start_tag.fa> -e <end_tag.fa> [options]
```

| Flag | Default | Description |
|---|---|---|
| `-h` | — | Show help message |
| `-i` | — | Input reads file (FASTA or FASTQ; `.fa`/`.fq`, optionally gzipped) |
| `-s` | — | FASTA file containing tag sequences marking the start of the amplicon reads |
| `-e` | — | FASTA file containing tag sequences marking the end of the amplicon reads |
| `-o` | — | Output directory |
| `-t INT` | 16 | Number of threads |
| `-r FLOAT` | 0.5 | Minimum fraction of the tag sequence that must be covered in the read |
| `-x INT` | 400 | Maximum allowed distance (bp) between the tag and the read end |

---

### `nanocargo_repseq` — Core module

Identifies representative sequences from tagged reads using a containment-based
clustering approach.

```
nanocargo_repseq -i <input.fa> -o <output_dir> [options]
```

**Required arguments:**

| Flag | Description |
|---|---|
| `-i` | Input FASTA file of reads containing tag sequences |
| `-o` | Output directory |

**Optional arguments:**

| Flag | Default | Description |
|---|---|---|
| `-h` | — | Show help message |
| `-t INT` | 16 | Number of threads |
| `-paf` | — | Input all-vs-all alignment in PAF format (ava-ont). If provided, the all-vs-all alignment step will be skipped |

**Read filtering parameters:**

| Flag | Default | Description |
|---|---|---|
| `-s FLOAT` | 0.99 | Minimum coverage fraction used to determine containment (fraction of the contained read covered by the alignment) |
| `-m FLOAT` | 0.8 | Minimum alignment identity used to determine containment relationships between reads for clustering |
| `-n INT` | 0 | Maximum number of internal alignments allowed |
| `-d INT` | 1000 | Distance threshold (bp) for defining internal alignments. An alignment is considered internal if its distance from the ends of the center read exceeds this value |

**Final-representative parameters:**

| Flag | Default | Description |
|---|---|---|
| `-p FLOAT` | 0.8 | Identity threshold used to merge representative reads into final representative reads. Higher values mean only very similar reads are merged |

---

### `nanocargo_polish`

Iteratively polishes representative sequences to improve accuracy. Only high-quality overlaps between reads and representative sequences are used for Racon polishing.

```
nanocargo_polish -i <input.fa> -f <reads.fastq> -o <output_dir> [options]
```

| Flag | Default | Description |
|---|---|---|
| `-h` | — | Show help message |
| `-i` | — | Input representative sequence FASTA file |
| `-f` | — | Input FASTQ file for Racon and optional Medaka polishing |
| `-o` | — | Output directory |
| `-t INT` | 16 | Number of threads |
| `-r INT` | 3 | Number of Racon polishing rounds |
| `-m MODEL` | none | Medaka model. Medaka polishing runs only if a model is provided |
| `--keep-intermediate` | false | Keep intermediate files generated during polishing |

---

### `nanocargo_pipeline` — Full pipeline

Runs all four modules in sequence with default parameters. See [Quick Start](#quick-start) for an example.

```
nanocargo_pipeline -i <reads.fastq.gz> -o <output_dir> -s <start_tag.fa> -e <end_tag.fa> [options]
nanocargo_pipeline -d <directory>      -o <output_dir> -s <start_tag.fa> -e <end_tag.fa> [options]
```

> **Note:** Use either `-i` or `-d`, but not both.

| Flag | Default | Description |
|---|---|---|
| `-h` | — | Show help message |
| `-i` | — | Input FASTQ file (.fastq.gz) |
| `-d` | — | Input directory containing multiple FASTQ files |
| `-s` | — | FASTA file containing tag sequences marking the start of the amplicon reads |
| `-e` | — | FASTA file containing tag sequences marking the end of the amplicon reads |
| `-o` | — | Output directory |
| `-t INT` | 16 | Number of threads |
| `-r INT` | 3 | Number of Racon polishing rounds |
| `-m MODEL` | none | Medaka model. Medaka polishing runs only if a model is provided |
| `--keep-intermediate` | false | Keep intermediate files generated during polishing |

---

## Configuration

The default parameters were optimized for Oxford Nanopore data basecalled with 
`dna_r10.4.1_e8.2_400bps_sup@v4.3.0` and may need adjustment for other 
basecalling models or dataset characteristics. If you need to tune the pipeline,
copy the annotated template to `config.yaml` and edit the values:

```bash
cp /path/to/nanocargo/scripts/config.example.yaml /path/to/nanocargo/scripts/config.yaml
```

Each parameter is documented inline in `config.example.yaml`. The most commonly adjusted parameters are summarized below.

| Module | Parameter | Default | Description |
|---|---|---|---|
| `alignment_filter` | `ID_LONG` | 0.8 | Min alignment identity for long reads |
| `alignment_filter` | `ID_SHORT` | 0.9 | Min alignment identity for short reads |
| `alignment_filter` | `SCOV_LONG` | 0.99 | Min target coverage for long-read containment |
| `alignment_filter` | `SCOV_SHORT` | 0.99 | Min target coverage for short-read containment |
| `alignment_filter` | `SHORT_LEN` | 4000 | Length cutoff (bp) separating short and long reads |
| `collapse_reads` | `SCOV_TH` | 0.99 | Min target coverage for containment (overridden by `-s`) |
| `collapse_reads` | `ID_TH` | 0.8 | Min identity for containment and endpoint classification (overridden by `-m`) |
| `collapse_reads` | `MA_TH` | 0 | Max internal alignments allowed in a center read (overridden by `-n`) |
| `collapse_reads` | `MA_DELTA` | 1000 | Distance threshold (bp) for internal alignment classification (overridden by `-d`) |
| `rep_superbin` | `ID_TH2` | 0.8 | Min identity for merging representatives into final representatives (overridden by `-p`) |
| `racon_polish` | `min_identity_long` | 0.75 | Min identity for long-read alignments passed to Racon |
| `racon_polish` | `min_identity_short` | 0.9 | Min identity for short-read alignments passed to Racon |
| `racon_polish` | `min_qcov_long` | 0.98 | Min query coverage for long-read alignments passed to Racon |
| `racon_polish` | `min_qcov_short` | 0.98 | Min query coverage for short-read alignments passed to Racon |

---

## License
This project is licensed under the [GNU General Public License v3.0](LICENSE).