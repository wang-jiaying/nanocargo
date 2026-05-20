[![DOI](https://zenodo.org/badge/1186161821.svg)](https://doi.org/10.5281/zenodo.19441193)
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

If you use nanocargo in your research, please cite:

Wang J. (2026). nanocargo (Version 0.1.0). GitHub. Zenodo. https://doi.org/10.5281/zenodo.19441193

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
 
### Preparing tag sequence files
 
Tag sequences define the expected boundaries of the amplicon and are used to
verify that reads span the full amplicon region. They correspond to the IS (or
other MGE) sequences flanking the accessory region of interest.
 
**General rule:**
Both `-s` (start tag) and `-e` (end tag) sequences must be provided in
**5'→3' orientation**, exactly as they appear in the amplicon read.
 
| Tag | Position in amplicon | Sequence to provide |
|---|---|---|
| Start tag (`-s`) | 5' end of amplicon | IS/MGE sequence at the 5' boundary of the accessory region, in 5'→3' direction |
| End tag (`-e`) | 3' end of amplicon | IS/MGE sequence at the 3' boundary of the accessory region, in 5'→3' direction (reverse complement if necessary — see below) |
 
---
 
#### Example: IS1071 composite transposon (IS1071-CT)
 
IS1071-CT consists of two IS1071 copies in inverted orientation flanking an
accessory region. Outward-facing primers bind within IS1071 and amplify the
accessory region between them. The resulting amplicon is flanked on both sides
by IS1071-derived sequence:
 
```
5'──[IS1071 →]──────────────────────[← IS1071]──3'
      ↑                                    ↑
   start tag                            end tag
 (5'→3' as-is)            (5'→3' = RC of the minus-strand IS1071)
 primer →                                ← primer
```
 
- **`starttag_IS1071.fasta`** — IS1071 sequence from the primer binding site to
  the start of the accessory region, in 5'→3' orientation
- **`endtag_IS1071.fasta`** — IS1071 sequence at the 3' end of the amplicon,
  written in 5'→3' orientation as it appears in the read (i.e. reverse
  complement of the reference IS1071 sequence on the bottom strand)
> **Important:** The end tag must be given in 5'→3' orientation **as it
> appears in the amplicon read** — not as it appears in the IS reference on the
> complementary strand. If the 3'-end IS copy is on the minus strand of your
> reference, take its reverse complement before providing it to `-e`.
 
---
 
#### For other IS elements or MGEs
 
The same logic applies to any composite transposon or MGE pair. To prepare your
tag files:
 
1. Identify the two IS/MGE copies flanking your accessory region of interest.
2. Locate the outward-facing primer binding sites within each IS/MGE copy.
3. Extract the IS/MGE sequence between the primer binding site and the
   IS/accessory junction (this is the sequence that will flank every amplicon
   read).
4. **Start tag (`-s`):** use this sequence directly in 5'→3' orientation.
5. **End tag (`-e`):** provide the sequence in 5'→3' orientation **as it appears
   in the amplicon read**. If the 3'-flanking IS copy is on the reverse strand
   of your reference, take the reverse complement of the extracted sequence
   before saving it to the FASTA file.
 
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
| `-M STR` | sup | Basecalling mode: `sup` (super-accuracy) or `fast`. Selects the corresponding parameter preset |

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
| `-M STR` | sup | Basecalling mode: `sup` (super-accuracy) or `fast`. Selects the corresponding parameter preset |
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

Two built-in parameter presets are provided for common basecalling modes:

| Mode | Flag | Config file | Recommended for |
|---|---|---|---|
| Super-accuracy (default) | `-M sup` or omit `-M` | `config.yaml` | Data basecalled with `sup` models |
| Fast | `-M fast` | `config_fast.yaml` | Data basecalled with `fast` or `hac` models |

If you need to further customize parameters beyond the presets, copy the
annotated template and edit directly:

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
