# pod5-demux

A tool for demultiplexing raw nanopore signals from POD5 files based on classified sequencing data (BAM/SAM/FASTQ).

---

## Table of Contents

1. [Tool Description](#tool-description)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Usage](#usage)
5. [Output Description](#output-description)
6. [Example Dataset](#example-dataset)
7. [Limitations](#limitations)

---

## Tool Description

`pod5-demux` splits POD5 files according to barcodes identified during basecalling. The inputs are:

- **classified sequencing data** in BAM, SAM, or FASTQ format
- **raw signal data** in POD5 format

The tool creates a set of output POD5 files, where each file (or folder) contains only the reads corresponding to a single sample (barcode).

Processing runs in parallel on multiple CPU cores using a **split–merge** strategy:
1. Each POD5 file is split into temporary parts by barcode.
2. The parts are merged into the final output files.

---

## Requirements

- Python >= 3.8
- Operating System: **Linux** (recommended), Windows (with limitations — see [Limitations](#limitations))

Dependencies are installed automatically during package installation:

| Package     | Purpose                                   |
|-------------|-------------------------------------------|
| `pod5`      | Reading and writing POD5 files            |
| `biopython` | Parsing FASTQ files                       |
| `pysam`     | Reading BAM/SAM files (Linux/macOS only)  |
| `typer`     | Command-line interface (CLI)              |

---

## Installation

### Installation from local source code

```bash
# 1. Clone or extract the source code into a directory
cd pod5-demux

# 2. Install the package (including dependencies)
pip install .
```

### Installation from GitHub

```bash
pip install git+https://github.com/BioSys-BUT/pod5_demux.git
```

After installation, the `pod5_demux` command will be available:

```bash
pod5_demux --help
```

---

## Usage

### Basic Syntax

```bash
pod5_demux -s <path_to_map> -p <path_to_pod5> -o <output_directory> [options]
```

### Parameters

| Parameter | Alias | Description | Default |
|---|---|---|---|
| `--seq` | `-s` | Path to the mapping file or directory (BAM/SAM/FASTQ) | — (required) |
| `--pod5` | `-p` | Path to the POD5 file or directory | — (required) |
| `--output` | `-o` | Output directory | `pod5_demux_output` |
| `--bc` | `-b` | Filter output to a specific barcode; also used as fallback for unannotated files (see below) | — (optional) |
| `--mode` | `-m` | Output mode: `folder` or `single_file` | `folder` |
| `--threads` | `-t` | Number of CPU threads | `8` |
| `--overwrite` | — | Overwrite the output directory if it already exists | `False` |

### Output Modes

**`folder`** (default) — preserves the original file structure, each barcode gets its own subfolder:

```
output/
├── barcode01/
│   ├── run1.pod5
│   └── run2.pod5
├── barcode02/
│   └── run1.pod5
└── unclassified/
│   ├── run1.pod5
│   └── unclassified.pod5
```

**`single_file`** — all reads belonging to a single barcode are merged into one file:

```
output/
├── barcode01.pod5
├── barcode02.pod5
└── unclassified.pod5
```

### Barcode Detection and the `--bc` Parameter

The tool determines the barcode for each read using the following priority order:

1. **Annotation inside the file** — BC/RG tag in BAM/SAM, or `barcode=barcodeXX` in the FASTQ header.
2. **Folder or file name** — if the input path (or its parent folder) is named `barcodeXX`, that name is used automatically. This covers the common case where a basecaller outputs one FASTQ file per barcode into a folder named after the barcode.
3. **`--bc` fallback** — if neither of the above yields a barcode (e.g. a plain FASTQ without annotations and with a non-standard filename), the value supplied via `--bc` is used as the barcode name for all reads in that file.

The `--bc` parameter also acts as an **output filter**: when specified, only reads assigned to that barcode are written to the output - all others are skipped. This is useful for quickly extracting a single sample from a large dataset without demultiplexing everything.

### Usage Examples

**Standard demultiplexing from a BAM map, output to per-barcode folders:**

```bash
pod5_demux \
  -s /data/basecalled/calls.bam \
  -p /data/raw/ \
  -o /data/demuxed/
```

**Demultiplexing from FASTQ files organised into per-barcode folders (auto-detection):**

```bash
# The tool detects barcodes from the folder names automatically — no --bc needed
pod5_demux \
  -s /data/basecalled/fastq/ \
  -p /data/raw/ \
  -o /data/demuxed/
```

Assuming the following input structure:
```
fastq/
├── barcode01/
│   └── reads.fastq
├── barcode02/
│   └── reads.fastq
└── barcode03/
    └── reads.fastq
```

**Demultiplexing from FASTQ into single files per barcode:**

```bash
pod5_demux \
  -s /data/basecalled/fastq/ \
  -p /data/raw/ \
  -o /data/demuxed/ \
  -m single_file
```

**Extracting only one specific barcode from a fully annotated BAM:**

```bash
pod5_demux \
  -s /data/basecalled/calls.bam \
  -p /data/raw/ \
  -o /data/demuxed_bc05/ \
  -b barcode05
```

**Assigning a barcode to a plain FASTQ file with no annotations and a non-standard name:**

```bash
# The file sample_A.fastq has no barcode in its header and is not named barcodeXX,
# so --bc is required to tell the tool which barcode to assign.
pod5_demux \
  -s /data/basecalled/sample_A.fastq \
  -p /data/raw/ \
  -o /data/demuxed_bc05/ \
  -b barcode05
```

**Using more CPU cores:**

```bash
pod5_demux \
  -s /data/basecalled/calls.bam \
  -p /data/raw/ \
  -o /data/demuxed/ \
  -t 16
```

**Overwriting an existing output directory:**

```bash
pod5_demux \
  -s /data/basecalled/calls.bam \
  -p /data/raw/ \
  -o /data/demuxed/ \
  --overwrite
```

---

## Output Description

Upon completion, the tool creates the following in the output directory:

- **POD5 files** divided by barcode (structured according to the selected mode)
- **`demux_stats.txt`** — a summary report containing:
  - input paths and run parameters
  - total number of processed reads
  - number of demultiplexed and unclassified reads with percentages
  - per-barcode read counts with recovery rates
  - map loading, split, merge, and total processing times

Example content of `demux_stats.txt`:

```
**************************************************
 DEMULTIPLEXING SUMMARY
**************************************************
Input map:               example_data/map/
Input POD5:              example_data/pod5
Output directory:        pod5_demux_output
Output mode:             folder
CPU cores:               8
--------------------------------------------------
Reads loaded from map:   380
Processed reads (POD5):  398
Demultiplexed:           379 (95.2 %)
Unclassified:            19 (4.8 %)
Unique barcodes in map:  3
Files created:           4
--------------------------------------------------
Map loading time:        0.10 s
Split time:              0.26 s
Merge time:              0.67 s
Total processing time:   1.03 s
--------------------------------------------------
Reads per barcode:
 barcode01:       127 (out of 127 from map -> 100.0 %)
 barcode02:       126 (out of 127 from map ->  99.2 %)
 barcode03:       126 (out of 126 from map -> 100.0 %)
 unclassified:    19
==================================================
```

---

## Example Dataset

A small example dataset is included in the `example_data/` directory in the repository root. It contains reads from a real Promethion run, split across 3 barcodes and 2 POD5 files, and can be used to verify that the tool is working correctly after installation.

### Structure

```
example_data/
├── pod5/
│   ├── run1.pod5
│   └── run2.pod5
├── map/
│   ├── reads.bam
│   ├── reads.sam
│   └── fastq/
│       ├── barcode01/
│       ├── barcode02/
│       └── barcode03/
└── expected_output/
    └── demux_stats.txt
```

### Running the example

```bash
# Using BAM map
pod5_demux \
  -s example_data/map/reads.bam \
  -p example_data/pod5/ \
  -o example_data/test_output_bam

# Using SAM map
pod5_demux \
  -s example_data/map/reads.sam \
  -p example_data/pod5/ \
  -o example_data/test_output_sam

# Using FASTQ map
pod5_demux \
  -s example_data/map/fastq/ \
  -p example_data/pod5/ \
  -o example_data/test_output_fastq
```

After running, compare the read counts in the output `demux_stats.txt` with `example_data/expected_output/demux_stats.txt` to verify the result.

---

## Limitations

- **Windows:** The `pysam` library is not available on Windows. The tool automatically detects this and switches to custom fallback parsers for both BAM and SAM (which are slower than pysam). The FASTQ format is fully supported on all platforms.
- **Undemultiplexed reads:** Reads whose UUID is not found in the mapping file are placed into the `unclassified` group.
- **Input format:** The tool automatically detects the format of the input mapping file. The input folder must contain files of only one format (BAM and FASTQ cannot be combined in the same folder).
- **Output directory:** If the specified output directory already exists and is not empty, the tool will exit with an error. Use `--overwrite` to clear it automatically.