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

| Package    | Purpose                                   |
|------------|-------------------------------------------|
| `pod5`     | Reading and writing POD5 files            |
| `biopython`| Parsing FASTQ files                       |
| `pysam`    | Reading BAM/SAM files (Linux/macOS only)  |
| `typer`    | Command-line interface (CLI)              |

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
pip install git+https://github.com/ncechova/pod5_demux.git
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

| Parameter | Alias | Description | Default Value |
|---|---|---|---|
| `--seq` | `-s` | Path to the mapping file or directory (BAM/SAM/FASTQ) | — (required) |
| `--pod5` | `-p` | Path to the POD5 file or directory | — (required) |
| `--output` | `-o` | Output directory | `pod5_demux_output` |
| `--bc` | `-b` | Forces the assignment of a specific barcode (e.g., `barcode05`) to all mapped reads | — (optional) |
| `--mode` | `-m` | Output mode: `folder` or `single_file` | `folder` |
| `--threads` | `-t` | Number of CPU threads | 8 |

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

**`single_file`** — all reads of a single barcode are merged into one file:

```
output/
├── barcode01.pod5
├── barcode02.pod5
└── unclassified.pod5
```

### Filtering and Enforcing a Specific Barcode (`--bc`)

If the `--bc <barcode_name>` parameter is used, the tool will not search for barcodes in the metadata of the mapping files. Instead, it **assigns the specified barcode to all reads in the map**. At the same time, only this specific sample is filtered out - all other reads in the original POD5 files are skipped.

This feature is ideal in situations where:
- The input file (e.g., FASTQ) already contains data from only one sample but lacks barcode metadata.
- You want to efficiently and quickly extract only the data belonging to one specific sample from the original large POD5 dataset without processing the rest.

### Usage Examples

**Demultiplexing with a BAM map, output to folders:**

```bash
pod5_demux \
  -s /data/basecalled/calls.bam \
  -p /data/raw/ \
  -o /data/demuxed/
```

**Demultiplexing with a FASTQ map, everything into a single file per barcode:**

```bash
pod5_demux \
  -s /data/basecalled/fastq/ \
  -p /data/raw/ \
  -o /data/demuxed/ \
  -m single_file
```

**Change the number of CPU cores:**

```bash
pod5_demux \
  -s /data/basecalled/calls.bam \
  -p /data/raw/ \
  -o /data/demuxed/ \
  -t 16
```

**Extracting and assigning one specific barcode:**

```bash
pod5_demux \
  -s /data/basecalled/sample_barcode05.fastq \
  -p /data/raw/ \
  -o /data/demuxed_single/ \
  -b barcode05
```

---

## Output Description

Upon completion, the tool creates the following in the output directory:

- **POD5 files** divided by barcode (structured according to the selected mode)
- **`demux_stats.txt`** - a summary report containing:
  - total number of processed reads
  - total processing time
  - output mode used
  - number of unique barcodes
  - number of output files created

Example content of `demux_stats.txt`:

```
Processed reads: 150432
Total time: 47.83 s
Mode: folder
Unique barcodes: 12
Created files: 24
```

---

## Example Dataset

*(To be added)*

---

## Limitations

- **Windows:** The `pysam` library is not available on Windows. The tool automatically detects this and switches to custom fallback parsers for both BAM and SAM (which are slower than pysam). The FASTQ format is fully supported on all platforms.
- **Undemultiplexed reads:** Reads whose UUID is not found in the mapping file are placed into the `unclassified` group.
- **Input format:** The tool automatically detects the format of the input mapping file. The input folder must contain files of only one format (BAM and FASTQ cannot be combined in the same folder).