# DNABarcodeDiscovery

A method for discovering DNA barcodes for identifying medicinal plants by combining Pi (nucleotide diversity) and sequence clusters.


### Requirements
This project depends on numpy, plotly and pymafft. pymafft is a Rust rewrite of the MAFFT algorithm, so the Rust toolchain must be installed beforehand to compile it. Python 3.10 or newer, and older than 3.14, is required.

### Usage

The workflow has four steps.

Step 1: run a multiple sequence alignment.

```bash
uv run barcode-discovery align raw.fasta -o aligned.fasta --strategy auto -t 10
```

Step 2: compute sequence diversity. This produces an interactive dot plot that helps you choose the best intervals to fuse:

```bash
# Compute sliding-window Pi; writes window_pi.tsv and an interactive HTML dot plot
uv run barcode-discovery scan aligned.fasta -o results --pi-threshold 0.05

# Optional: show a synchronized gene annotation track below the plot
# (arrows indicate strand direction; the track zooms together with the plot)
uv run barcode-discovery scan aligned.fasta -o results \
    --gff3 ref.gff3 --reference-id NC_000000
```

The coordinates in `--gff3` belong to the reference species. Using the gap pattern of the `--reference-id` sequence in the alignment, the program maps gene coordinates to alignment columns. Only features with `type=gene` are used, and the gene name is taken from `Name`.

Step 3: generate a config file template.

```bash
# Generate the config template
uv run barcode-discovery init-config -o config.toml
```

Step 4: fill in `config.toml`, mainly the intervals you selected, and pass it to the evaluation step.

```bash
# Hover over the plot to read window start positions, fill in config.toml, then run:
uv run barcode-discovery evaluate config.toml
```

### Config file format

```toml
# barcode-discovery evaluate config
# Run: barcode-discovery evaluate config.toml

# MAFFT-aligned FASTA used in the scan step
fasta = "data/salvia_cp_26spe_align.fas"

# Window Pi table written by the scan step
window_pi = "results/window_pi.tsv"

# Species ID used to report real coordinates (must match the FASTA header ID after '>',
# or the part of the header before the first whitespace)
reference_id = "NC_020431.1"

# (Optional) GFF3 annotation of the reference genome; it must describe the same genome as reference_id.
# When set, evaluate annotates each interval as geneA(10bp)+geneA-geneB(180bp)+geneB(20bp).
# Leave empty if no annotation is needed.
gff3 = "data/S.miltiorrhiza_chloroplast.gff3"
# (Optional) seqid in the GFF3 that corresponds to the reference genome, if the file has several sequences
gff3_seqid = ""

# Output directory
out = "results"

# Sequence cluster ratio threshold (%); a candidate above this value is marked PASS
ratio_threshold = 80.0

# ---- Candidate barcodes ----
# Each [[candidate]] is one candidate barcode; write as many as you need.
# intervals: one or more intervals, each = [start of the first window, start of the last window]
#   Take the start coordinates from the "start" column of window_pi.tsv (or the "window start"
#   shown when hovering over the interactive plot).
#   The windows within an interval are fused into [start of first window, end of last window].
#   Multiple intervals are concatenated into a single barcode (overlapping or adjacent intervals are merged automatically).
[[candidate]]
name = "candidate_1"
intervals = [
    #[84, 246],
]

#[[candidate]]
#name = "candidate_2"
#intervals = [
    #[111976, 112056],
#]

```

### Output files

- `window_pi.tsv`: table of sliding-window Pi values
- `pi_sliding_window.html`: interactive plot of sliding-window Pi values
- `barcode_candidates.tsv`: table of candidate barcodes
