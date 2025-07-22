# ChIP-seq Peak Replicate Analyzer

A modular toolkit for analyzing C# Custom thresholds and file patterns
python differential_enrichment.py cond1/ cond2/ \
  --file-pattern ".*_TCF4_.*\.csv" \
  --min-neg-log-pval 3 \
  --max-distance 1000 \
  --min-replicates 2q peaks across replicates to find conserved peaks and perform differential enrichment analysis between conditions.

## Overview

This toolkit addresses the common problem in ChIP-seq analysis where you have replicated assays on the same tissue/condition and want to find agreement between them. The analysis uses a k-partite graph approach where peaks from different replicates are connected if they overlap within a specified genomic distance. The toolkit then identifies connected components that span the minimum required number of replicates.

## Features

- **Single Tissue Analysis**: Find conserved peaks across replicates of a single condition
- **Differential Enrichment**: Compare peak conservation between two conditions/tissues
- **Flexible Thresholds**: Configurable p-value filtering, distance thresholds, and minimum replicate requirements
- **Efficient Processing**: Chromosome-based optimization for large datasets
- **Comprehensive Logging**: Detailed progress tracking and performance metrics

## Installation

### Dependencies

The toolkit requires the following Python packages:

```bash
pip install pandas numpy networkx tqdm
```

### Files

The toolkit consists of three main files:

- `utils.py` - Core utility functions for peak analysis
- `single_tissue_analysis.py` - Script for analyzing peaks within a single tissue/condition
- `differential_enrichment.py` - Script for comparing peaks between two conditions

## Usage

### Single Tissue Analysis

Analyze peaks from multiple replicates of a single tissue/condition:

```bash
# Basic usage - analyze all CSV files in current directory
python single_tissue_analysis.py .

# Require at least 2 of 3 replicates with custom distance threshold
python single_tissue_analysis.py . --min-replicates 2 --max-distance 500

# Filter peaks by p-value and specify custom file pattern
python single_tissue_analysis.py data/ --min-neg-log-pval 3 --file-pattern ".*_peaks\.csv"
```

#### Arguments

- `directory` - Directory containing peak files (default: current directory)
- `--file-pattern` - Regex pattern to match files (default: `.*\.csv`)
- `--min-neg-log-pval` - Minimum -log10(p-value) threshold (default: 2, meaning p-value <= 0.01)
- `--max-distance` - Maximum distance for peak overlap (default: 1000)
- `--min-replicates` - Minimum replicates required (default: all replicates)
- `--output-dir` - Output directory (default: input directory)
- `--output-prefix` - Prefix for output files (default: conserved_peaks)

### Differential Enrichment Analysis

Compare peak conservation between two conditions (e.g., XX vs XY chromosomes):

```bash
# Basic comparison between two directories
python differential_enrichment.py XX_data/ XY_data/ --condition1-name XX --condition2-name XY

# Custom thresholds and file patterns
python differential_enrichment.py cond1/ cond2/ \
  --file-pattern ".*_TCF4_.*\.csv" \
  --min-pvalue 0.05 \
  --max-distance 500 \
  --min-replicates 2
```

#### Arguments

- `condition1_dir` - Directory containing condition 1 peak files
- `condition2_dir` - Directory containing condition 2 peak files
- `--condition1-name` - Name for condition 1 (default: condition1)
- `--condition2-name` - Name for condition 2 (default: condition2)
- `--file-pattern` - Regex pattern to match files (default: `.*\.csv`)
- `--min-neg-log-pval` - Minimum -log10(p-value) threshold (default: 2, meaning p-value <= 0.01)
- `--max-distance` - Maximum distance for peak overlap (default: 1000)
- `--min-replicates` - Minimum replicates required (default: all replicates)
- `--output-dir` - Output directory (default: current directory)

## Input Data Format

Peak files should be CSV files with the following required columns:

- `PeakID` - Unique identifier for each peak
- `Chrom` - Chromosome name
- `Start` - Peak start position
- `End` - Peak end position

Optional columns:

- `Pvalue` - Statistical significance as -log10(p-value) (for filtering)
- `Symbol` - Associated gene symbol

## Output Files

### Single Tissue Analysis

- `conserved_peaks_[N]of[M]_replicates.csv` - Conserved peaks across N of M replicates

### Differential Enrichment Analysis

- `[condition1]_conserved_[N]of[M]_peaks.csv` - Conserved peaks in condition 1
- `[condition2]_conserved_[N]of[M]_peaks.csv` - Conserved peaks in condition 2
- `[condition1]_unique_[N]of[M]_peaks.csv` - Peaks unique to condition 1
- `[condition2]_unique_[N]of[M]_peaks.csv` - Peaks unique to condition 2

### Output Columns

- `cluster_id` - Unique identifier for the peak cluster
- `chrom` - Chromosome
- `start` - Cluster start position (minimum of all peaks)
- `end` - Cluster end position (maximum of all peaks)
- `length` - Cluster length
- `avg_pvalue` - Average p-value of peaks in cluster
- `max_pvalue` - Maximum p-value of peaks in cluster
- `num_replicates` - Number of replicates represented
- `peaks` - List of original peak IDs in cluster
- `genes` - List of associated genes

## Algorithm Details

### Peak Graph Construction

1. **Node Creation**: Each peak becomes a node in a graph, labeled with its replicate
2. **Edge Creation**: Peaks from different replicates are connected if they:
   - Are on the same chromosome
   - Overlap or are within `max_distance` of each other
3. **Optimization**: Peaks are grouped by chromosome for efficient comparison

### Conservation Analysis

1. **Connected Components**: Find all connected components in the peak graph
2. **Replicate Filtering**: Keep only components that span at least `min_replicates`
3. **Cluster Analysis**: Calculate summary statistics for each conserved cluster

### Uniqueness Analysis

1. **All Peak Extraction**: Extract all peaks from condition 2 (not just conserved)
2. **Chromosome Grouping**: Group peaks by chromosome for efficient comparison
3. **Overlap Testing**: Test each conserved peak from condition 1 against all condition 2 peaks
4. **Unique Identification**: Peaks with no overlaps are considered unique

## Example Workflow

### For the Original TCF4 XX/XY Analysis

```bash
# Analyze XX and XY samples separately, then compare
python differential_enrichment.py . . \
  --condition1-name XX \
  --condition2-name XY \
  --file-pattern "X[XY]_TCF4_\d\.csv" \
  --min-neg-log-pval 2 \
  --max-distance 1000 \
  --min-replicates 3
```

This will:
1. Load XX_TCF4_*.csv and XY_TCF4_*.csv files
2. Find peaks conserved across all 3 replicates in each condition
3. Identify peaks unique to XX (no overlap with ANY XY peak)
4. Identify peaks unique to XY (no overlap with ANY XX peak)

## Performance Considerations

- **Memory Usage**: Large datasets are processed chromosome by chromosome
- **Progress Tracking**: Detailed progress bars for long-running operations
- **Logging**: Comprehensive timing and performance metrics
- **Optimization**: Inner progress bars only shown for computationally intensive chromosomes

## Legacy Compatibility

The original `peak_grouper.py` script is maintained for backward compatibility but the new modular structure is recommended for new analyses.
Analyze technical replicates of ChIP-seq peak data to identify well-conserved peaks
