#!/usr/bin/env python3
"""
Single Tissue Analysis

This script analyzes ChIP-seq peaks from multiple replicates of a single tissue/condition
to find conserved peaks across replicates. It uses a k-partite graph approach where
peaks from different replicates are connected if they overlap within a specified distance,
then finds connected components that span the minimum required number of replicates.

Usage:
    python single_tissue_analysis.py [directory] [options]

Example:
    python single_tissue_analysis.py . --min-neg-log-pval 3 --max-distance 1000 --min-replicates 2
"""

import os
import argparse
import time
import logging
from utils import (
    load_peak_files, 
    build_peak_graph, 
    find_conserved_peaks, 
    analyze_conserved_peaks,
    save_results
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def main():
    """Main function for single tissue analysis."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Analyze ChIP-seq peaks to find conserved peaks across replicates',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze all CSV files in current directory, require all replicates
  python single_tissue_analysis.py .
  
  # Require at least 2 of 3 replicates with custom distance threshold
  python single_tissue_analysis.py . --min-replicates 2 --max-distance 1000
  
  # Filter peaks by p-value and specify custom file pattern
  python single_tissue_analysis.py data/ --min-neg-log-pval 3 --file-pattern ".*_peaks\.csv"
        """
    )
    
    parser.add_argument(
        'directory', 
        nargs='?', 
        default='.', 
        help='Directory containing the CSV files (default: current directory)'
    )
    parser.add_argument(
        '--file-pattern', 
        type=str, 
        default=r'.*\.csv', 
        help='Regex pattern to match peak files (default: .*\\.csv)'
    )
    parser.add_argument(
        '--min-neg-log-pval', 
        type=float, 
        default=2, 
        help='Minimum -log10(p-value) threshold for filtering peaks (default: 2, meaning p-value <= 0.01)'
    )
    parser.add_argument(
        '--max-distance', 
        type=int, 
        default=1000, 
        help='Maximum distance for considering peaks as overlapping (default: 1000)'
    )
    parser.add_argument(
        '--min-replicates', 
        type=int, 
        help='Minimum number of replicates required for a conserved peak (default: all replicates)'
    )
    parser.add_argument(
        '--output-dir', 
        help='Output directory for results (default: input directory)'
    )
    parser.add_argument(
        '--output-prefix', 
        default='conserved_peaks',
        help='Prefix for output files (default: conserved_peaks)'
    )
    
    args = parser.parse_args()
    
    # Set defaults
    directory = args.directory
    output_dir = args.output_dir if args.output_dir else directory
    
    total_start_time = time.time()
    logger.info(f"Starting single tissue analysis")
    logger.info(f"Input directory: {directory}")
    logger.info(f"File pattern: {args.file_pattern}")
    logger.info(f"Min -log10(p-value): {args.min_neg_log_pval}")
    logger.info(f"Max distance: {args.max_distance}")
    
    # Load files
    replicate_dfs = load_peak_files(
        directory, 
        file_pattern=args.file_pattern, 
        min_neg_log_pval=args.min_neg_log_pval
    )
    
    if not replicate_dfs:
        logger.error("No files found matching the specified pattern")
        return
    
    total_replicates = len(replicate_dfs)
    min_replicates = args.min_replicates if args.min_replicates else total_replicates
    
    # Validate min_replicates
    if min_replicates < 1 or min_replicates > total_replicates:
        logger.error(f"min-replicates must be between 1 and {total_replicates}")
        return
    
    logger.info(f"Found {total_replicates} replicate files")
    logger.info(f"Requiring conservation across at least {min_replicates} replicates")
    
    # Build peak graph
    peak_graph = build_peak_graph(replicate_dfs, args.max_distance)
    
    # Find conserved peaks
    conserved_clusters = find_conserved_peaks(peak_graph, min_replicates, total_replicates)
    
    # Analyze conserved peaks
    conserved_results = analyze_conserved_peaks(conserved_clusters, peak_graph, min_replicates)
    
    logger.info(f"Found {len(conserved_clusters)} conserved peak clusters")
    
    # Prepare output
    rep_str = f"{min_replicates}of{total_replicates}" if min_replicates < total_replicates else "all"
    suffix = f"_{rep_str}_replicates"
    
    results_dict = {
        args.output_prefix: conserved_results
    }
    
    # Save results
    save_results(results_dict, output_dir, suffix)
    
    # Calculate and log total runtime
    total_elapsed = time.time() - total_start_time
    logger.info(f"Analysis complete in {total_elapsed:.2f} seconds")
    
    # Print summary
    print("\n=== SINGLE TISSUE ANALYSIS SUMMARY ===")
    print(f"Total replicates analyzed: {total_replicates}")
    print(f"Conserved peaks (min {min_replicates} replicates): {len(conserved_results)}")
    print(f"Output directory: {output_dir}")
    print(f"Total runtime: {total_elapsed:.2f} seconds")
    
    return {
        'conserved_peaks': conserved_results,
        'total_replicates': total_replicates,
        'min_replicates': min_replicates
    }


if __name__ == "__main__":
    main()
