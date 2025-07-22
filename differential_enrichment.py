#!/usr/bin/env python3
"""
Differential Enrichment Analysis

This script performs differential enrichment analysis between two conditions/tissues
(e.g., XX vs XY chromosomes). It finds conserved peaks within each condition and then
identifies peaks that are unique to each condition by comparing against all peaks
from the other condition.

Usage:
    python differential_enrichment.py condition1_dir condition2_dir [options]

Example:
    python differential_enrichment.py XX_samples/ XY_samples/ --condition1-name XX --condition2-name XY
"""

import os
import argparse
import time
import logging
from collections import defaultdict
from tqdm import tqdm
import pandas as pd
from utils import (
    load_peak_files,
    build_peak_graph,
    find_conserved_peaks,
    analyze_conserved_peaks,
    extract_all_peaks,
    save_results
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def find_condition_unique_peaks(condition1_results, condition2_all_peaks, 
                               condition1_name, condition2_name, max_distance=1000):
    """
    Find peaks unique to condition1 by comparing against ALL peaks from condition2.
    
    Args:
        condition1_results: DataFrame of conserved peaks from condition1
        condition2_all_peaks: List of all peak dictionaries from condition2
        condition1_name: Name of condition1 for logging
        condition2_name: Name of condition2 for logging
        max_distance: Maximum distance for considering peaks as overlapping
    
    Returns:
        pandas.DataFrame: Peaks unique to condition1
    """
    start_time = time.time()
    logger.info(f"Finding {condition1_name} unique peaks (comparing against ALL {condition2_name} peaks)...")
    
    # Group condition2 peaks by chromosome for efficient comparison
    condition2_peaks_by_chrom = defaultdict(list)
    for peak in condition2_all_peaks:
        condition2_peaks_by_chrom[peak['chrom']].append(peak)
    
    logger.info(f"Extracted {len(condition2_all_peaks)} total {condition2_name} peaks for comparison")
    
    # Find peaks unique to condition1
    unique_peaks = []
    with tqdm(total=len(condition1_results), desc=f"Finding {condition1_name} unique peaks") as pbar:
        for i, condition1_row in condition1_results.iterrows():
            unique = True
            chrom = condition1_row['chrom']
            
            # Only compare against condition2 peaks on the same chromosome
            for condition2_peak in condition2_peaks_by_chrom.get(chrom, []):
                if (max(condition1_row['start'], condition2_peak['start']) <= 
                    min(condition1_row['end'], condition2_peak['end']) + max_distance):
                    unique = False
                    break
            
            if unique:
                unique_peaks.append(condition1_row)
            pbar.update(1)
    
    elapsed = time.time() - start_time
    logger.info(f"Found {len(unique_peaks)} peaks unique to {condition1_name} in {elapsed:.2f} seconds")
    
    return pd.DataFrame(unique_peaks)


def main():
    """Main function for differential enrichment analysis."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Perform differential enrichment analysis between two conditions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare XX vs XY samples with default settings
  python differential_enrichment.py XX_data/ XY_data/ --condition1-name XX --condition2-name XY
  
  # Custom file patterns and thresholds
  python differential_enrichment.py cond1/ cond2/ \\
    --file-pattern ".*_TCF4_.*\.csv" \\
    --min-neg-log-pval 3 \\
    --max-distance 1000 \\
    --min-replicates 2
        """
    )
    
    parser.add_argument(
        'condition1_dir', 
        help='Directory containing condition 1 peak files'
    )
    parser.add_argument(
        'condition2_dir', 
        help='Directory containing condition 2 peak files'
    )
    parser.add_argument(
        '--condition1-name', 
        default='condition1',
        help='Name for condition 1 (used in output files, default: condition1)'
    )
    parser.add_argument(
        '--condition2-name', 
        default='condition2',
        help='Name for condition 2 (used in output files, default: condition2)'
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
        default='.',
        help='Output directory for results (default: current directory)'
    )
    
    args = parser.parse_args()
    
    total_start_time = time.time()
    logger.info(f"Starting differential enrichment analysis")
    logger.info(f"Condition 1 ({args.condition1_name}): {args.condition1_dir}")
    logger.info(f"Condition 2 ({args.condition2_name}): {args.condition2_dir}")
    logger.info(f"File pattern: {args.file_pattern}")
    logger.info(f"Min -log10(p-value): {args.min_neg_log_pval}")
    logger.info(f"Max distance: {args.max_distance}")
    
    # Load files for both conditions
    condition1_dfs = load_peak_files(
        args.condition1_dir, 
        file_pattern=args.file_pattern, 
        min_neg_log_pval=args.min_neg_log_pval
    )
    
    condition2_dfs = load_peak_files(
        args.condition2_dir, 
        file_pattern=args.file_pattern, 
        min_neg_log_pval=args.min_neg_log_pval
    )
    
    if not condition1_dfs or not condition2_dfs:
        logger.error("No files found in one or both directories")
        return
    
    # Set min_replicates defaults
    condition1_total = len(condition1_dfs)
    condition2_total = len(condition2_dfs)
    condition1_min = args.min_replicates if args.min_replicates else condition1_total
    condition2_min = args.min_replicates if args.min_replicates else condition2_total
    
    # Validate min_replicates
    if condition1_min < 1 or condition1_min > condition1_total:
        logger.error(f"min-replicates must be between 1 and {condition1_total} for condition 1")
        return
    if condition2_min < 1 or condition2_min > condition2_total:
        logger.error(f"min-replicates must be between 1 and {condition2_total} for condition 2")
        return
    
    logger.info(f"Condition 1 ({args.condition1_name}): {condition1_total} replicates, "
                f"requiring {condition1_min} for conservation")
    logger.info(f"Condition 2 ({args.condition2_name}): {condition2_total} replicates, "
                f"requiring {condition2_min} for conservation")
    
    # Process condition 1 peaks
    logger.info(f"Processing {args.condition1_name} peaks...")
    condition1_graph = build_peak_graph(condition1_dfs, args.max_distance)
    condition1_conserved = find_conserved_peaks(condition1_graph, condition1_min, condition1_total)
    condition1_results = analyze_conserved_peaks(condition1_conserved, condition1_graph, condition1_min)
    logger.info(f"Found {len(condition1_conserved)} conserved {args.condition1_name} peaks")
    
    # Process condition 2 peaks
    logger.info(f"Processing {args.condition2_name} peaks...")
    condition2_graph = build_peak_graph(condition2_dfs, args.max_distance)
    condition2_conserved = find_conserved_peaks(condition2_graph, condition2_min, condition2_total)
    condition2_results = analyze_conserved_peaks(condition2_conserved, condition2_graph, condition2_min)
    logger.info(f"Found {len(condition2_conserved)} conserved {args.condition2_name} peaks")
    
    # Extract all peaks for uniqueness analysis
    logger.info("Extracting all peaks for differential analysis...")
    all_condition1_peaks = extract_all_peaks(condition1_dfs)
    all_condition2_peaks = extract_all_peaks(condition2_dfs)
    
    # Find peaks unique to each condition
    condition1_unique = find_condition_unique_peaks(
        condition1_results, all_condition2_peaks, 
        args.condition1_name, args.condition2_name, args.max_distance
    )
    
    condition2_unique = find_condition_unique_peaks(
        condition2_results, all_condition1_peaks,
        args.condition2_name, args.condition1_name, args.max_distance
    )
    
    # Define output filenames
    rep_str1 = f"{condition1_min}of{condition1_total}" if condition1_min < condition1_total else "all"
    rep_str2 = f"{condition2_min}of{condition2_total}" if condition2_min < condition2_total else "all"
    
    # Prepare results dictionary
    results_dict = {
        f"{args.condition1_name}_conserved_{rep_str1}_peaks": condition1_results,
        f"{args.condition2_name}_conserved_{rep_str2}_peaks": condition2_results,
        f"{args.condition1_name}_unique_{rep_str1}_peaks": condition1_unique,
        f"{args.condition2_name}_unique_{rep_str2}_peaks": condition2_unique
    }
    
    # Save results
    save_results(results_dict, args.output_dir)
    
    # Calculate total runtime
    total_elapsed = time.time() - total_start_time
    logger.info(f"Analysis complete in {total_elapsed:.2f} seconds")
    
    # Print summary
    print("\n=== DIFFERENTIAL ENRICHMENT ANALYSIS SUMMARY ===")
    print(f"{args.condition1_name} conserved peaks (min {condition1_min} replicates): {len(condition1_results)}")
    print(f"{args.condition2_name} conserved peaks (min {condition2_min} replicates): {len(condition2_results)}")
    print(f"{args.condition1_name} unique peaks (no overlap with ANY {args.condition2_name} peak): {len(condition1_unique)}")
    print(f"{args.condition2_name} unique peaks (no overlap with ANY {args.condition1_name} peak): {len(condition2_unique)}")
    print(f"Output directory: {args.output_dir}")
    print(f"Total runtime: {total_elapsed:.2f} seconds")
    print(f"Note: Uniqueness check was performed against ALL peaks, not just conserved peaks")
    
    return {
        f'{args.condition1_name}_conserved': condition1_results,
        f'{args.condition2_name}_conserved': condition2_results,
        f'{args.condition1_name}_unique': condition1_unique,
        f'{args.condition2_name}_unique': condition2_unique
    }


if __name__ == "__main__":
    main()
