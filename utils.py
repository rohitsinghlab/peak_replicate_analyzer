"""
Utility functions for ChIP-seq peak analysis.

This module contains shared functions for loading data, building peak graphs,
and analyzing conserved peaks across replicates.
"""

import os
import pandas as pd
import numpy as np
from collections import defaultdict
import networkx as nx
import time
import re
from tqdm import tqdm
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def load_peak_files(directory, file_pattern=None, min_neg_log_pval=2, prefix_groups=None):
    """
    Load peak CSV files from a directory and optionally group them by filename prefix.
    
    Args:
        directory (str): Directory containing the CSV files
        file_pattern (str): Regex pattern to match files (default: r'.*\.csv')
        min_neg_log_pval (float): Minimum -log10(p-value) threshold for filtering peaks (default: 2, meaning p-value <= 0.01)
        prefix_groups (list): List of prefixes to group files by (e.g., ['XX', 'XY'])
    
    Returns:
        dict or list: If prefix_groups provided, returns dict with groups as keys.
                     Otherwise returns list of (filename, dataframe) tuples.
    """
    if file_pattern is None:
        file_pattern = r'.*\.csv'
    
    start_time = time.time()
    logger.info(f"Loading files from {directory} with min -log10(p-value) {min_neg_log_pval}...")
    
    # Find all files matching the pattern
    pattern = re.compile(file_pattern)
    all_files = [f for f in os.listdir(directory) if pattern.match(f)]
    
    logger.info(f"Found {len(all_files)} files matching pattern: {file_pattern}")
    
    loaded_data = []
    
    # Required columns for peak analysis
    required_columns = {'PeakID', 'Chrom', 'Start', 'End', 'Pvalue'}
    
    for filename in tqdm(all_files, desc="Loading files"):
        file_path = os.path.join(directory, filename)
        try:
            df = pd.read_csv(file_path)
            
            # Check if all required columns are present
            missing_columns = required_columns - set(df.columns)
            if missing_columns:
                logger.warning(f"Skipping {filename}: missing required columns {missing_columns}")
                continue
            
            # Filter by p-value if Pvalue column exists
            # Note: Pvalue column contains -log10(p-value), so higher values = more significant
            if 'Pvalue' in df.columns:
                original_len = len(df)
                df = df[df['Pvalue'] >= min_neg_log_pval]
                filtered_len = len(df)
                logger.info(f"Filtered {filename}: {original_len} -> {filtered_len} peaks "
                           f"(removed {original_len - filtered_len} peaks below -log10(p-value) {min_neg_log_pval})")
            
            loaded_data.append((filename, df))
            
        except Exception as e:
            logger.warning(f"Skipping {filename}: error reading file - {e}")
            continue
    
    # Check if we loaded any valid files
    if not loaded_data:
        logger.error(f"No valid peak files found in {directory} matching pattern {file_pattern}")
        logger.error(f"Required columns: {required_columns}")
        return [] if not prefix_groups else {group: [] for group in prefix_groups or []}
    
    # Group by prefixes if requested
    if prefix_groups:
        grouped_data = {group: [] for group in prefix_groups}
        for filename, df in loaded_data:
            for group in prefix_groups:
                if filename.startswith(group):
                    grouped_data[group].append((filename, df))
                    break
        
        elapsed = time.time() - start_time
        group_counts = {group: len(files) for group, files in grouped_data.items()}
        logger.info(f"Loaded and grouped files: {group_counts} in {elapsed:.2f} seconds")
        return grouped_data
    
    elapsed = time.time() - start_time
    logger.info(f"Loaded {len(loaded_data)} files in {elapsed:.2f} seconds")
    return loaded_data


def do_peaks_overlap(peak1, peak2, max_distance=1000):
    """
    Check if two peaks overlap or are within max_distance of each other.
    
    Args:
        peak1, peak2: Peak objects with 'Start' and 'End' attributes
        max_distance (int): Maximum distance for considering peaks as overlapping
    
    Returns:
        bool: True if peaks overlap or are within max_distance
    """
    return max(peak1['Start'], peak2['Start']) <= min(peak1['End'], peak2['End']) + max_distance


def build_peak_graph(replicate_dfs, max_distance=1000):
    """
    Build a k-partite graph where each part represents a replicate.
    Peaks from different replicates are connected if they overlap.
    
    Args:
        replicate_dfs: List of (replicate_name, dataframe) tuples
        max_distance (int): Maximum distance for considering peaks as overlapping
    
    Returns:
        networkx.Graph: Graph with peaks as nodes and overlaps as edges
    """
    start_time = time.time()
    logger.info("Building peak graph...")
    
    G = nx.Graph()
    
    # Add all peaks as nodes, labeled with their replicate
    total_peaks = sum(len(df) for _, df in replicate_dfs)
    with tqdm(total=total_peaks, desc="Adding nodes") as pbar:
        for i, (replicate_name, df) in enumerate(replicate_dfs):
            for _, peak in df.iterrows():
                node_id = f"{replicate_name}:{peak['PeakID']}"
                G.add_node(node_id, replicate=i, peak_data=peak)
                pbar.update(1)
    
    logger.info(f"Added {total_peaks} peaks as nodes in {time.time() - start_time:.2f} seconds")
    
    # Group nodes by chromosome for more efficient comparison
    chrom_to_nodes = defaultdict(list)
    for node_id, data in G.nodes(data=True):
        chrom = data['peak_data']['Chrom']
        chrom_to_nodes[chrom].append((node_id, data))
    
    # Connect overlapping peaks from different replicates
    edge_start_time = time.time()
    edge_count = 0
    total_chroms = len(chrom_to_nodes)
    
    with tqdm(total=total_chroms, desc="Processing chromosomes") as pbar_chroms:
        for chrom, nodes in chrom_to_nodes.items():
            # Calculate total comparisons for this chromosome
            total_comparisons = sum(1 for i in range(len(nodes)) 
                                   for j in range(i+1, len(nodes)) 
                                   if nodes[i][1]['replicate'] != nodes[j][1]['replicate'])
            
            # Only show inner progress bar if there are enough comparisons
            if total_comparisons > 1000:
                inner_pbar = tqdm(total=total_comparisons, desc=f"Chrom {chrom}", leave=False)
            else:
                inner_pbar = None
                
            # Compare peaks within this chromosome
            for i, (node1_id, node1_data) in enumerate(nodes):
                for j in range(i+1, len(nodes)):
                    node2_id, node2_data = nodes[j]
                    # Only connect peaks from different replicates
                    if node1_data['replicate'] != node2_data['replicate']:
                        if do_peaks_overlap(node1_data['peak_data'], node2_data['peak_data'], max_distance):
                            G.add_edge(node1_id, node2_id)
                            edge_count += 1
                        if inner_pbar:
                            inner_pbar.update(1)
            
            if inner_pbar:
                inner_pbar.close()
            pbar_chroms.update(1)
    
    elapsed = time.time() - edge_start_time
    logger.info(f"Added {edge_count} edges in {elapsed:.2f} seconds")
    logger.info(f"Total graph building time: {time.time() - start_time:.2f} seconds")
    
    return G


def find_conserved_peaks(G, min_replicates=3, total_replicates=3):
    """
    Find clusters of peaks that contain peaks from at least min_replicates.
    
    Args:
        G: NetworkX graph with peak nodes and overlap edges
        min_replicates (int): Minimum number of replicates required for conservation
        total_replicates (int): Total number of replicates (for logging)
    
    Returns:
        list: List of conserved peak clusters (connected components)
    """
    start_time = time.time()
    logger.info(f"Finding peaks conserved across at least {min_replicates} of {total_replicates} replicates...")
    
    conserved_clusters = []
    
    # Find connected components (potential conserved peaks)
    components = list(nx.connected_components(G))
    logger.info(f"Found {len(components)} connected components")
    
    for component in tqdm(components, desc="Analyzing components"):
        # Check which replicates are represented in this component
        replicate_sets = set()
        for node in component:
            replicate = G.nodes[node]['replicate']
            replicate_sets.add(replicate)
        
        # If the component contains peaks from at least min_replicates, it's conserved
        if len(replicate_sets) >= min_replicates:
            conserved_clusters.append(component)
    
    elapsed = time.time() - start_time
    logger.info(f"Found {len(conserved_clusters)} conserved clusters in {elapsed:.2f} seconds")
    
    return conserved_clusters


def analyze_conserved_peaks(clusters, G, min_replicates):
    """
    Analyze the conserved peak clusters and extract summary statistics.
    
    Args:
        clusters: List of conserved peak clusters
        G: NetworkX graph with peak data
        min_replicates (int): Minimum number of replicates (for logging)
    
    Returns:
        pandas.DataFrame: Summary statistics for each conserved peak cluster
    """
    start_time = time.time()
    logger.info(f"Analyzing {len(clusters)} conserved peak clusters...")
    
    results = []
    
    for i, cluster in enumerate(tqdm(clusters, desc="Analyzing clusters")):
        # Extract peak data from the cluster
        peaks = [G.nodes[node]['peak_data'] for node in cluster]
        
        # Calculate aggregate values
        chrom = peaks[0]['Chrom']  # All peaks in a cluster are on the same chromosome
        start = min(peak['Start'] for peak in peaks)
        end = max(peak['End'] for peak in peaks)
        
        # Handle p-values if present
        if 'Pvalue' in peaks[0]:
            avg_pvalue = np.mean([peak['Pvalue'] for peak in peaks])
            max_pvalue = max(peak['Pvalue'] for peak in peaks)
        else:
            avg_pvalue = np.nan
            max_pvalue = np.nan
        
        # Count replicates represented
        replicate_sets = set()
        for node in cluster:
            replicate = G.nodes[node]['replicate']
            replicate_sets.add(replicate)
        num_replicates = len(replicate_sets)
        
        # Get nearby genes (could be multiple)
        genes = set()
        for peak in peaks:
            if 'Symbol' in peak and pd.notna(peak['Symbol']):
                genes.add(peak['Symbol'])
        
        # Compile result
        result = {
            'cluster_id': i,
            'chrom': chrom,
            'start': start,
            'end': end,
            'length': end - start,
            'avg_pvalue': avg_pvalue,
            'max_pvalue': max_pvalue,
            'num_replicates': num_replicates,
            'peaks': [peak['PeakID'] for peak in peaks],
            'genes': list(genes)
        }
        results.append(result)
    
    elapsed = time.time() - start_time
    logger.info(f"Analyzed conserved peaks in {elapsed:.2f} seconds")
    
    return pd.DataFrame(results)


def extract_all_peaks(replicate_dfs):
    """
    Extract all peaks from a set of dataframes, regardless of conservation.
    
    Args:
        replicate_dfs: List of (filename, dataframe) tuples
    
    Returns:
        list: List of peak dictionaries with basic genomic coordinates
    """
    all_peaks = []
    for filename, df in replicate_dfs:
        for _, peak in df.iterrows():
            all_peaks.append({
                'chrom': peak['Chrom'],
                'start': peak['Start'],
                'end': peak['End'],
                'file': filename,
                'peak_id': peak['PeakID']
            })
    return all_peaks


def save_results(results_dict, output_dir, suffix=""):
    """
    Save analysis results to CSV files.
    
    Args:
        results_dict (dict): Dictionary with result names as keys and DataFrames as values
        output_dir (str): Output directory for saving files
        suffix (str): Optional suffix for filenames
    """
    logger.info(f"Saving results to {output_dir}")
    
    for name, df in results_dict.items():
        filename = f"{name}{suffix}.csv"
        filepath = os.path.join(output_dir, filename)
        df.to_csv(filepath, index=False)
        logger.info(f"Saved {name}: {len(df)} entries to {filename}")
