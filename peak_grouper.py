import os
import pandas as pd
import numpy as np
from collections import defaultdict
import networkx as nx
import time
import argparse
from tqdm import tqdm
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

def load_files(directory, min_pvalue=0):
    """Load all TCF4 CSV files from the directory that match X?_TCF4_?.csv pattern and group by XX and XY."""
    import re
    
    start_time = time.time()
    logger.info(f"Loading files from {directory} with min p-value {min_pvalue}...")
    
    xx_dfs = []
    xy_dfs = []
    
    # Use regex to match the specific pattern X?_TCF4_?.csv
    file_pattern = re.compile(r'X[XY]_TCF4_\d\.csv')
    all_files = [f for f in os.listdir(directory) if file_pattern.match(f)]
    
    logger.info(f"Found {len(all_files)} files matching the pattern X?_TCF4_?.csv")
    
    for filename in tqdm(all_files, desc="Loading files"):
        file_path = os.path.join(directory, filename)
        df = pd.read_csv(file_path)
        
        # Filter by p-value before further processing
        original_len = len(df)
        df = df[df['Pvalue'] >= min_pvalue]
        filtered_len = len(df)
        
        logger.info(f"Filtered {filename}: {original_len} -> {filtered_len} peaks (removed {original_len - filtered_len} peaks below p-value {min_pvalue})")
        
        if filename.startswith('XX'):
            xx_dfs.append((filename, df))
        elif filename.startswith('XY'):
            xy_dfs.append((filename, df))
    
    elapsed = time.time() - start_time
    logger.info(f"Loaded {len(xx_dfs)} XX files and {len(xy_dfs)} XY files in {elapsed:.2f} seconds")
    
    return xx_dfs, xy_dfs


def do_peaks_overlap(peak1, peak2, max_distance=1000):
    """Check if two peaks overlap or are within max_distance of each other."""
    # Check for overlap or proximity (assuming peaks are already from the same chromosome)
    return max(peak1['Start'], peak2['Start']) <= min(peak1['End'], peak2['End']) + max_distance

def build_peak_graph(replicate_dfs, max_distance=1000):
    """
    Build a tripartite graph where each part represents a replicate.
    Peaks from different replicates are connected if they overlap.
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
    """Find clusters of peaks that contain peaks from at least min_replicates."""
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
    """Analyze the conserved peak clusters."""
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
        avg_pvalue = np.mean([peak['Pvalue'] for peak in peaks])
        max_pvalue = max(peak['Pvalue'] for peak in peaks)
        
        # Count replicates represented
        replicate_sets = set()
        for node in cluster:
            replicate = G.nodes[node]['replicate']
            replicate_sets.add(replicate)
        num_replicates = len(replicate_sets)
        
        # Get nearby genes (could be multiple)
        genes = set()
        for peak in peaks:
            if pd.notna(peak['Symbol']):
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

def extract_all_peaks(dfs):
    """Extract all peaks from a set of dataframes, regardless of conservation."""
    all_peaks = []
    for filename, df in dfs:
        for _, peak in df.iterrows():
            all_peaks.append({
                'chrom': peak['Chrom'],
                'start': peak['Start'],
                'end': peak['End'],
                'file': filename,
                'peak_id': peak['PeakID']
            })
    return all_peaks

def find_unique_peaks(xx_results, xy_dfs, xx_dfs, xy_results, max_distance=1000):
    """
    Find peaks unique to XX or XY by comparing against ALL peaks from the other set.
    This is a more stringent definition of uniqueness, now optimized by chromosome.
    """
    start_time = time.time()
    logger.info("Finding unique peaks (comparing against ALL peaks)...")
    
    # Extract all peaks from XY and XX for comparison
    logger.info("Extracting all peaks for comparison...")
    all_xy_peaks = extract_all_peaks(xy_dfs)
    all_xx_peaks = extract_all_peaks(xx_dfs)
    
    logger.info(f"Extracted {len(all_xy_peaks)} total XY peaks and {len(all_xx_peaks)} total XX peaks for uniqueness check")
    
    # Group all peaks by chromosome for more efficient comparison
    xy_peaks_by_chrom = defaultdict(list)
    for peak in all_xy_peaks:
        xy_peaks_by_chrom[peak['chrom']].append(peak)
    
    xx_peaks_by_chrom = defaultdict(list)
    for peak in all_xx_peaks:
        xx_peaks_by_chrom[peak['chrom']].append(peak)
    
    # Find peaks unique to XX (conserved XX peaks that don't overlap with ANY XY peak)
    logger.info("Finding XX unique peaks (comparing against ALL XY peaks)...")
    xx_unique = []
    with tqdm(total=len(xx_results), desc="Finding XX unique peaks") as pbar:
        for i, xx_row in xx_results.iterrows():
            unique = True
            chrom = xx_row['chrom']
            # Only compare against XY peaks on the same chromosome
            for xy_peak in xy_peaks_by_chrom.get(chrom, []):
                if max(xx_row['start'], xy_peak['start']) <= min(xx_row['end'], xy_peak['end']) + max_distance:
                    unique = False
                    break
            if unique:
                xx_unique.append(xx_row)
            pbar.update(1)
    
    # Find peaks unique to XY (conserved XY peaks that don't overlap with ANY XX peak)
    logger.info("Finding XY unique peaks (comparing against ALL XX peaks)...")
    xy_unique = []
    with tqdm(total=len(xy_results), desc="Finding XY unique peaks") as pbar:
        for i, xy_row in xy_results.iterrows():
            unique = True
            chrom = xy_row['chrom']
            # Only compare against XX peaks on the same chromosome
            for xx_peak in xx_peaks_by_chrom.get(chrom, []):
                if max(xy_row['start'], xx_peak['start']) <= min(xy_row['end'], xx_peak['end']) + max_distance:
                    unique = False
                    break
            if unique:
                xy_unique.append(xy_row)
            pbar.update(1)
    
    elapsed = time.time() - start_time
    logger.info(f"Found {len(xx_unique)} peaks unique to XX and {len(xy_unique)} peaks unique to XY in {elapsed:.2f} seconds")
    logger.info(f"Note: Uniqueness check was performed against ALL peaks, not just conserved peaks")
    
    return pd.DataFrame(xx_unique), pd.DataFrame(xy_unique)

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process TCF4 ChIP-seq peaks to find conserved peaks across replicates')
    parser.add_argument('directory', help='Directory containing the CSV files')
    parser.add_argument('--min-pvalue', type=float, default=0, help='Minimum p-value threshold for filtering peaks')
    parser.add_argument('--max-distance', type=int, default=1000, help='Maximum distance for considering peaks as overlapping')
    parser.add_argument('--min-replicates', type=int, default=3, help='Minimum number of replicates required for a conserved peak (default: 3, all replicates)')
    parser.add_argument('--output-dir', help='Output directory for results (defaults to input directory)')
    
    args = parser.parse_args()
    directory = args.directory
    min_pvalue = args.min_pvalue
    max_distance = args.max_distance
    min_replicates = args.min_replicates
    output_dir = args.output_dir if args.output_dir else directory
    
    # Validate min_replicates
    if min_replicates < 1 or min_replicates > 3:
        logger.error("min-replicates must be between 1 and 3")
        return
    
    total_start_time = time.time()
    logger.info(f"Starting TCF4 peak analysis with min p-value {min_pvalue}, max distance {max_distance}, and min replicates {min_replicates}")
    
    # Load files
    xx_dfs, xy_dfs = load_files(directory, min_pvalue)
    
    # Process XX peaks
    logger.info("Processing XX peaks...")
    xx_graph = build_peak_graph(xx_dfs, max_distance)
    xx_conserved = find_conserved_peaks(xx_graph, min_replicates, len(xx_dfs))
    xx_results = analyze_conserved_peaks(xx_conserved, xx_graph, min_replicates)
    logger.info(f"Found {len(xx_conserved)} conserved XX peaks (min {min_replicates} replicates)")
    
    # Process XY peaks
    logger.info("Processing XY peaks...")
    xy_graph = build_peak_graph(xy_dfs, max_distance)
    xy_conserved = find_conserved_peaks(xy_graph, min_replicates, len(xy_dfs))
    xy_results = analyze_conserved_peaks(xy_conserved, xy_graph, min_replicates)
    logger.info(f"Found {len(xy_conserved)} conserved XY peaks (min {min_replicates} replicates)")
    
    # Find unique peaks - now comparing against ALL peaks, not just conserved ones
    xx_unique, xy_unique = find_unique_peaks(xx_results, xy_dfs, xx_dfs, xy_results, max_distance)
    
    # Define output filenames based on min_replicates
    rep_str = f"{min_replicates}of3" if min_replicates < 3 else "all"
    
    # Save results
    logger.info(f"Saving results to {output_dir}")
    xx_results.to_csv(os.path.join(output_dir, f'XX_conserved_{rep_str}_peaks.csv'), index=False)
    xy_results.to_csv(os.path.join(output_dir, f'XY_conserved_{rep_str}_peaks.csv'), index=False)
    xx_unique.to_csv(os.path.join(output_dir, f'XX_unique_{rep_str}_peaks.csv'), index=False)
    xy_unique.to_csv(os.path.join(output_dir, f'XY_unique_{rep_str}_peaks.csv'), index=False)
    
    total_elapsed = time.time() - total_start_time
    logger.info(f"Analysis complete in {total_elapsed:.2f} seconds")
    
    # Print summary
    print("\n=== SUMMARY ===")
    print(f"XX conserved peaks (min {min_replicates} replicates): {len(xx_results)}")
    print(f"XY conserved peaks (min {min_replicates} replicates): {len(xy_results)}")
    print(f"XX unique peaks (no overlap with ANY XY peak): {len(xx_unique)}")
    print(f"XY unique peaks (no overlap with ANY XX peak): {len(xy_unique)}")
    print(f"Total runtime: {total_elapsed:.2f} seconds")
    
    return {
        'xx_conserved': xx_results,
        'xy_conserved': xy_results,
        'xx_unique': xx_unique,
        'xy_unique': xy_unique
    }

if __name__ == "__main__":
    main()
