#!/usr/bin/env python3
"""
Single Target MAGMA Analysis Script
Purpose: Run MAGMA analysis for all traits within a specific target
Author: Chief Scientist
Date: 2024
"""

import os
import sys
import subprocess
import csv
import logging
import re
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
import argparse

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(f'magma_single_target_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def load_trait_n_values(csv_file_path: str) -> dict:
    """
    Load trait N values from CSV file.
    
    Args:
        csv_file_path: Path to traits_n.csv file
        
    Returns:
        Dictionary mapping trait names to their N values
    """
    trait_n_mapping = {}
    
    try:
        with open(csv_file_path, 'r', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                trait_name = row.get('Traits', '').strip()
                n_value = row.get('N', '').strip()
                
                if trait_name and n_value:
                    try:
                        trait_n_mapping[trait_name] = int(n_value)
                    except ValueError:
                        logger.warning(f"Invalid N value for trait '{trait_name}': {n_value}")
                        
        logger.info(f"Loaded N values for {len(trait_n_mapping)} traits")
        return trait_n_mapping
        
    except Exception as e:
        logger.error(f"Failed to load trait N values: {e}")
        raise

def extract_phenotype_from_filename(filename: str) -> str:
    """
    Extract phenotype name from input filename.
    
    Args:
        filename: Input filename (e.g., 'st.overlap.4magma.txt')
        
    Returns:
        Phenotype name extracted from filename (e.g., 'st')
    """
    # Remove the suffix '.overlap.4magma.txt'
    phenotype = filename.replace('.overlap.4magma.txt', '')
    return phenotype

def find_trait_files(target_dir: str, pattern: str = "*.overlap.4magma.txt") -> list:
    """
    Find all MAGMA input files for a specific target.
    
    Args:
        target_dir: Target directory path
        pattern: File pattern to match
        
    Returns:
        List of tuples: (phenotype_name, full_file_path)
    """
    trait_files = []
    target_path = Path(target_dir)
    
    # Check if directory exists
    if not target_path.exists():
        logger.error(f"Target directory does not exist: {target_dir}")
        return []
    
    # Look for files in 01.magma_input subdirectory
    input_dir = target_path / "01.magma_input"
    if not input_dir.exists():
        logger.error(f"Input directory does not exist: {input_dir}")
        return []
    
    # Find all matching files
    for trait_file in input_dir.glob(pattern):
        if trait_file.is_file():
            # Extract phenotype name from filename
            phenotype_name = extract_phenotype_from_filename(trait_file.name)
            trait_files.append((phenotype_name, str(trait_file)))
    
    logger.info(f"Found {len(trait_files)} trait files in target directory")
    return trait_files

def run_magma_for_trait(target_name: str,
                       phenotype_name: str,
                       magma_file_path: str, 
                       n_value: int, 
                       output_dir: Path) -> bool:
    """
    Run MAGMA analysis for a single trait within a target.
    
    Args:
        target_name: Target directory name (e.g., adhd2017)
        phenotype_name: Phenotype name (e.g., st)
        magma_file_path: Path to the MAGMA input file
        n_value: Sample size for the trait
        output_dir: Output directory path
        
    Returns:
        True if successful, False otherwise
    """
    # Construct paths
    magma_executable = os.environ.get("FDRREG_MAGMA", "magma")
    bfile_path = os.environ.get("FDRREG_MAGMA_BFILE", "")
    gene_annot_path = os.environ.get("FDRREG_MAGMA_ANNOTATION", "")
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Output file path
    output_file = output_dir / f"{phenotype_name}"
    
    # Construct MAGMA command
    cmd = [
        magma_executable,
        "--bfile", bfile_path,
        "--pval", str(magma_file_path),
        f"use=snpid,p.decor",
        f"N={n_value}",
        "--gene-annot", gene_annot_path,
        "--out", str(output_file)
    ]
    
    # Log the command
    cmd_str = " ".join(cmd)
    logger.info(f"Running MAGMA for {target_name}/{phenotype_name}: {cmd_str}")
    
    try:
        # Run the command
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Log successful completion
        logger.info(f"Completed MAGMA analysis for {target_name}/{phenotype_name}")
        if result.stdout:
            logger.debug(f"STDOUT: {result.stdout[:500]}...")
            
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"MAGMA failed for {target_name}/{phenotype_name}: {e}")
        logger.error(f"STDERR: {e.stderr[:500] if e.stderr else 'No stderr'}")
        return False
        
    except Exception as e:
        logger.error(f"Unexpected error for {target_name}/{phenotype_name}: {e}")
        return False

def main():
    """
    Main function to run MAGMA analysis for a single target.
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Single Target MAGMA Analysis')
    parser.add_argument('--target', 
                       required=True,
                       help='Target directory name (e.g., adhd2017)')
    parser.add_argument('--base-dir', 
                       default=os.environ.get("FDRREG_RESULTS_DIR", "data/pipeline"),
                       help='Base directory containing target directories')
    parser.add_argument('--n-file',
                       default=os.environ.get("FDRREG_TRAIT_N", ""),
                       help='Path to traits_n.csv file')
    parser.add_argument('--max-workers', 
                       type=int, 
                       default=4,
                       help='Maximum number of parallel workers (default: 4)')
    parser.add_argument('--dry-run', 
                       action='store_true',
                       help='Print commands without executing')
    parser.add_argument('--log-level',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       default='INFO',
                       help='Set logging level')
    parser.add_argument('--skip-missing-n',
                       action='store_true',
                       help='Skip traits with missing N values instead of stopping')
    parser.add_argument('--output-dir',
                       default=None,
                       help='Custom output directory (default: base_dir/target/04.magma_output)')
    
    args = parser.parse_args()
    
    # Update logging level
    logger.setLevel(getattr(logging, args.log_level))
    
    # Start time
    start_time = datetime.now()
    logger.info(f"Starting single target MAGMA analysis for: {args.target}")
    logger.info(f"Start time: {start_time}")
    
    try:
        # Step 1: Load trait N values
        trait_n_mapping = load_trait_n_values(args.n_file)
        
        # Step 2: Define target path
        target_dir = Path(args.base_dir) / args.target
        
        # Step 3: Find all trait files in the target
        trait_files = find_trait_files(target_dir)
        
        if not trait_files:
            logger.warning(f"No trait files found in target directory: {target_dir}")
            return
        
        # Step 4: Process each trait
        successful = 0
        failed = 0
        skipped = 0
        
        # Prepare tasks for parallel execution
        tasks = []
        for phenotype_name, magma_file_path in trait_files:
            # Try to match phenotype name to N value
            n_value = None
            
            # Direct match
            if phenotype_name in trait_n_mapping:
                n_value = trait_n_mapping[phenotype_name]
                logger.debug(f"Direct match: {phenotype_name} -> N={n_value}")
            else:
                # Try alternative matches
                # Remove year suffixes like 2016, 2017, etc.
                phenotype_base = re.sub(r'\d{4}$', '', phenotype_name)
                if phenotype_base in trait_n_mapping:
                    n_value = trait_n_mapping[phenotype_base]
                    logger.debug(f"Base match: {phenotype_name} -> {phenotype_base} -> N={n_value}")
                else:
                    # Try partial matches
                    for csv_trait in trait_n_mapping.keys():
                        if csv_trait in phenotype_name or phenotype_name in csv_trait:
                            n_value = trait_n_mapping[csv_trait]
                            logger.debug(f"Partial match: {phenotype_name} -> {csv_trait} -> N={n_value}")
                            break
            
            if n_value is None:
                logger.warning(f"No N value found for phenotype: {phenotype_name} (file: {magma_file_path})")
                skipped += 1
                if args.skip_missing_n:
                    continue
                else:
                    logger.error(f"Stopping execution due to missing N value for {phenotype_name}")
                    logger.error(f"Add this trait to traits_n.csv or use --skip-missing-n to skip")
                    return
                
            if args.dry_run:
                # Just print the command
                if args.output_dir:
                    output_dir = Path(args.output_dir)
                else:
                    output_dir = target_dir / "04.magma_output"
                
                output_file = output_dir / f"{phenotype_name}"
                cmd_str = f"{os.environ.get('FDRREG_MAGMA', 'magma')} --bfile {os.environ.get('FDRREG_MAGMA_BFILE', '')} --pval {magma_file_path} use=snpid,p.decor N={n_value} --gene-annot {os.environ.get('FDRREG_MAGMA_ANNOTATION', '')} --out {output_file}"
                print(f"DRY RUN: {cmd_str}")
                continue
                
            # Set output directory
            if args.output_dir:
                output_dir = Path(args.output_dir)
            else:
                output_dir = target_dir / "04.magma_output"
                
            tasks.append((args.target, phenotype_name, magma_file_path, n_value, output_dir))
        
        if args.dry_run:
            logger.info(f"Dry run complete. Would process {len(tasks)} traits for target: {args.target}")
            return
        
        # Execute tasks in parallel
        if tasks:
            logger.info(f"Processing {len(tasks)} traits for target: {args.target} with {args.max_workers} workers")
            
            with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
                # Submit all tasks
                future_to_trait = {}
                for target_name, phenotype_name, magma_file_path, n_value, output_dir in tasks:
                    future = executor.submit(
                        run_magma_for_trait,
                        target_name,
                        phenotype_name,
                        magma_file_path,
                        n_value,
                        output_dir
                    )
                    future_to_trait[future] = f"{target_name}/{phenotype_name}"
                
                # Process completed tasks
                for future in as_completed(future_to_trait):
                    trait_info = future_to_trait[future]
                    try:
                        success = future.result()
                        if success:
                            successful += 1
                            logger.info(f"Successfully completed: {trait_info}")
                        else:
                            failed += 1
                            logger.error(f"Failed: {trait_info}")
                    except Exception as e:
                        logger.error(f"Exception for {trait_info}: {e}")
                        failed += 1
        
        # Summary
        end_time = datetime.now()
        duration = end_time - start_time
        
        logger.info(f"\n{'='*60}")
        logger.info(f"SINGLE TARGET ANALYSIS COMPLETE")
        logger.info(f"{'='*60}")
        logger.info(f"Target: {args.target}")
        logger.info(f"Total trait files found: {len(trait_files)}")
        logger.info(f"Successfully processed: {successful}")
        logger.info(f"Failed: {failed}")
        logger.info(f"Skipped (missing N values): {skipped}")
        logger.info(f"Total duration: {duration}")
        logger.info(f"Log file: magma_single_target_*.log")
        logger.info(f"{'='*60}\n")
        
        # Show failed tasks summary if any
        if failed > 0:
            logger.info("Some tasks failed. Check the log file for details.")
            logger.info("Failed tasks may need manual attention.")
        
    except Exception as e:
        logger.error(f"Single target analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
