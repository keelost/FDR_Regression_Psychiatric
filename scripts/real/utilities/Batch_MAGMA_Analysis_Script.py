#!/usr/bin/env python3
"""
Batch MAGMA Analysis Script - Final Version
Purpose: Automatically run MAGMA analysis for multiple traits based on input file names
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
        logging.FileHandler('magma_batch_analysis.log'),
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
        # Try with utf-8-sig to handle BOM in Excel-generated CSVs
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
    For example: 'st.overlap.4magma.txt' -> 'st'
                 'ever.st.overlap.4magma.txt' -> 'ever.st'
                 'sleep.dura.overlap.4magma.txt' -> 'sleep.dura'
    
    Args:
        filename: Input filename
        
    Returns:
        Phenotype name extracted from filename
    """
    # Remove the suffix '.overlap.4magma.txt'
    phenotype = filename.replace('.overlap.4magma.txt', '')
    return phenotype

def find_magma_input_files(base_dir: str, pattern: str = "*.overlap.4magma.txt") -> list:
    """
    Find all MAGMA input files in the specified directory structure.
    
    Args:
        base_dir: Base directory to search in
        pattern: File pattern to match
        
    Returns:
        List of tuples: (trait_dir_name, phenotype_name, full_file_path)
    """
    magma_files = []
    base_path = Path(base_dir)
    
    # Look for directories matching <base>/*/01.magma_input/.
    for trait_dir in base_path.glob("*/01.magma_input/"):
        if trait_dir.is_dir():
            # Extract trait name from parent directory
            trait_name = trait_dir.parent.name
            
            # Find all matching files in this directory
            for magma_file in trait_dir.glob(pattern):
                if magma_file.is_file():
                    # Extract phenotype name from filename
                    phenotype_name = extract_phenotype_from_filename(magma_file.name)
                    
                    # Create a combined trait name: directory_traitname
                    # e.g., adhd2017 + st = adhd2017_st
                    combined_trait = f"{trait_name}_{phenotype_name}"
                    
                    magma_files.append((trait_name, phenotype_name, combined_trait, str(magma_file)))
    
    logger.info(f"Found {len(magma_files)} MAGMA input files")
    return magma_files

def run_magma_analysis(trait_name: str, 
                      phenotype_name: str,
                      combined_trait: str,
                      magma_file_path: str, 
                      n_value: int, 
                      output_base_dir: str) -> bool:
    """
    Run MAGMA analysis for a single trait.
    
    Args:
        trait_name: Directory name (e.g., adhd2017)
        phenotype_name: Phenotype extracted from filename (e.g., st)
        combined_trait: Combined trait name for output (e.g., adhd2017_st)
        magma_file_path: Path to the MAGMA input file
        n_value: Sample size for the trait
        output_base_dir: Base directory for output
        
    Returns:
        True if successful, False otherwise
    """
    # Construct paths
    magma_executable = os.environ.get("FDRREG_MAGMA", "magma")
    bfile_path = os.environ.get("FDRREG_MAGMA_BFILE", "")
    gene_annot_path = os.environ.get("FDRREG_MAGMA_ANNOTATION", "")
    
    # Create output directory path: .../trait_name/04.magma_output/
    output_dir = Path(output_base_dir) / trait_name / "04.magma_output"
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
    logger.info(f"Running MAGMA for {trait_name}/{phenotype_name}: {cmd_str}")
    
    try:
        # Run the command
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Log successful completion
        logger.info(f"Completed MAGMA analysis for {trait_name}/{phenotype_name}")
        if result.stdout:
            logger.debug(f"STDOUT: {result.stdout[:500]}...")
            
        return True
        
    except subprocess.CalledProcessError as e:
        logger.error(f"MAGMA failed for {trait_name}/{phenotype_name}: {e}")
        logger.error(f"STDERR: {e.stderr[:500] if e.stderr else 'No stderr'}")
        return False
        
    except Exception as e:
        logger.error(f"Unexpected error for {trait_name}/{phenotype_name}: {e}")
        return False

def main():
    """
    Main function to orchestrate batch MAGMA analysis.
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Batch MAGMA Analysis')
    parser.add_argument('--base-dir', 
                       default=os.environ.get("FDRREG_RESULTS_DIR", "data/pipeline"),
                       help='Base directory containing trait subdirectories')
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
    
    args = parser.parse_args()
    
    # Update logging level
    logger.setLevel(getattr(logging, args.log_level))
    
    # Start time
    start_time = datetime.now()
    logger.info(f"Starting batch MAGMA analysis at {start_time}")
    
    try:
        # Step 1: Load trait N values
        trait_n_mapping = load_trait_n_values(args.n_file)
        
        # Step 2: Find all MAGMA input files
        magma_files = find_magma_input_files(args.base_dir)
        
        if not magma_files:
            logger.warning("No MAGMA input files found. Check directory structure.")
            return
        
        # Step 3: Process each trait
        successful = 0
        failed = 0
        skipped = 0
        
        # Prepare tasks for parallel execution
        tasks = []
        for trait_name, phenotype_name, combined_trait, magma_file_path in magma_files:
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
                output_dir = Path(args.base_dir) / trait_name / "04.magma_output"
                output_file = output_dir / f"{phenotype_name}"
                cmd_str = f"{os.environ.get('FDRREG_MAGMA', 'magma')} --bfile {os.environ.get('FDRREG_MAGMA_BFILE', '')} --pval {magma_file_path} use=snpid,p.decor N={n_value} --gene-annot {os.environ.get('FDRREG_MAGMA_ANNOTATION', '')} --out {output_file}"
                print(f"DRY RUN: {cmd_str}")
                continue
                
            tasks.append((trait_name, phenotype_name, combined_trait, magma_file_path, n_value, args.base_dir))
        
        if args.dry_run:
            logger.info(f"Dry run complete. Would process {len(tasks)} traits")
            return
        
        # Execute tasks in parallel
        if tasks:
            logger.info(f"Processing {len(tasks)} traits with {args.max_workers} workers")
            
            with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
                # Submit all tasks
                future_to_trait = {}
                for trait_name, phenotype_name, combined_trait, magma_file_path, n_value, output_base_dir in tasks:
                    future = executor.submit(
                        run_magma_analysis,
                        trait_name,
                        phenotype_name,
                        combined_trait,
                        magma_file_path,
                        n_value,
                        output_base_dir
                    )
                    future_to_trait[future] = f"{trait_name}/{phenotype_name}"
                
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
        logger.info(f"BATCH ANALYSIS COMPLETE")
        logger.info(f"{'='*60}")
        logger.info(f"Total MAGMA input files found: {len(magma_files)}")
        logger.info(f"Successfully processed: {successful}")
        logger.info(f"Failed: {failed}")
        logger.info(f"Skipped (missing N values): {skipped}")
        logger.info(f"Total duration: {duration}")
        logger.info(f"Log file: magma_batch_analysis.log")
        logger.info(f"{'='*60}\n")
        
        # Show failed tasks summary if any
        if failed > 0:
            logger.info("Some tasks failed. Check the log file for details.")
            logger.info("Failed tasks may need manual attention.")
        
    except Exception as e:
        logger.error(f"Batch analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
