#!/usr/bin/env python3

import os
import sys

def main():
    
    index_path = sys.argv[1]
    
    # Check if directory exists
    if not os.path.exists(index_path):
        raise ValueError(f'Error: STAR index directory not found: {index_path}')

    required_files = [
        "chrLength.txt",
        "chrNameLength.txt", 
        "chrName.txt",
        "chrStart.txt",
        "exonGeTrInfo.tab",
        "exonInfo.tab",
        "geneInfo.tab",
        "Genome",
        "genomeParameters.txt",
        "Log.out",
        "SA",
        "SAindex",
        "sjdbInfo.txt",
        "sjdbList.fromGTF.out.tab",
        "sjdbList.out.tab",
        "transcriptInfo.tab"
    ]
    
    # Get all files in the index directory
    found_files = os.listdir(index_path)
    missing_files = [f for f in required_files if f not in found_files]
    
    # Exit if missing files present
    if len(missing_files) > 0:
        raise ValueError(f'Missing files in STAR index: {missing_files}')

if __name__ == "__main__":
    main()
    
