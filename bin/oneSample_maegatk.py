#!/usr/bin/env python3

"""
maegatk for one cell (MAESTER): read filtering, fgbio UMI consensus, re-mapping, allele counts.
Usage: oneSample_maegatk.py <cell.bam> <cell> <mt.fa> [options]
Writes <cell>.{A,C,G,T,coverage,depth}.txt in the working directory.
The steps are those of https://github.com/caleblareau/maegatk and the MAESTER paper (Miller et al., 2022).
"""

import os
import argparse
import pysam


##


parser = argparse.ArgumentParser(description='maegatk for one cell')
parser.add_argument('bam')
parser.add_argument('sample', help='cell name, used as prefix of the outputs')
parser.add_argument('fasta', help='MT FASTA (a single contig, bwa-indexed)')
parser.add_argument('--umi_tag', default='UB')
parser.add_argument('--NHmax', type=int, default=2)
parser.add_argument('--NMmax', type=int, default=15)
parser.add_argument('--min_reads', type=int, default=3, help='reads per UMI family for the fgbio consensus')
parser.add_argument('--base_qual', type=int, default=0)
parser.add_argument('--alignment_quality', type=int, default=0)
parser.add_argument('--proper_paired', action='store_true')
args = parser.parse_args()

inputbam = args.bam
sample = args.sample
fasta_file = args.fasta

# Contig name and length from the FASTA
with open(fasta_file) as fh:
    _lines = [x.strip() for x in fh if x.strip()]
_names = [x[1:].split()[0] for x in _lines if x.startswith('>')]
if len(_names) != 1:
    raise ValueError(f'{fasta_file}: expected a single contig, found {len(_names)}')
mito_genome = _names[0]
mito_length = str(sum(len(x) for x in _lines if not x.startswith('>')))

umi_barcode = args.umi_tag
base_qual = str(args.base_qual)
alignment_quality = str(args.alignment_quality)
NHmax = str(args.NHmax)
NMmax = str(args.NMmax)
min_reads = str(args.min_reads)
proper_paired = str(bool(args.proper_paired))


# Software paths
python = "python"
fgbio = 'fgbio'

# Script locations: oneSample_maegatk.py and other script are placed in the same folder
script_dir = os.path.dirname(os.path.abspath(__file__))
filtclip_py = os.path.join(script_dir, "filterClipBam.py")
sumstatsBP_py = os.path.join(script_dir, "sumstatsBP.py")

# Prepare intermediate filenames: skip of the original script .replace synthax needed by Snakemake 
filtlog = 'rmlog.log'
filtlog = 'filter.log'
temp_bam0 = 'temp0.bam'
temp_bam1 = 'temp1.bam' 
temp_sam15 = 'sam15.sam'
temp_bam15 = 'bam15.bam'
temp_bam2 = 'bam2.bam'
temp_fastq = 'temp0.fastq'
outputbam = 'consensus.bam' # The last one, actually used for pileup by "sumstatsBP.py" 
prefixSM = sample
outputdepth = sample + ".depth.txt"


##


# From here on, the code is exactly as it is in https://github.com/caleblareau/maegatk

# 1) Filter bam files
pycall = " ".join([python, filtclip_py, inputbam, filtlog, mito_genome, proper_paired, NHmax, NMmax]) + " > " + temp_bam0
os.system(pycall)

# 2) Sort the filtered bam file
fgcallone =  fgbio + " GroupReadsByUmi -s Identity -e 0 -i " + temp_bam0 + " -o " + temp_bam1 + " -t " + umi_barcode
os.system('echo "'+fgcallone+'"')
os.system(fgcallone)

# 2.5) Modify the UB tag
samtoolscall = 'samtools view -H ' + temp_bam1 + '> ' + temp_sam15 + '; samtools view ' + temp_bam1 + '| awk \'OFS="\t" {$13=$13""$4; print $0}\' >> ' + temp_sam15 + '; samtools view -b ' + temp_sam15 + '> ' + temp_bam15
os.system('echo "'+samtoolscall+'"')
os.system(samtoolscall)

# 3) Call consensus reads
fgcalltwo = fgbio + " CallMolecularConsensusReads -t "+umi_barcode+" -i "+temp_bam15+" -o " + temp_bam2 +" -M " + min_reads
os.system(fgcalltwo)

# 4) Convert consensus bam to fastq
# bedtools_call = "bedtools bamtofastq -i "+ temp_bam2 +" -fq " + temp_fastq # Bedtools stopped working for some reason, replacing it with samtools fastq
samtoolscall2 = 'samtools fastq -T cM ' + temp_bam2 + " | sed 's/\tcM:i:/_/g' > " + temp_fastq
os.system(samtoolscall2)

# 5) Remap + sort bam files
bwa_call = "bwa mem " + fasta_file + " " + temp_fastq + " | samtools sort -o "+ outputbam +" -"
os.system(bwa_call)
pysam.index(outputbam)

# 6) Get allele counts per sample / base pair and per-base quality scores
alleleCountcall = " ".join([python, sumstatsBP_py, outputbam, prefixSM, mito_genome, mito_length, base_qual, sample, fasta_file, alignment_quality])
os.system(alleleCountcall)

# 7) Get depth from the coverage sparse matrix
with open(prefixSM + ".coverage.txt", 'r') as coverage:
	depth = 0
	for row in coverage:
		s = row.split(",")
		depth += int(s[2].strip())

with open(outputdepth, 'w') as d:
	d.write(sample + "\t" + str(round(float(depth)/float(mito_length),2)) + "\n")


## 