#!/usr/bin/python

import os
from itertools import product

path_wd = '/data/cossa2_dare/MI_TO_bench'
path_pipeline = '/data/cossa2/phylo_inference'
path_output = os.path.join(path_wd, 'results/LAST_ALL_TOGETHER/tune_mito_maegatk')

cmd = f'nextflow run {path_pipeline}/main.nf -c {path_wd}/nextflow.config -params-file {path_wd}/params/LAST_ALL_TOGETHER/params_tune_mito_maegatk.json -profile conda_garr -entry tune'

MIN_CELL_NUMBER = [1,10]
DISTANCES = ['weighted_jaccard', 'jaccard']
combos = list(product(MIN_CELL_NUMBER, DISTANCES))

for i, (min_cell_number, distance) in enumerate(combos):
    
    print(f'{i}/{len(combos)}')
    print(i, min_cell_number, distance)

    outdir = os.path.join(path_output, f'{min_cell_number}_{distance}')
    new_cmd = cmd + f' --outdir {outdir} --min_cell_number {min_cell_number} --distance_metric {distance}'
    
    os.system(new_cmd)
