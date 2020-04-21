#!/bin/bash
set -eo pipefail 

snakemake --restart-times 1 -w 60 -j 200 --drmaa " -l mfree={resources.mem}G -pe serial {threads} -l h_rt=128:00:00 -V -cwd -S /bin/bash " --drmaa-log-dir logs -k $@ 

#snakemake --report results/results.html


