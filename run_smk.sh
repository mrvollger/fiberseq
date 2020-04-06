#!/bin/bash

snakemake --restart-times 1 -w 60 -j 500 --drmaa " -l mfree={resources.mem}G -pe serial {threads} -l h_rt=128:00:00 -V -cwd -S /bin/bash " --drmaa-log-dir logs $@ 

