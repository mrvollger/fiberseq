#!/bin/bash

snakemake -w 60 -j 500 -p --drmaa " -l mfree={resources.mem}G -pe serial {threads} -l h_rt=128:00:00 -V -cwd -S /bin/bash " --drmaa-log-dir temp/logs $@ 

