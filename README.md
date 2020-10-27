# reference free fiber-seq calls

This has a pipline for calling accesibility per ZMW. 

This pipeline uses snakemake and must have a config file (`call_methyl.yaml`) with these options speficied:
```
outdir: "methyl_calls_CD4_stim_DS76351/"
ccs: "data/ccs/CD4_stim_DS76351.ccs.bam"
zmw: "data/ccs/CD4_stim_DS76351.ccs.bam.zmws"
subreads: "data/subreads/CD4_stim_DS76351.subreads.bam"
```

It also uses a special fork of minimap2 for a `--zmw-hit-only` option. To install that cd into `software` and follow the commands in `software/commands.sh`.

### This pipeline does have references to software specifc to the Eichler cluster, so it will not work elsewhere without modification. 




