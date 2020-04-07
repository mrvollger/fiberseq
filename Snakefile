import os 
import sys
import pysam 
import tempfile 
import numpy as np

SMKDIR = os.path.dirname(workflow.snakefile) 
SMRTBIN = "/net/eichler/vol26/projects/sequencing/pacbio/smrt-link/smrtcmds/bin"
shell.prefix(f"source {SMKDIR}/env/env.cfg; set -eo pipefail; ")

#
# define a tempdir 
#
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
else:
    TMPDIR = tempfile.gettempdir()

#
# inputs
#
configfile: "call_methyl.yaml"
subreads = os.path.abspath( config["subreads"] )
ccs = os.path.abspath( config["ccs"] )
zmw = os.path.abspath( config["zmw"] )


#
# outputs
#
workdir: config["outdir"] # all results are relative to workdir

sys.stderr.write( "{}\n{}\n{}\n{}\n".format(subreads, ccs, zmw, config["outdir"]))

#
# run parameters
#
N_BATCHES = 1000 # values much over 1000 make for a slow dag
BATCHES = [ "{:06}".format(x) for x in range(N_BATCHES) ]
THREADS = 4 
DEBUG=False

wildcard_constraints:
	B = "\d+",

def tempd(f):
	if(DEBUG): return(f)
	return(temp(f))



#
# final outputs
#
rule all:
	input:
		gff = expand("results/gff/{B}.gff.gz", B=BATCHES),
		pkl = "results/calls.csv.pkl",
		bam = "results/subreads_to_ccs.bam",

#
# create ~equal size batches of subreads
#
zmw_batch_fmt = "temp/zmw_batch/{B}.txt"
rule make_batches:
	input:
		zmw = zmw,
	output:
		zmws = temp( expand(zmw_batch_fmt, B=BATCHES) ),
	resources:
		mem = 8, 
	threads: 1 
	run:
		ZMWS = { int(x.strip()) for x in open(zmws)}
		LEN = len(ZMWS)
		counter = 0
		opos = 0
		out = open(output["zmws"][opos], "w+")
		for idx, line in enumerate(open(input["zmw"])):
			out.write(line)
			counter += 1
			if(counter >= LEN/N_BATCHES and opos < N_BATCHES - 1 ): 
				counter = 0; opos += 1
				out.close(); out = open(output["zmws"][opos], "w+")
		out.close()

#
# split ccs and subreads by zmw batches
#
def get_zmw(wc):
	return(zmw_batch_fmt.format(B=str(wc.B)))
		
rule ref_ccs:
	input:
		zmw = get_zmw,
		bam = ccs,
		pbi = ccs + ".pbi",
	output:
		ref = temp("temp/refs/{B}.fasta"),
		fai = temp("temp/refs/{B}.fasta.fai"),
		bam = temp("temp/refs/{B}.bam"),
		pbi = temp("temp/refs/{B}.bam.pbi"),
	resources:
		mem = 4, 
	threads: 1 
	shell:"""
{SMRTBIN}/bamsieve --whitelist {input.zmw} {input.bam} {output.bam}
samtools fasta {output.bam} > {output.ref}
samtools faidx {output.ref}
"""

rule subreads:
	input:
		zmw = get_zmw,
		bam = subreads,
		pbi = subreads + ".pbi",
	output:
		bam = temp("temp/subreads/{B}.bam"),
		pbi = temp("temp/subreads/{B}.bam.pbi"),
		fastq = temp("temp/subreads/{B}.fastq"),
	resources:
		mem = 4, 
	threads: 1 
	shell:"""
{SMRTBIN}/bamsieve --whitelist {input.zmw} {input.bam} {output.bam}
samtools fastq {output.bam} > {output.fastq}
"""

#
# align subreads from a zmw to the ccs of that zmw in batches
#
rule align:
	input:
		ref = rules.ref_ccs.output.ref,
		fai = rules.ref_ccs.output.fai,
		bam = rules.subreads.output.bam,
		pbi = rules.subreads.output.pbi,
		fastq = rules.subreads.output.fastq,
	output:
		bam = temp("temp/align/subreads_to_ccs.{B}.bam"),
	resources:
		mem = 8,
	threads: 1
	shell:"""
{SMKDIR}/software/minimap2/minimap2 \
	--zmw-hit-only --eqx -Y -ax map-pb -t {threads} \
	{input.ref} {input.fastq} | \
	samtools view -b -F 4 | \
	samtools sort -m {resources.mem}G -@ {threads} - | \
	pbbamify {input.ref} {input.bam}  >  {output.bam}
"""

rule index_align:
	input:
		bam = rules.align.output.bam,
	output:
		bai = temp(rules.align.output.bam + ".bai"),
		pbi = temp(rules.align.output.bam + ".pbi"),
	threads: 1
	resources:
		mem = 2,
	shell:"""
samtools index {input.bam}
pbindex {input.bam}
"""

#
# call m6A modifications per ZMW batch 
#
rule call_m6A:
	input:
		ref = rules.ref_ccs.output.ref,
		fai = rules.ref_ccs.output.fai,
		bam = rules.align.output.bam,
		bai = rules.index_align.output.bai,
		pbi = rules.index_align.output.pbi,
	output:
		csv = temp("temp/mods/{B}.csv"),
		gff = temp("temp/mods/{B}.gff"),
	resources:
		mem = 4,
	threads: THREADS
	shell:"""
{SMRTBIN}/ipdSummary {input.bam} \
	--reference {input.ref} \
	--identify m6A --methylFraction -j {threads} \
	--methylMinCov 5 \
	--gff {output.gff} --csv {output.csv}
"""			





CSVTYPES = {
	"refName" : "category",
	"tpl" : np.uint32,
	"strand" : np.uint8,
	"base" : "category",
	"score" : np.uint16,
	"tMean" : np.float32,
	"tErr" : np.float32,
	"modelPrediction": np.float32,
	"ipdRatio": np.float32,
	"coverage": np.uint16, 
	"frac": np.float32,
	"fracLow": np.float32,
	"fracUp": np.float32
}
			

rule csv_pkl:
	input:
		csv = rules.call_m6A.output.csv, 
	output:
		pkl = temp("temp/mods/{B}.csv.pkl"),
	resources:
		mem = 4, 
	threads: 1
	run:
		import pandas as pd
		tdf = pd.read_csv(input["csv"], engine="c", dtype=CSVTYPES)  
		tdf = tdf.loc[ ~ tdf.frac.isna()]
		tdf.to_pickle(output["pkl"])
	



#
# merge results
#
rule pkl_merge:
	input:
		pkls = expand(rules.csv_pkl.output.pkl, B=BATCHES),
	output:
		pkl = "results/calls.csv.pkl",
	resources:
		mem = 32,
	threads: 1
	run:
		import pandas as pd
		dfs = []
		for f in input["pkls"]:
			dfs.append(  pd.read_pickle(f)  )
			sys.stderr.write(f"\r{f}")
		sys.stderr.write("\n")
		df = pd.concat(dfs, ignore_index=False)
		df.refName = df.refName.astype("category")
		df.astype(CSVTYPES, copy=False)
		print(df.dtypes)
		df.to_pickle(output["pkl"])


rule compress_gff:
	input:
		gff = rules.call_m6A.output.gff,
	output:
		gz = "results/gff/{B}.gff.gz",
	resources:
		mem = 4,
	threads: 1
	shell: """
gzip -c {input.gff} > {output.gz}
"""
		

rule bam_merge:
	input:
		bams = expand(rules.align.output.bam, B=BATCHES)
	output:
		bam = "results/subreads_to_ccs.bam",
		fofn = temp("temp/subreads_to_ccs.fofn"),
	resources:
		mem = 1,
	threads: 4
	run:
		fofn = open(output["fofn"], "w+")
		for bam in input["bams"]: fofn.write( os.path.abspath(str(bam)) + "\n" )
		fofn.close()
		shell("head {output.fofn}")
		shell("samtools merge -@ {threads} -b {output.fofn} {output.bam}")
	



