import os 
import sys
import pysam 
import tempfile 

SMKDIR = os.path.dirname(workflow.snakefile) 
SMRTBIN = "/net/eichler/vol26/projects/sequencing/pacbio/smrt-link/smrtcmds/bin"
shell.prefix(f"source {SMKDIR}/env/env.cfg; set -eo pipefail; ")

# define a tempdir 
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
else:
    TMPDIR = tempfile.gettempdir()



#subread = "data/subreads/CD4_stim_DS75687.subreads.bam"
subreads = "data/test_data/subreads.bam"
#ccs = "data/ccs/CD4_stim_DS75687.ccs.bam"
ccs = "data/test_data/ccs.bam"
#zmws = "data/ccs/CD4_stim_DS75687.ccs.bam.zmws"
zmws = "data/test_data/test.zmws"
# final out file
out = "temp/test.subreads_to_ccs.bam"
# number of batches to split the process into
N_BATCHES = 200
BATCHES = list(range(N_BATCHES))
THREADS=4 # threads per batch per job

wildcard_constraints:
	B = "\d+",

DEBUG=False
def tempd(f):
	if(DEBUG): return(f)
	return(temp(f))


rule all:
	input:
		out = expand("temp/mods/{B}.gff", B=BATCHES),

#
# create ~equal size batches of subreads
#
zmw_batch_fmt = "temp/zmw_batch/{B}.txt"
rule make_batches:
	input:
		zmw = zmws,
	output:
		zmws = tempd( expand(zmw_batch_fmt, B=BATCHES) ),
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
		ref = tempd("temp/refs/{B}.fasta"),
		fai = tempd("temp/refs/{B}.fasta.fai"),
		bam = temp("temp/refs/{B}.bam"),
		pbi = temp("temp/refs/{B}.bam.pbi"),
	resources:
		mem = 8, 
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
		mem = 8, 
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
		bam = "temp/align/subreads_to_ccs.{B}.bam"
	resources:
		mem=8
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
		bai = rules.align.output.bam + ".bai",
		pbi = rules.align.output.bam + ".pbi",
	threads: 1
	resources:
		mem=8
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
		csv = "temp/mods/{B}.cvs",
		gff = "temp/mods/{B}.gff",
	resources:
		mem=8
	threads: THREADS
	shell:"""
{SMRTBIN}/ipdSummary {input.bam} \
	--reference {input.ref} \
	--identify m6A --methylFraction -j {threads} \
	--methylMinCov 5 \
	--gff {output.gff} --csv {output.csv}
"""			


