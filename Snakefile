import os 
import sys
import pysam 
SMKDIR = os.path.dirname(workflow.snakefile) 
shell.prefix(f"source {SMKDIR}/env/env.cfg; set -eo pipefail; ")

#subread = "data/subreads/CD4_stim_DS75687.subreads.bam"
subread = "data/test_data/subreads.bam"
#ccs = "data/ccs/CD4_stim_DS75687.ccs.bam"
ccs = "data/test_data/ccs.bam"
#zmws = "data/ccs/CD4_stim_DS75687.ccs.bam.zmws"
zmws = "data/test_data/test.zmws"
# final out file
out = "test.subreads_to_ccs.bam"

ZMWS = { int(x.strip()) for x in open(zmws)}
TYPES= ["ccs", "subread"]

wildcard_constraints:
	TYPE = "|".join(TYPES),
	ZMW = "\d+",


rule all:
	input:
		out

def get_bam(wc):
	if( str(wc.TYPE) == "subread"): 
		return(subread)
	return(ccs)

rule split_bam_by_zmw:
	input:
		bam = get_bam,
	output:
		bams = temp(expand("temp/{{TYPE}}/{ZMW}.bam", ZMW=ZMWS)),
	run:
		print(input["bam"])
		bam = pysam.AlignmentFile(str(input["bam"]), check_sq=False)
		obam = None	
		preZMW = None
		counter = 0
		for rec in bam.fetch(until_eof=True):
			# check for zmw	
			if(not rec.has_tag("zm")):
				continue
			# update out file in nessisary
			curZMW = rec.get_tag("zm")
			if(preZMW != curZMW):
				if(obam is not None): obam.close()
				oname = f"temp/{wildcards.TYPE}/{curZMW}.bam"; 
				obam = pysam.AlignmentFile(oname, "wb", template=bam)
			
				counter += 1
				if(counter % 1000 == 0):
					sys.stderr.write("\rSplitting {} on ZMW: {:.02%}".format(wildcards["TYPE"], counter/len(ZMWS)))
			
			# write record 
			obam.write(rec)
			preZMW=curZMW
		obam.close()

rule pbmm2:
	input:
		rbam = "temp/ccs/{ZMW}.bam",
		bam = "temp/subread/{ZMW}.bam",
	output:
		fasta = temp("temp/aln/{ZMW}.fasta"),
		bam = temp("temp/aln/{ZMW}.bam"),
	threads: 1
	shell:"""
samtools fasta {input.rbam} > {output.fasta} 
pbmm2 align --preset SUBREAD -j {threads} --sort --no-bai {output.fasta} {input.bam} {output.bam}
"""

rule merge:
	input:
		bams = expand("temp/aln/{ZMW}.bam", ZMW=ZMWS),
	output:
		fofn = temp("temp/aln/bam.fofn"),
		bam = out,
	threads: 16
	run:
		f = open(output["fofn"], "w+")
		for bam in input["bams"]: f.write(bam + "\n")
		f.close()
		shell("samtools merge -@ {threads} -b {output.fofn} {output.bam}")
