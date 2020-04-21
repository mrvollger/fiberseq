#!/usr/bin/env python
import argparse
import os 
import sys
import pysam
import re
import pandas as pd

# global var for inputs
args=None 

# MIN_QUAL_TBL
QUAL_TBL = {
	12:13,
	13:14,
	14:14,
	15:15,
	16:16,
	17:16,
	18:17,
	19:17,
	20:18,
	21:18,
	22:19,
	23:19,
	24:20,
	25:21,
	26:21,
	27:22,
	28:22,
	29:23,
	30:23,
	31:23,
	32:23,
	33:23,
	34:23,
	35:23,
	36:23,
	37:23
}

def min_qual(np):
	if( np < 12):
		return(13)
	if(np > 27):
		return(23)
	return( QUAL_TBL[np] )

def qual_to_char(num):
    return( chr( min(int(num) + 33, 126)) )

def make_diff_str(l):
    pre = 0; rtn = ""
    total = 0
    for x in l:
        inc = x-pre
        rtn += ",{}".format(inc); pre = x+1
        total += inc+1
    assert total == l[-1] +1 
    return(rtn)

# NEED to have someone validate this logic based on :
# https://404-3666509-gh.circle-artifacts.com/0/root/project/pdfs/SAMtags.pdf
def get_mod_tags(gff, read_name, aln, mod_types=["m6A"]):
	pat = re.compile("([^=;]+)=(\d+\.*\d*);*")    
	qualf = "" # phred values for forward strand modification calls
	qualr = "" # phred values for reverse strand modification calls
	posf = []
	posr = []

	positions = []
	qualities = []

	for rec in gff.fetch(reference=read_name):
		#if(rec.feature in mod_types):
			#ms = re.findall(pat, rec.attributes); atts = { atr: float(val)  for atr, val in ms }
			#if("identificationQv" in atts):
				#qual_char = qual_to_char(atts["identificationQv"]) 
		start = rec.start
		strand = rec.strand
		base = aln.query_sequence[rec.start]

		if( (rec.score >= min_qual(aln.get_tag("np")))
				and (strand == "+" and base == "A") or (strand == "-" and base == "T") ):
			qual_char = qual_to_char(rec.score) 
			
			# add to total list
			positions.append(start)
			qualities.append(rec.score)
			
			# add to strand specific lists 
			if(strand == "+"): 
				qualf += qual_char
				posf.append(start)
			elif(strand == "-"):
				qualr += qual_char
				posr.append(start)
			else:
				raise("The strand must be indicated.")
		
	#
	# Modify sequence with "W"
	#
	vals = set(posf + posr)
	qual = aln.query_qualities
	seq = ""
	for idx, bp in enumerate(aln.query_sequence):
		if(idx in vals): 
			seq += "W"
		else:
			seq += bp
	aln.query_sequence = seq
	aln.query_qualities = qual


	#
	# make tags 
	#
	MM = ""; MP = ""
	if(len(posf) > 0):
		MM += "A+a{};".format(make_diff_str(posf))
		MP += "{}".format(qualf)
	if(len(posr) > 0):
		MM += "T-a{};".format(make_diff_str(posr))
		MP += " {}".format(qualr)

	if(len(MM) == 0):
		return(None)
	else:
		return( (MM, MP, "153,255,204", positions, qualities) )


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("ccs", help="positional input, must be ccs bam")
	parser.add_argument("gff", help="positional input, must be a tabix gff file with tbi")
	parser.add_argument("outbam", help="positional output, must be bam")
	parser.add_argument("pkl", help="output pkl file")
	parser.add_argument("-n", "--np", help="minimum number of passes to mark", type=int, default=5)
	parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
	args = parser.parse_args()

	bam = pysam.AlignmentFile(args.ccs, check_sq=False)
	gff =  pysam.TabixFile(args.gff, parser=pysam.asGFF3())
	obam = pysam.AlignmentFile(args.outbam, "wb", template=bam)

	contigs = set(gff.contigs)
	tmp_d = { "name":[], "np":[], "length":[], "positions":[], "qualities":[] }	
	for idx, aln in enumerate(bam.fetch(until_eof=True)):
		if( (aln.query_name in contigs) and (aln.get_tag("np") >= args.np) ):
			new_tags = get_mod_tags(gff, aln.query_name, aln)
			if(new_tags is not None):
				aln.set_tag("MM", new_tags[0])
				aln.set_tag("MP", new_tags[1])
				aln.set_tag("yc", new_tags[2])
				
				tmp_d["name"].append(aln.query_name) 
				tmp_d["np"].append(aln.get_tag("np")) 
				tmp_d["length"].append(aln.query_length) 
				tmp_d["positions"].append(new_tags[3]) 
				tmp_d["qualities"].append(new_tags[4]) 

		# wirte (modified) aln to out
		obam.write(aln)	
		sys.stderr.write(f"\rReads done: {idx+1}")
	
	df = pd.DataFrame(tmp_d)
	df.to_pickle(args.pkl)
	
	sys.stderr.write(f"\rDone\n")
	obam.close()
	bam.close()
	gff.close()


