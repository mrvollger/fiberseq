#!/usr/bin/env python
import argparse
import os 
import sys
import pysam
import re

# global var for inputs
args=None 

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
def get_mod_tags(gff, read_name, mod_types=["m6A"]):
    pat = re.compile("([^=;]+)=(\d+\.*\d*);*")    
    qualf = "" # phred values for forward strand modification calls
    qualr = "" # phred values for reverse strand modification calls
    posf = []
    posr = []
    for rec in gff.fetch(reference=read_name):
        if(rec.feature in mod_types):
            ms = re.findall(pat, rec.attributes); atts = { atr: float(val)  for atr, val in ms }
            if("identificationQv" in atts):
                qual_char = qual_to_char(atts["identificationQv"]) 
                start = rec.start
                strand = rec.strand
                if(strand == "+"): 
                    qualf += qual_char
                    posf.append(start)
                elif(strand == "-"):
                    qualr += qual_char
                    posr.append(start)
                else:
                    raise("The strand must be indicated.")
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
        return({"MM":MM, "MP":MP})


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("ccs", help="positional input, must be ccs bam")
	parser.add_argument("gff", help="positional input, must be a tabix gff file with tbi")
	parser.add_argument("outbam", help="positional output, must be bam")
	parser.add_argument("-n", "--np", help="minimum number of passes to mark", type=int, default=5)
	parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
	args = parser.parse_args()

	bam = pysam.AlignmentFile(args.ccs, check_sq=False)
	gff =  pysam.TabixFile(args.gff, parser=pysam.asGFF3())
	obam = pysam.AlignmentFile(args.outbam, "wb", template=bam)

	contigs = set(gff.contigs)
	
	for aln in bam.fetch(until_eof=True):
		# skip reads that had no modifications 
		if( (aln.query_name not in contigs) or (aln.get_tag("np") < args.np) ):
			obam.write(aln)	
			continue 
		
		new_tags = get_mod_tags(gff, read_name=aln.query_name)
		if(new_tags is not None):
			for tag in new_tags: aln.set_tag(tag, new_tags[tag])
		
		obam.write(aln)	

	obam.close()
	bam.close()
	gff.close()


