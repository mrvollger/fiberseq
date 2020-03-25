#!/usr/bin/env python
import argparse
import os
import sys
import pysam

# global var for inputs
args=None 


def split_bam_by_zmw(bamf, outfmt, whitelist = None):
		sys.stderr.write(f"Reading {bamf}\n")
		bam = pysam.AlignmentFile(bamf, check_sq=False)
		obam = None	
		preZMW = None
		counter = 0
		for rec in bam.fetch(until_eof=True):
			# check for zmw	
			if(not rec.has_tag("zm")):
				continue
			curZMW = rec.get_tag("zm")
			
			# continue if zmw not in the whitelist
			if(whitelist is not None and curZMW not in whitelist): 
				continue

			# update out file in nessisary
			if(preZMW != curZMW):
				if(obam is not None): obam.close()
				oname = outfmt.format(ZMW=curZMW)
				
				# throw an error if it already exists 
				if(not args.force and os.path.exists(oname)): 
					raise NameError(f'{oname} already exists. Use --force to overwrite.')
				
				obam = pysam.AlignmentFile(oname, "wb", template=bam)
			
				counter += 1
				if(counter % 1000 == 0):
					sys.stderr.write("\rSplitting {} on ZMW: {}, Completed: {}".format(bamf, curZMW, counter))
			
			# write record 
			obam.write(rec)
			preZMW=curZMW
		
		obam.close()
		sys.stderr.write("\n")



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("infile", help="positional bam file")
	parser.add_argument("outformat", help="format for output ZMW files, must contains {ZMW}. e.g. out/{ZMW}.bam")
	parser.add_argument("--whitelist", help="list of ZMWs to be split", default = None)
	parser.add_argument("-n", "--number", help="numeric option", type=int, default=5)
	parser.add_argument('-f', "--force", help="overwrite existing bam files",  action="store_true", default=False)
	args = parser.parse_args()

	sys.stderr.write(f"Splitting {args.infile} into {args.outformat} by ZMW.\n")
	assert "{ZMW}" in args.outformat, "Missing {ZMW} in the output format."
	assert args.outformat[-4:] == ".bam", "Missing .bam ext in the output format."
	os.makedirs( os.path.dirname(args.outformat), exist_ok=True)		
	
	whitelist = None
	if(args.whitelist is not None):	
		whitelist = { int(x.strip()) for x in open(args.whitelist) }
		sys.stderr.write("Making {} ZMW bams from whitelist\n".format(len(whitelist)))

	split_bam_by_zmw(args.infile, args.outformat, whitelist=whitelist)
