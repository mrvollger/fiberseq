#!/usr/bin/env python
import argparse
import os 
import sys
import pysam


# global var for inputs
args=None 

def get_zmw(s):
	mv, zmw = s.split("/")[0:2]
	return(mv, zmw)

def validate_bam(f):
	bam = pysam.AlignmentFile(f)
	count = 0
	for rec in bam.fetch(until_eof=True):
		mv, zmw = get_zmw(rec.query_name)
		rmv, rzmw = get_zmw(rec.reference_name)
		if(mv == rmv or zmw == rzmw ):
			count += 1
			sys.stderr.write(f"\rValid alignments: {count}")
		else:
			sys.stdout.write(f"Invalid zmw alignment: {ref.reference_name} {ref.query_name}\n")

	sys.stderr.write(f"\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("infile", help="positional input")
	parser.add_argument("-s", "--string", help="string option")
	parser.add_argument("-n", "--number", help="numeric option", type=int, default=5)
	parser.add_argument("-l", "--list", nargs="*", help="list with zero or more entries")
	parser.add_argument("-l2", "--list2", nargs="+", help="list one or more entries")
	parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
	args = parser.parse_args()
	validate_bam(args.infile)

