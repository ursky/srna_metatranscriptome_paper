#!/usr/bin/env python
import sys
genes={}
for line in open("img_annotation.master"):
        cut=line.strip().split("\t")
        genes[cut[3]]=line.strip()

for line in open(sys.argv[1]):
	if line.strip() in genes:
		print genes[line.strip()]
	else:
		print "NA"
