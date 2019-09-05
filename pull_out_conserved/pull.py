#!/usr/bin/env python
import sys
soi = {}
for line in open(sys.argv[1]):
	soi[line.strip()]=None


seq={}
for line in open("srnas.fa"):
	if line[0]==">":
		name=line[1:-1]
	else:
		seq[name]=line.strip()


for line in open("srnas.gff"):
	cut=line.strip().split("\t")
	info=cut[8].split(";")
	for f in info:
		if "transcript_id" in f:
			ID = f.split('"')[1]
	if ID in soi:
		name=cut[0]+":"+str(int(cut[3])-1)+"-"+cut[4]+"("+cut[6]+")"
		print ID+"\t"+seq[name]



