#!/usr/bin/env python
import sys


brite2func={}
for line in open("brite2function.tab"):
        cut=line.strip().split("\t")
        if "Drug Dev" in cut[1] or "Human Dis" in cut[1] or "Organismal Sys" in cut[1]: continue
        brite2func[cut[0]]=cut[1]

ko2brite={}
for line in open("ko2brite.tab"):
	cut=line.strip().split("\t")
	brites=cut[1].split(";")
	good_brites=[]
	for i in range(len(brites)):
		brite=brites[i]
		if brite[2:] in brite2func:
			good_brites.append(brite)
	cut[1]=";".join(good_brites)
	if cut[1]!="":
		ko2brite[cut[0]]=cut[1]


GaIDs={}
for line in open("img_annotation.gff"):
	cut=line.strip().split("\t")
	fields = cut[8].split(";")
	for f in fields:
		if f.startswith("ID"):
			ID = f.split("=")[1]
		if f.startswith("locus_tag"):
			gene = f.split("=")[1]
	GaIDs[gene]=ID



print "\t".join(["#Locus", "ID", "GaID", "Gene Product", "BRITE Pathways", "Functions"])
for line in open("img_annotation.products"):
	cut=line.strip().split("\t")
	cut.append(GaIDs[cut[0]])
	if "KO" in cut[2]:
		if cut[2].split(":")[1] in ko2brite:
			cut.append(ko2brite[cut[2].split(":")[1]])
			brites = ko2brite[cut[2].split(":")[1]].split(";")
			functions=[]
			for brite in brites:
				functions.append(brite2func[brite[2:]])
			cut.append("|".join(functions))
	else:
		cut.append("NA")
		cut.append("NA")
	ID=cut[2]; name=cut[1]
	cut[1]=ID; cut[2]=name
	print "\t".join(cut)



