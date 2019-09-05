#!/usr/bin/env python 
# ./script Archaea 1
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


bins_ct = 100




def overlap(srna_st, srna_fi, gene_info, array):
	st = gene_info[0]
	fi = gene_info[1]
	strand = gene_info[2]
	gene_len = fi-st

	for i in range(bins_ct*3):
		n = i*1.0/(bins_ct*3) * gene_len*3 + st - gene_len
		if n>=srna_st and n<=srna_fi:
			if strand=="+":
				array[i]+=1
			elif strand=="-":
				array[bins_ct*3-i-1]+=1
			else:
				print "impossible"
				quit()
	return array


def start_position(srna_st, srna_fi, gene_info, array):
	st = gene_info[0]
	fi = gene_info[1]
	gene_len = fi-st

	strand = gene_info[2]
	if strand=="+":
		fiveprime = srna_fi
	else:
		fiveprime = srna_st



	relative_st = (fiveprime-st)*bins_ct/gene_len + bins_ct
	if relative_st>=3*bins_ct:
		relative_st=3*bins_ct-1
	elif relative_st<=0:
		relative_st=0

	if strand=="+":
		array[relative_st]+=1
	elif strand=="-":
		array[3*bins_ct-relative_st-1]+=1
	else:
		print "impossible"
		quit()
	return array
	
	
###################### MAIN ########################
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)
fig, axs = plt.subplots(2, 3, figsize=(12,8))


taxonomy={}
for line in open("contig_taxonomy.tab"):
        cut=line.strip().split("\t")
	if len(cut)<2: continue
        taxonomy[cut[0]]=cut[1]


labels=[['A','B','C'],['D','E','F']]

for y,taxon in enumerate(["Bacteria", "Archaea"]):
	for x,plot in enumerate([1,2,3]):
		label = labels[y][x]
		print label
		genes={}
		for line in open("img_annotation.gff"):
			cut=line.strip().split("\t")
			contig = cut[0]
			if contig not in taxonomy:
				continue
			if taxon not in taxonomy[contig]:
				continue
	
			strand = cut[6]
			st = int(cut[3])
			fi = int(cut[4])
			ID = cut[8].split(";")[0]
			genes[ID] = (st, fi, strand)


		overlaps = listofzeros = [0]*bins_ct*3
		starts = listofzeros = [0]*bins_ct*3
		for line in open("small_antisense_ncRNAs.gff"):
			if "antisense_to_gene" not in line:
				continue
			cut=line.strip().split("\t")
			st = int(cut[3])
			fi = int(cut[4])
			strand = cut[6]
			ID = cut[8].split(";")[-2].split('"')[1]
			if ID not in genes:
				continue
			gene_info = genes[ID]

			if plot==1:
				overlaps = overlap(st, fi, gene_info, overlaps)
			elif plot==2:
				starts = start_position(st, fi, gene_info, starts)
			elif plot==3:
				overlaps = overlap(gene_info[0], gene_info[1], (st,fi,strand), overlaps)
	

		if plot==1 or plot==3:
			data=overlaps
		if plot==2:
			data=starts
		
		
		axs[y,x].scatter(range(-bins_ct,2*bins_ct), data, c='b', alpha=0.4, s=8)
		axs[y,x].annotate(label, xy=(-0.12, 1.0), xycoords="axes fraction", fontsize=20)
		axs[y,x].set_xticks([0, bins_ct])
		axs[y,x].set_xticklabels(["Start (5')", "End (3')"], rotation='horizontal', fontsize=14)
		axs[y,x].grid(linestyle="--", alpha=0.8)
		axs[y,x].spines['right'].set_visible(False)
		axs[y,x].spines['top'].set_visible(False)
		
		if label=="A":
			axs[y,x].set_title("asRNA overlaps on genes")
		elif label=="B":
			axs[y,x].set_title("asRNA start sites on genes")
		elif label=="C":
			axs[y,x].set_title("Gene overlaps on asRNAs")
		
		if label=="A":
			axs[y,x].set_ylabel("Sequence enrichment in Bacteria")
		elif label=="D":
			axs[y,x].set_ylabel("Sequence enrichment in Archaea")

		


		
plt.tight_layout(w_pad=0.5, h_pad=1)
plt.savefig("figure.png", dpi=600)

