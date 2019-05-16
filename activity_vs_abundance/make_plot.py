#!/usr/bin/env python
import sys
import pickle
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

def load_file(filename):
	out={}
	for line in open(filename):
		if "NODE" not in line:
			continue
		cut=line.strip().split("\t")
		contig=cut[0]
		length = int(contig.split("_")[3])
		if length>=3000:
			out[contig]=float(cut[1])
	return out


def load_contig_abundances():
	out = {}
	for timepoint in ["1","9","13","17"]:
		for day in ["a", "b"]:
			if timepoint=="9" and day=="a":
				continue
			for replicate in ["1","2","3"]:
				sample = timepoint+day+'-'+replicate
				filename = sample+".quant"
				print "loading from "+filename
				data = load_file(filename)
				out[sample]=data
	return out


def load_srna_abundances(filename):
	print "loading from "+filename
	out = {}
	for i,line in enumerate(open(filename)):
		cut = line.strip().split("\t")
		if i==0:
			head=cut
			for sample in head[1:]:
				out[sample]={}
		else:
			for j,val in enumerate(cut):
				if j==0:
					srna=".".join(val.split(".")[:2])
				else:
					if val=="NA":
						val=0
					sample=head[j]
					out[sample][srna]=float(val)
	return out


def load_gtf(filename):
	print "loading from "+filename
	out={}
	for line in open(filename):
		cut=line.strip().split("\t")
		contig=cut[0]
		info=cut[8].split(";")
		for f in info:
			if "transcript_id" in f:
				srna = ".".join(f.split('"')[1].split(".")[:-1])
				continue
		out[srna]=contig
	return out


def load_taxa(filename):
	taxa={}
	for line in open(filename):
		cut=line.strip().split("\t")
		if len(cut)!=2:
			continue
		taxa[cut[0]]=cut[1].split(";")
	return taxa


def get_averages(srnas, contig_abundances, srna_abundances, timepoint, taxa):
	data = {}
	print "calculating averages for "+timepoint
	for srna in srnas:
		data[srna]=([],[])
		contig=srnas[srna]
		for sample in contig_abundances:
			if not sample.startswith(timepoint):
				continue
			contig_abund = contig_abundances[sample][contig]
			srna_abund = srna_abundances[sample][srna]
			data[srna][0].append(contig_abund)
			data[srna][1].append(srna_abund)
	out_contigs=[]
	out_srnas=[]
	colors=[]
	labels={"Euryarchaeota":"blue", "Bacteroidetes":"red", "Cyanobacteria":"cyan", "Other":"black"}
	
	for srna in data:
		contig = srnas[srna]
		if contig in taxa:
			taxon=taxa[contig]
			if "Eukaryota" in taxon:
				continue
			if len(taxon)<2:
				c="black"
			elif taxon[1] in labels:
				c=labels[taxon[1]]
			else:
				c="black"
		else:
			c="black"

		ave_contig = np.mean(data[srna][0])
		ave_srna = np.mean(data[srna][1])
		out_contigs.append(ave_contig)
		out_srnas.append(ave_srna)
		colors.append(c)
	return out_contigs, out_srnas, colors



######## MAIN ########

if "reload" in sys.argv:
	contig_abundances = load_contig_abundances()
	print "saving contig_abundances.pkl"
	output = open('contig_abundances.pkl', 'wb')
	pickle.dump(contig_abundances, output)
	output.close()
else:
	print "loading contig_abundances.pkl"
	pkl_file = open('contig_abundances.pkl', 'rb')
	contig_abundances = pickle.load(pkl_file)
	pkl_file.close()


srna_abundances = load_srna_abundances("all_meta-sRNAs_expression.txt")
intergenic = load_gtf("small_intergenic_ncRNAs.gff")
antisense = load_gtf("small_antisense_ncRNAs.gff")
taxa = load_taxa("contig_taxonomy.tab")

##############  PLOTTING ###############
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)
fig, axs = plt.subplots(2, 4, figsize=(12,6))
times={"1":"1am","9":"9am","13":"1pm","17":"5pm"}
for y,timepoint in enumerate(["1","9","13","17"]):
	for x,srna_type in enumerate([intergenic, antisense]):
		print x, y
		ax = axs[x,y]
		ax.loglog([0.001, 10000], [0.1, 1000000], '--', alpha=0.5)
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_xlim(0.001, 10000)
		ax.set_ylim(0.1, 1000000)
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.grid(linestyle='--', alpha=0.5)

		if y==0 and x==0:
			ax.set_ylabel("sRNA expression (tpm)")
			ax.annotate("A", xy=(-0.3, 1.02), xycoords="axes fraction", fontsize=18)
		if y==0 and x==1:
			ax.set_ylabel("sRNA expression (tpm)")
			ax.annotate("B", xy=(-0.3, 1.02), xycoords="axes fraction", fontsize=18)
		if x==0:
			ax.set_title(times[timepoint])
		if x==1:
			ax.set_xlabel("Total contig activity (tpm)")

		contig_averages,srna_average,colors = get_averages(srna_type, contig_abundances, srna_abundances, timepoint, taxa)
		ax.scatter(contig_averages, srna_average, c=colors, alpha=0.1, s=5)


print "making legend..."
labels={"Euryarchaeota":"blue", "Bacteroidetes":"red", "Cyanobacteria":"cyan", "Other":"black"}
ax = fig.add_axes([0.0, 0.0, 1, 0.05])
ax.axis("off")
legend_elements = [Patch(facecolor=labels["Euryarchaeota"], edgecolor='k', label="Euryarchaeota", linewidth=1),
        Patch(facecolor=labels["Bacteroidetes"], edgecolor='k', label="Bacteroidetes", linewidth=1),
        Patch(facecolor=labels["Cyanobacteria"], edgecolor='k', label="Cyanobacteria", linewidth=1),
        Patch(facecolor=labels["Other"], edgecolor='k', label="Other", linewidth=1)]
ax.legend(handles=legend_elements, loc="lower center", framealpha=1, frameon=True, facecolor='w', ncol=4, columnspacing=1, handlelength=1, prop={'size': 12})


plt.tight_layout(w_pad=0.5, h_pad=1)
plt.subplots_adjust(bottom=0.18)
plt.savefig("figure.png", dpi=600)






