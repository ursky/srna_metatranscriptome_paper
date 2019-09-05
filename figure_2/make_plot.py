#!/usr/bin/env python
import sys
import pickle
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
from scipy import stats

def load_file(filename, contigs1, contigs2):
	print "loading data from "+filename
	out={}
	for line in open(filename):
		if "NODE" not in line:
			continue
		cut=line.strip().split("\t")
		contig=cut[0]
		length = int(contig.split("_")[3])
		if contig not in contigs1 and contig not in contigs2:
			continue
		out[contig]=float(cut[1])
	return out


def load_contig_abundances(contigs1, contigs2):
	out = {}
	for timepoint in ["01","09","13","17"]:
		for day in ["a", "b"]:
			if timepoint=="09" and day=="a":
				continue
			for replicate in ["1","2","3"]:
				sample = timepoint+day+'-'+replicate
				filename = "CONTIG_TOTAL_EXPRESSION/"+sample+".quant.counts"
				data = load_file(filename, contigs1, contigs2)
				out[sample]=data
	for timepoint in ["9am", "9pm"]:
		for replicate in ["1","2","3","4","5","6","7","8","9","10","11","12"]:
			sample = timepoint+'-'+replicate
			filename = "CONTIG_TOTAL_EXPRESSION/"+sample+".quant.counts"
			data = load_file(filename, contigs1, contigs2)
			out[sample]=data
	return out


def load_srna_abundances(filename):
	print "loading from "+filename
	out = {}; contigs={}
	for i,line in enumerate(open(filename)):
		cut = line.strip().split("\t")
		if i==0:
			head=cut
			for sample in head[2:]:
				if "sRNA_ID" in sample:
					continue
				out[sample]={}
		else:
			for j,val in enumerate(cut):
				if j==0:
					srna=val
					contigs[val]=None
				elif j==1:
					if "Ga" not in val:
						srna+="-"+val
						out[srna]={}
				elif j==2 and "STRG" in val:
					srna+="-"+val
					out[srna]={}
				else:
					if val=="NA" or val=="#DIV/0!":
						val="0"
					sample=rename_sample(head[j])
					out[srna][sample]=float(val)
	return out, contigs


def rename_sample(sample):
	if "2017" not in sample:
		new_sample=sample[1:]
		if "1b" in sample or "1a" in sample or "9b" in sample:
			new_sample="0"+new_sample
		new_sample = new_sample[:-1] +"-" + new_sample[-1]
	else:
		cut = sample.split("-")[1].split("m")
		new_sample = cut[0]+"m"+"-"+cut[1]
	return sample[0]+new_sample

def load_taxa(filename):
	taxa={}
	for line in open(filename):
		cut=line.strip().split("\t")
		if len(cut)!=2:
			continue
		taxa[cut[0]]=cut[1].split(";")
	return taxa


def get_averages(contig_abundances, srna_abundances, year, taxa, labels):
	data = {}
	print "calculating averages for "+year
	for srna in srna_abundances:
		data[srna]=([],[])
		contig=srna.split("-")[0]
		for sample in srna_abundances[srna]:
			if year!="ALL":
				if year=="2016" and "2017" in sample:
					continue
				if year=="2017" and "am" not in sample and "pm" not in sample:
					continue
			if sample[0]=="g":
				continue
		
			contig_abund = contig_abundances[sample[1:]][contig]
			srna_abund = srna_abundances[srna][sample]
			data[srna][0].append(contig_abund)
			data[srna][1].append(srna_abund)
	out_contigs=[]
	out_srnas=[]
	colors=[]
	
	for srna in data:
		contig = srna.split("-")[0]
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




def plot_antisense(filename, year, taxa, labels, ax):
	print "plotting from "+filename
	x=[]; y=[]; colors=[]
	for i,line in enumerate(open(filename)):
		cut = line.strip().split("\t")
		if i==0:
			head=cut
		else:
			srnas=[]; genes=[]
			for j,val in enumerate(cut):
				if j==0:
					contig = val
					if contig in taxa:
						taxon=taxa[contig]
						if "Eukaryota" in taxon:
							c="w"
						if len(taxon)<2:
							c="black"
						elif taxon[1] in labels:
							c=labels[taxon[1]]
						else:
							c="black"
					else:
						c="black"
				elif j==1:
					continue
				elif j==2:
					continue
				else:
					if val=="NA" or val=="#DIV/0!":
						val="0"
					sample=head[j]
					if year!="ALL":
						if year=="2016" and ("am" in sample or "pm" in sample):
							continue
						if year=="2017" and "2017" not in sample:
							continue

					val = float(val)
					if sample[0]=="s":
						srnas.append(val)
					elif sample[0]=="g":
						genes.append(val)
			
			t_srnas=[]; t_genes=[]
                        for i,srna in enumerate(srnas):
                                gene=genes[i]
                                if gene<1 and srna<1:
                                        continue
                                t_srnas.append(srna)
                                t_genes.append(gene)

			srna = np.mean(t_srnas)
			gene = np.mean(t_genes)
			#print "\t".join(cut[:3] + [str(srna), str(gene), ";".join(taxon)])
			x.append(gene)
			y.append(srna)
			colors.append(c)

	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(0.001,100)
	ax.set_ylim(1, 5000)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.grid(linestyle='--', alpha=0.5)

	ax.set_ylabel("Mean sRNA expression (tpm)")
	ax.set_xlabel("Mean gene expression (tpm)")

	ax.set_title(year)
	ax.scatter(x, y, alpha=0.4, c=colors, s=8)



def plot_correlations(filename, year, taxa, labels, contig_abundances, ax):
	print "plotting from "+filename
	x=[]; y=[]; colors=[]
	for i,line in enumerate(open(filename)):
		cut = line.strip().split("\t")
		if i==0:
			head=cut
		else:
			srnas=[]; genes=[]
			for j,val in enumerate(cut):
				if j==0:
					contig = val
					if contig in taxa:
						taxon=taxa[contig]
						if "Eukaryota" in taxon:
							c="w"
						if len(taxon)<2:
							c="black"
						elif taxon[1] in labels:
							c=labels[taxon[1]]
						else:
							c="black"
					else:
						c="black"
				elif j==1:
					continue
				elif j==2:
					continue
				else:
					if val=="NA" or val=="#DIV/0!":
						val="0"
					sample=head[j]
					if year!="ALL":
						if year=="2016" and ("am" in sample or "pm" in sample):
							continue
						if year=="2017" and "2017" not in sample:
							continue

					val = float(val)
					alt_sample_name = rename_sample(sample)[1:]
					total_activity = contig_abundances[alt_sample_name][contig]
					#if total_activity<1:
					#	continue

					if sample[0]=="s":
						srnas.append(val)
					elif sample[0]=="g":
						genes.append(val)

			t_srnas=[]; t_genes=[]
			for i,srna in enumerate(srnas):
				gene=genes[i]
				if gene<1 and srna<1:
					continue
				t_srnas.append(srna)
				t_genes.append(gene)
			if len(t_genes)<6:
				continue
			
			if "Eukaryota" in taxon:
				continue
			
			test = stats.pearsonr(t_srnas, t_genes)
			if test[1]<0.01 and test[0]<0:
				c="b"
				print cut[2], contig, taxon, np.mean(srnas), np.mean(genes)
			else:
				c="lightgray"

			#if cut[2]=="STRG.142101" or cut[2]=="STRG.83785":
				#print cut[1], cut[2]
				#print_corr(srnas, genes)
			srna = np.mean(t_srnas)
			gene = np.mean(t_genes)
			x.append(srna)
			y.append(test[0])
			colors.append(c)

	ax.set_xscale('log')
	ax.set_xlim(1, 10000)
	ax.set_ylim(-1.1, 1.1)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.grid(linestyle='--', alpha=0.5)

	ax.set_ylabel("sRNA to gene Pearson correlation")
	ax.set_xlabel("Mean sRNA expression (tpm)")

	ax.set_title(year)
	ax.scatter(x, y, alpha=0.4, s=8, c=colors)


def print_corr(srnas, genes):
	print "\n"
	for i,srna in enumerate(srnas):
		gene = genes[i]
		#print str(srna)+"\t"+str(gene)
	print "\n"

######## MAIN ########
######################
######################
######################


antisense,contigs1 = load_srna_abundances("2016vs2017_all_asRNA_and_genes_non-normalized_TPM_table.txt")
intergenic,contigs2 = load_srna_abundances("2016vs2017_all_itsRNAs_non-normalized_TPM_table.txt")
taxa = load_taxa("contig_taxonomy.tab")


if "reload" in sys.argv:
	contig_abundances = load_contig_abundances(contigs1, contigs2)
	print "saving contig_abundances.pkl"
	output = open('contig_abundances.pkl', 'wb')
	pickle.dump(contig_abundances, output)
	output.close()
else:
	print "loading contig_abundances.pkl"
	pkl_file = open('contig_abundances.pkl', 'rb')
	contig_abundances = pickle.load(pkl_file)
	pkl_file.close()




##############  PLOTTING ###############
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)
fig, axs = plt.subplots(2, 2, figsize=(8,8))
labels={"Euryarchaeota":"magenta", "Bacteroidetes":"gold", "Cyanobacteria":"cyan", "Other":"black"}


##############  PLOTTING UPPER PANELS ###############
ax = axs[0,0]
plot_antisense("2016vs2017_all_asRNA_and_genes_normalized_TPM_table.txt", "ALL", taxa, labels, ax)
ax.set_title("Gene vs asRNA expression")

ax = axs[0,1]
plot_correlations("2016vs2017_all_asRNA_and_genes_normalized_TPM_table.txt", "ALL", taxa, labels, contig_abundances, ax)
ax.set_title("Gene vs asRNA correlation")

##############  PLOTTING LOWER PANELS ###############

for y,srna_type in enumerate(["intergenic", "antisense"]):
	print "plotting "+srna_type
	ax = axs[1,y]
	ax.loglog([0.001, 10000], [0.1, 1000000], '--', alpha=0.5)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(0.001, 1000)
	ax.set_ylim(5, 100000)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.grid(linestyle='--', alpha=0.5)

	ax.set_ylabel("sRNA expression (tpm)")
	ax.set_xlabel("Total contig activity (tpm)")
	if srna_type=="intergenic":
		contig_averages,srna_average,colors = get_averages(contig_abundances, intergenic, "ALL", taxa, labels)
		ax.set_title("itsRNA vs total activity")
	if srna_type=="antisense":
		contig_averages,srna_average,colors = get_averages(contig_abundances, antisense, "ALL", taxa, labels)
		ax.set_title("asRNA vs total activity")

	ax.scatter(contig_averages, srna_average, c=colors, alpha=0.4, s=8)


print "making legend..."
ax = fig.add_axes([0.0, 0.0, 1, 0.05])
ax.axis("off")
legend_elements = [Patch(facecolor=labels["Euryarchaeota"], edgecolor='k', label="Euryarchaeota", linewidth=1),
        Patch(facecolor=labels["Bacteroidetes"], edgecolor='k', label="Bacteroidetes", linewidth=1),
        Patch(facecolor=labels["Cyanobacteria"], edgecolor='k', label="Cyanobacteria", linewidth=1),
        Patch(facecolor=labels["Other"], edgecolor='k', label="Other", linewidth=1)]
ax.legend(handles=legend_elements, loc="lower center", framealpha=1, frameon=True, facecolor='w', ncol=4, columnspacing=1, handlelength=1, prop={'size': 12})


# annotations
axs[0,0].annotate("A", xy=(-0.19, 1.05), xycoords="axes fraction", fontsize=20)
axs[0,1].annotate("B", xy=(-0.19, 1.05), xycoords="axes fraction", fontsize=20)
axs[1,0].annotate("C", xy=(-0.19, 1.05), xycoords="axes fraction", fontsize=20)
axs[1,1].annotate("D", xy=(-0.19, 1.05), xycoords="axes fraction", fontsize=20)


plt.tight_layout(w_pad=0.5, h_pad=1)
plt.subplots_adjust(bottom=0.15)
plt.savefig("figure.png", dpi=600)






