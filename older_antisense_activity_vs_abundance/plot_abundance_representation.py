#!/usr/bin/env python
#./plot_abundance_representation.py small_antisense_ncRNAs.gff all_meta-sRNAs_sig_correlation.txt
import sys
import matplotlib.pyplot as plt


tot_expr={}
for line in open("ALL.quant.counts"):
	if line[0]!="N":
		 continue
	cut=line.strip().split("\t")
	tot_expr[cut[0]]=float(cut[1])


taxa={}
for line in open("contig_taxonomy.tab"):
	cut=line.strip().split("\t")
	if len(cut)!=2:
		continue
	taxa[cut[0]]=cut[1].split(";")


sig_srnas={}
for line in open(sys.argv[2]):
	sig_srnas[line.strip()]=None




fig, ax = plt.subplots(figsize=(10,8))
ax.set_xscale('log')
ax.set_yscale('log')

labels={"Euryarchaeota":"blue", "Bacteroidetes":"magenta", "Cyanobacteria":"cyan", "Other":"black"}

data={}
sig_data={}

print "\t".join([  "srna", "srna_tpm", "contig_tpm", "taxon", "color" ])
for line in open(sys.argv[1]):
	cut=line.strip().split("\t")
	contig=cut[0]
	contig_tpm=tot_expr[contig]
	if contig in taxa:
		taxon=taxa[contig]
		if len(taxon)<2:
			label="Other"
		elif taxon[1] in labels:
			label=taxon[1]
		else:
			label="Other"
	else:
		label="Other"

	info=cut[8].split(";")
	for f in info:
		if "TPM" in f:
			srna_tpm=float(f.split('"')[1])
		elif "transcript_id" in f:
			srna = ".".join(f.split('"')[1].split(".")[:-1])
			if srna in sig_srnas:
				sig=True
			else:
				sig=False


	if label not in data:
		data[label]=([], [])
		sig_data[label]=([], [])

	if sig==True:
		sig_data[label][0].append(contig_tpm)
		sig_data[label][1].append(srna_tpm)
	else:
		data[label][0].append(contig_tpm)
		data[label][1].append(srna_tpm)
	
	#print "\t".join([  srna, str(srna_tpm), str(contig_tpm), label, labels[label]  ])
	if srna_tpm>contig_tpm*1 and label=="Euryarchaeota":
		print srna


for label in data:
	ax.scatter(data[label][0], data[label][1], label=label, c=labels[label], alpha=0.1)
for label in sig_data:
	ax.scatter(sig_data[label][0], sig_data[label][1], c=labels[label], alpha=1)


ax.legend()
ax.set_ylabel("sRNA expression (tpm)")
ax.set_xlabel("Total contig expression (tpm)")
plt.show()
