#!/usr/bin/env python2
print "loading modules..."
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import operator

def load_tpms(filename):
	print"load TPM data from "+filename+"..."
	TPMs=[]
	for line in open(filename):
		if line[0]=="#":
			continue
		cut=line.strip().split("\t")
		if len(cut)!=9:
			continue
		info=cut[8].split(";")
		for field in info:
			if "TPM" in field:
				TPM=float(field.split('"')[1])
				break
		TPMs.append(TPM)

	print "sorting " +str(len(TPMs)) +" TPM values..."
	TPMs.sort(reverse=True)
	counts=[]
	for i in range(len(TPMs)):
		counts.append(i)
	return counts, TPMs


#MAIN
#font = {'family': 'arial', 'weight': 'normal', 'size': 12}
#plt.rc('font', **font)
fig = plt.figure(figsize=(6,5))
ax = fig.add_subplot(111)

count_maxs=[]
tpm_maxs=[]
	

counts,TPMs = load_tpms("small_antisense_ncRNAs.gff")
count_maxs.append(len(TPMs))
tpm_maxs.append(TPMs[0])
ax.scatter(counts, TPMs, s=10, alpha=1, label="asRNAs")


counts,TPMs = load_tpms("small_intergenic_ncRNAs.gff")
count_maxs.append(len(TPMs))
tpm_maxs.append(TPMs[0])
ax.scatter(counts, TPMs, s=10, alpha=1, label="itsRNAs")



print "plotting..."
if len(sys.argv)>3:
	# make legend
	handles, labels = ax.get_legend_handles_labels()
	ax.legend(handles[::-1], labels[::-1])
	hl = sorted(zip(handles, labels),
		key=operator.itemgetter(1))
	handles2, labels2 = zip(*hl)
	ax.legend(handles2, labels2)
else:
	ax.legend()

ax.set_xlabel("log10 Small ncRNA count", fontsize=15)
ax.set_ylabel("log10 Transcript expression (TPM)", fontsize=15)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1, max(count_maxs)*2)
ax.set_ylim(TPMs[-1], 3000)
plt.grid(b=True, which='both', color='k', alpha=0.1, linestyle='--')

plt.tight_layout()
plt.savefig("figure.png", dpi=300)
#plt.show()


