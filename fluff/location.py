import HTSeq
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from fluff.plot import create_grid_figure

def get_location_data(interval, tracks, fragmentsize=200):
	chrom,start,end = interval
	window = HTSeq.GenomicInterval(chrom, start, end, "+")
	profiles = []
	for track in tracks:
		bamfile = HTSeq.BAM_Reader(track)
		profile = zeros(end - start, dtype="i")
		for almnt in bamfile[window]:
			almnt.iv.length = fragmentsize
			profile[almnt.iv.start - start:almnt.iv.end - start] += 1	
		profiles.append(profile)
	return profiles	

tracks = [
	"/home/simon/prj/xenopus/xenTro3b/bam/H3K27me3_stage9.bam",
	"/home/simon/prj/xenopus/xenTro3b/bam/H3K27me3_stage12.bam",
	"/home/simon/prj/xenopus/xenTro3b/bam/H3K4me3_stage12.bam",
	"/home/simon/prj/xenopus/xenTro3b/bam/H3K4me1_stage9.bam",
]


chrom = "scaffold_9"
start = 20787459
end = 20811523
profiles = get_location_data((chrom, start, end), tracks)

genestart = 20790000
geneend = 20800000
exonstarts = [0,9000]
exonsizes = [1000,1000]
genestrand = "-"

fig, axes = create_grid_figure(len(tracks) + 1, 1, plotwidth=5, plotheight=0.2)
groups = {0:[0,1]}
colors = ["#e41a1c","#e41a1c","#4daf4a","#377eb8"]
for i, profile in enumerate(profiles):
	ax = axes[i][0]
	#ax.plot(range(len(profiles[i])), profiles[i])
	ax.fill_between(range(start, end), zeros(len(profile)), profile, edgecolor='face', facecolor=colors[i], linewidth=0)
	if i == 2:
		ax.fill_between(range(start, end), zeros(len(profiles[i+ 1])), profiles[i + 1], edgecolor='face', facecolor=colors[i + 1], alpha=0.5, linewidth=0)
		
	for loc, spine in ax.spines.iteritems():
		spine.set_color('none')

for group, indices in groups.items():
	group_ymax = 0
	for i in indices:
		ymin, ymax = axes[i][0].get_ylim()
		if ymax > group_ymax:
			group_ymax = ymax
	for i in indices:
		axes[i][0].set_ylim(0, group_ymax)

ax = axes[len(tracks)][0]
for loc, spine in ax.spines.iteritems():
	spine.set_color('none')
ax.set_ylim(-2,2)
ax.axhline(0, (genestart - start)/float(end - start), (geneend - start)/ float(end - start), color="black")
step = fig.get_figwidth() / 400 
#print start,end

for i in arange((genestart -start)/ float(end-start) + step, (geneend - start)/ float(end-start) - step, step):
	#print i , i + (step /2 )
	#arr = plt.Arrow(i, 0, step, 0, width=1, aa=False)
	if genestrand == "+":
		arr = FancyArrowPatch((i,0),(i + step,0),arrowstyle='->',mutation_scale=30)
	else:
		arr = FancyArrowPatch((i + step,0),(i,0),arrowstyle='->',mutation_scale=30)
		
	ax.add_patch(arr)
	
	#arr.set_facecolor('g')

for exonstart,exonsize in zip(exonstarts, exonsizes):
	ax.axhspan(-0.8,0.8, (genestart + exonstart - start)/float(end - start), (genestart + exonstart + exonsize - start)/ float(end - start), color="black")

fig.savefig("test.png", dpi=300)
