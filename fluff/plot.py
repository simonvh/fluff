from numpy import *
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter,NullLocator
from matplotlib.font_manager import fontManager, FontProperties
from matplotlib.patches import FancyArrowPatch
from fluffio import *
import sys
#from scipy.stats import scoreatpercentile

DEFAULT_COLORS = ["#e41a1c","#4daf4a","#377eb8"]
FONTSIZE = 8

def coverage_plot(ax, x, data, color="red"):
	"""
	ax = matplotlib axes instance
	x = x-axis coordinates
	data = profile data
	color = color in any way matplotlib accepts
	"""
	
	# Might change this into an argument for the function
	percs = [10, 25]
	alphas = [0.1, 0.4]

	# Convert to numpy array
	vals = array(data)

	# Get the median
	m = median(vals, axis=0)

	# Draw the minimum percentiles
	lines = [array([scoreatpercentile(vals[:,i], perc) for i in range(len(vals[0]))]) for perc in percs] + [m]
	for (line_min, line_max), alpha in zip([(lines[i], lines[i + 1]) for i in range(len(percs))], alphas):
		ax.fill_between(x, line_min, line_max, facecolor=color, alpha=alpha, edgecolor='face')	

	# Draw the maximum percentiles
	lines = [m] + [array([scoreatpercentile(vals[:,i], 100 - perc) for i in range(len(vals[0]))]) for perc in percs[::-1]] 
	for (line_min, line_max), alpha in zip([(lines[i], lines[i + 1]) for i in range(len(percs))], alphas[::-1]):

		ax.fill_between(x, line_min, line_max, facecolor=color, alpha=alpha, edgecolor='face')	
		
		# Draw the median
		ax.plot(x, m, color="black", alpha=0.95, linewidth=0.8)

def create_grid_figure(nrows, ncolumns, plotwidth=2.0, plotheight=2.0, pad=0.1, padleft=0.1, padright=0.1, padtop=0.1, padbottom=0.1, clean=True):
	wsize = padleft + (ncolumns * plotwidth) + (pad * (ncolumns - 1)) + padright
	hsize = padtop + (nrows * plotheight) + (pad * (nrows - 1)) + padbottom

	fig = plt.figure(figsize=(wsize, hsize))
	wpadfraction = pad / wsize
	hpadfraction = pad / hsize
	wplotsize = plotwidth / wsize
	hplotsize = plotheight / hsize 

	axes = {}
	# Create all the subplots
	for row in range(nrows):
		axes[row] = {}
		for col in range(ncolumns):
			axes[row][col] = plt.subplot(nrows, ncolumns, row * ncolumns + col + 1)
			
			# No labels, ticks, etc.
			if clean:
				for ax in [axes[row][col].xaxis, axes[row][col].yaxis]:
					ax.set_major_formatter(NullFormatter())
					ax.set_major_locator(NullLocator())
	
	# Resize all the subplots
	for row in range(nrows):
		for col in range(ncolumns):
			x0 = (padleft / wsize ) + (wplotsize + wpadfraction) * col
			x1 = wplotsize
			y0 = (padbottom / hsize) + (nrows - row - 1) * (hplotsize + hpadfraction) 
			y1 = hplotsize
			coords = [x0, y0, x1, y1]
			axes[row][col].set_position(coords)

			for s in axes[row][col].spines.values():
				s.set_linewidth(0.8)


	return fig, axes

def profile_screenshot(fname, intervals, tracks, colors=None, scalegroups=[], scale=True, annotation=None, bgmode="color"):
	# Colors
	if not colors:
		colors = DEFAULT_COLORS

	# Sizes
	plotwidth = 5
	plotheight = 0.3
	padh = 0
	padw = 0.1
	padleft = 0.1
	padright = 0.1
	padtop = 0.1
	padbottom = 0.1
	clean = True

	# Genomic scale
	scale_height = 0.1
	# Annotation track height
	
	annotation_height = 0.0
	gene_tracks = []
	if annotation:
		for interval in intervals:
			gene_tracks.append(load_annotation(interval, annotation))
		max_tracks =  max([len(x.keys()) for x in gene_tracks])
		annotation_height = 0.2 * max_tracks

	#
	ncolumns = len(intervals)
	nrows = len(tracks)

	wsize = padleft + (ncolumns * plotwidth) + (padw * (ncolumns - 1)) + padright
	# Profile plots
	hsize = padtop + (nrows * plotheight) + (padh * (nrows - 1)) + padbottom
	hsize += scale_height + padh + annotation_height + padh

	fig = plt.figure(figsize=(wsize, hsize))
	padtopfraction = padtop / hsize
	padbottomfraction = padbottom / hsize
	wpadfraction = padw / wsize
	hpadfraction = padh / hsize
	wplotsize = plotwidth / wsize
	hplotsize = plotheight / hsize 
	hannotation = annotation_height / hsize
	hscale = scale_height / hsize
	
	for int_num, interval in enumerate(intervals):
		all_axes = []
		
		# Scale ax
		ax = plt.axes([padleft / wsize + (wplotsize + wpadfraction) * int_num, 1 - hscale - padtopfraction, wplotsize, hscale], axisbg="white")
		all_axes.append(ax)
		
		# Profile axes
		for i in range(nrows):
			bgcol = "white"
			if bgmode == "stripes":
				bgcol = {0:"white",1:(0.95,0.95,0.95)}[i % 2]
			elif bgmode == "color":
				bgcol = colors[i % len(colors)]

			ax = plt.axes([padleft / wsize + (wplotsize  + wpadfraction) * int_num, 1 - hscale - padtopfraction - (i + 1) * (hplotsize + hpadfraction), wplotsize, hplotsize], axisbg=bgcol)
			if bgmode == "color":
				ax.patch.set_alpha(0.07)
			all_axes.append(ax)

		# Annotation axes
		ax = plt.axes([padleft / wsize + (wplotsize + wpadfraction) * int_num, 0 + padbottomfraction, wplotsize, hannotation], axisbg="white")
		all_axes.append(ax)
	
		# No labels, ticks, etc.
		for axes in all_axes[1:]:
			for ax in [axes.xaxis, axes.yaxis]:
				ax.set_major_formatter(NullFormatter())
				ax.set_major_locator(NullLocator())
		
		for axes in all_axes:
			for s in axes.spines.values():	
				s.set_color('none')
		
		chrom,start,end = interval
		
		# Format the genomic scale
		ax = all_axes[0]
		ax.yaxis.set_major_formatter(NullFormatter())
		ax.xaxis.set_major_formatter(NullFormatter())
		ax.yaxis.set_major_locator(NullLocator())
		ax.set_xlim(start, end)
		font = FontProperties(size=FONTSIZE / 1.25, family=["Nimbus Sans L", "Helvetica", "sans-serif"])
		for s in [s for s in ax.xaxis.get_ticklocs()[:-1] if s > start and s < end][:-1]:
			plt.text((s - start) / (end - start) + 0.01, 0.5, str(int(s)), horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, fontproperties=font)

		# Load the actual data
		profiles = load_profile(interval, tracks)

		# Plot the profiles
		color_index = 0
		track_maxes = []
		for i, profile_group in enumerate(profiles):
			ax = all_axes[i + 1]
			
			if type(profile_group) == type([]):
				maxes = []
				for profile in profile_group:
					if len(profile_group) > 1:
						alpha = 0.4
					else:
						alpha = 1
					ax.fill_between(range(start, end), zeros(len(profile)), profile, edgecolor='face', facecolor=colors[color_index % len(colors)], linewidth=0, alpha=alpha)
					color_index += 1
					maxes.append(max(profile) * 1.1)
				track_maxes.append(max(maxes))
				ax.set_ylim(0,max(maxes))
			else:
				# Single track
				profile = profile_group
				ax.fill_between(range(start, end), zeros(len(profile)), profile, edgecolor='face', facecolor=colors[color_index % len(colors)], linewidth=0)
				color_index += 1
				track_maxes.append(max(profile) * 1.1)
				ax.set_ylim(0, max(profile) * 1.1)
		
		if scalegroups and len(scalegroups) > 0:
			for group in scalegroups:
				for ax_id in group:
					all_axes[ax_id].set_ylim(0, max([track_maxes[i - 1] for i in group]))

		# Plot the gene annotation
		if annotation:
			ax = all_axes[-1]
			ax.set_ylim(- 4 * max_tracks, 0)
			for track_id, genes in gene_tracks[int_num].items():
				for gene in genes:
					h_gene = -4 * track_id - 2
				
					genestart = gene[1]
					geneend = gene[2]
					exonstarts = [int(x) for x in gene[11].split(",") if x]
					exonsizes =  [int(x) for x in gene[10].split(",") if x]
					genestrand = gene[5]
		
					ax.axhline(h_gene, (genestart - start)/float(end - start), (geneend - start)/ float(end - start), color="black")
					step = plotwidth/ 400.0
					for i in arange((genestart -start)/ float(end-start) + step, (geneend - start)/ float(end-start) - step, step):
						if genestrand == "+":
							arr = FancyArrowPatch((i,h_gene),(i + step,h_gene),arrowstyle='->',mutation_scale=5)
						else:
							arr = FancyArrowPatch((i + step,h_gene),(i,h_gene),arrowstyle='->',mutation_scale=5)
						ax.add_patch(arr)
				
					for exonstart,exonsize in zip(exonstarts, exonsizes):
						ax.axhspan(h_gene - 0.6, h_gene + 0.6, (genestart + exonstart - start)/float(end - start), (genestart + exonstart + exonsize - start)/ float(end - start), color="black")

	
	
	plt.savefig(fname, dpi=300)

if __name__ == "__main__":
	
#	import matplotlib.pyplot as plt
#	
#	data = [
#		[10, 7, 5, 7, 10],
#		[9, 7, 5, 7, 9],
#		[8, 6, 5, 6, 8],
#		[7, 6, 5, 6, 7],
#		[6, 5, 5, 5, 6],
#		[5, 5, 5, 5, 5],
#		[4, 4, 5, 4, 4],
#		[3, 4, 5, 4, 3],
#		[2, 3, 5, 3, 2],
#		[1, 3, 5, 3, 1],
#	]
#	x = arange(5)	
#	plt.figure()
#	
#	ax = plt.subplot(111)
#	
#	coverage_plot(ax, x, data)
#	plt.show()

	tracks = [
		"/home/simon/prj/xenopus/xenTro3b/bam/EZH2_stage9.bam",
		"/home/simon/prj/xenopus/xenTro3b/bam/Jarid2_stage9.bam",
		"/home/simon/prj/xenopus/xenTro3b/bam/H3K27me3_stage9.bam",
		"/home/simon/prj/xenopus/xenTro3b/bam/H3K27me3_stage12.bam",
		"/home/simon/prj/xenopus/xenTro3b/bam/H3K27me3_stage16.bam",
		"/home/simon/prj/xenopus/xenTro3b/bam/H3K27me3_stage30.bam",
		"/home/simon/prj/xenopus/xenTro3b/bam/H3K4me1_stage9.bam",
		"/home/simon/prj/xenopus/xenTro3b/bam/H3K4me3_stage9.bam",
		"/home/simon/prj/xenopus/xenTro3b/bam/H3K4me3_stage12.bam",
	]
	intervals = [
		["scaffold_3", 112047153, 112091309],
		["scaffold_1", 20081849, 20104636],	
		["scaffold_1", 126564001, 126606826],
#		["scaffold_9", 20687459, 20911523],
	]

	
	colors = ["#a65628", "#ff7f00", "#e41a1c","#e41a1c","#e41a1c", "#e41a1c", "#377eb8", "#4daf4a", "#4daf4a"]
	annotation = "/home/simon/prj/xenopus/xenTro3b/annotation/Xentr7_2_Stable_name.bed"
	scalegroups = [[1],[2],[3,4,5,6],[7],[8,9]]
	profile_screenshot("test.png", intervals, tracks, colors=colors, annotation=annotation, scalegroups=scalegroups, bgmode="color")
