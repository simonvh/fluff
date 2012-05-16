import HTSeq
import sys
from numpy import zeros,min,max

def get_profile(interval, track, fragmentsize=200):
	chrom,start,end = interval
	window = HTSeq.GenomicInterval(chrom, start, end, "+")
	bamfile = HTSeq.BAM_Reader(track)
	profile = zeros(end - start, dtype="i")
	for almnt in bamfile[window]:
		almnt.iv.length = fragmentsize
		profile[almnt.iv.start - start:almnt.iv.end - start] += 1
	return profile

def load_profile(interval, tracks, fragmentsize=200):
	profiles = []
	for track_group in tracks:
		if type(track_group) == type([]):
			profile_group = []
			for track in track_group:
				profile = get_profile(interval, track, fragmentsize)
				profile_group.append(profile)
		else:
				track = track_group
				profile_group = get_profile(interval, track, fragmentsize)
		profiles.append(profile_group)
	return profiles

def get_free_track(overlap, start, end, max_end):
	for i,track in enumerate(overlap):
		if max(track[start:end]) == 0:
			track[start:end] += 1
			return overlap, i
	overlap.append(zeros(max_end, dtype="i"))
	overlap[-1][start:end] += 1
	return overlap, len(overlap) - 1

def load_annotation(interval, fname):
	genes = []
	chrom, start, end = interval
	for line in open(fname):
		if not line.startswith("#") and not line.startswith("track"):
			vals = line.strip().split("\t")
			if len(vals) != 12:
				sys.stderr.write("Need BED 12 format for annotation\n")
				sys.exit(1)
			for i in [1,2,6,7]:
				vals[i] = int(vals[i])
			if vals[0] == chrom:
				if (vals[1] > start and vals[1] < end) or (vals[2] > start and vals[2] < end):
					sys.stderr.write("Adding %s\n" % vals[3])
					genes.append(vals)
	
	if len(genes) == 0:
		return []

	min_start = min([gene[1] for gene in genes])
	max_end = max([gene[2] for gene in genes])
	overlap = []
	gene_tracks = {}
	for gene in genes:
		overlap,i = get_free_track(overlap, gene[1]  - min_start, gene[2] - min_start, max_end - min_start)	
		if gene_tracks.has_key(i):
			gene_tracks[i].append(gene)
		else:
			gene_tracks[i] = [gene]
	return gene_tracks
