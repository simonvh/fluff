# Define colors and color utilities for plots
#
# Copyright (c) 2012-2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This script is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License

from matplotlib.colors import colorConverter, LinearSegmentedColormap

COLOR_MAP = {"red":"#e41a1c","blue":"#377eb8","green":"#4daf4a","orange":"#ff7f00","brown":"#a65628","purple":"#984ea3","yellow":"#ffff33"}
DEFAULT_COLORS = "#e41a1c,#377eb8,#4daf4a,#984ea3,#ff7f00,#ffff33,#a65628"

def create_colormap(col1, col2):
	c1 = colorConverter.to_rgb(col1)
	c2 = colorConverter.to_rgb(col2)
	
	cdict = {
		'red': ((0.,c1[0], c1[0]),(1.,c2[0], c2[0])),
		'green': ((0.,c1[1], c1[1]),(1.,c2[1], c2[1])),
		'blue': ((0.,c1[2], c1[2]),(1.,c2[2], c2[2]))
	}
	return LinearSegmentedColormap('custom', cdict, 256)


