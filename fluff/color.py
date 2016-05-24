# Define colors and color utilities for plots
#
# Copyright (c) 2012-2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This script is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License

import colorbrewer
from matplotlib.colors import colorConverter, LinearSegmentedColormap, cnames

COLOR_MAP = {
    "red": "#e41a1c",
    "blue": "#377eb8",
    "green": "#4daf4a",
    "purple": "#984ea3",
    "orange": "#ff7f00",
    "yellow": "#ebbc57",
    "brown": "#a65628",
    "pink": "#f781bf",
    "grey": "#999999",
    "white": "#ffffff",
    "black": "#000000",
}


def is_pal(name):
    return name in colorbrewer.__dict__


def get_pal(name, n=None):
    ns = colorbrewer.__dict__[name].keys()
    if n > max(ns) or not n:
        n_index = max(ns)
    elif n < min(ns):
        n_index = min(ns)
    else:
        n_index = n

    pal = colorbrewer.__dict__[name][n_index][:n]
    for i in range(len(pal)):
        pal[i] = [x / 255.0 for x in pal[i]]

    return pal


DEFAULT_COLORS = get_pal("Set1")


def parse_colors(colors):
    if type("") == type(colors):
        colors = [x.strip() for x in colors.split(",")]

    parsed = []
    for c in colors:
        if type("") == type(c):
            # Named color
            if c in COLOR_MAP:
                parsed.append(COLOR_MAP[c])
            elif c in cnames:
                parsed.append(cnames[c])
            elif is_pal(c):
                # c is a Colorbrewer palette name
                parsed += get_pal(c)
            elif len(c.split(":")) == 2:
                p, n = c.split(":")
                if is_pal(p):
                    parsed += get_pal(p, int(n))
                else:
                    raise ValueError("%s is not a valid colorbrewer name!" % p)
            # Hex code?
            else:
                if c.startswith("#"):
                    parsed.append(c)
                else:
                    try:
                        int(c, 16)
                        parsed.append("#" + c)
                    except:
                        raise ValueError("Unknown color definition %s" % c)
        else:
            # c is not a strint, assume it's already a valid color
            parsed.append(c)
    return parsed

def create_colormap(*args):
    col = [colorConverter.to_rgb(c) for c in args]
    
    step = 1.0 / (len(col) - 1)
    
    cdict = {
            'red': [[i * step, col[i][0], col[i][0]] for i in range(len(col))],
            'green': [[i * step, col[i][1], col[i][1]] for i in range(len(col))],
            'blue': [[i * step, col[i][2], col[i][2]] for i in range(len(col))]
            }
    
    return LinearSegmentedColormap('custom', cdict, 256)
