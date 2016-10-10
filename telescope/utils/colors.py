# -*- coding: utf-8 -*-
__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"


def c2str(rgb):
    """ Represent a color tuple as a string

    Args:
        rgb (:obj:`tuple` of int): Color representation as (R, G, B), where
            each color is an integer between 0-255.

    Returns:
        string: Comma separated string with values for each channel

    """
    return '%d,%d,%d' % rgb


DARK2_PALETTE = {
    'teal':       ( 27, 158, 119), #1b9e77
    'vermilion':  (217,  95,   2), #d95f02
    'purple':     (117, 112, 179), #757063
    'magenta':    (231,  41, 138), #e7298a
    'green':      (102, 166,  30), #66a61e
    'yellow':     (230, 171,   2), #e6ab02
    'brown':      (166, 118,  29), #a6761d
    'gray':       (102, 102, 102), #666666
}


GREENS = [ #(50,168,133),
           #(73,177,146),
           #(95,187,160),
           (118,197,173),
           #(141,206,187),
           (164,216,201),
           #(189,226,214),
           (209,236,228),
           (232,245,241),
           #(255,255,255),
]
