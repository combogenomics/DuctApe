#!/usr/bin/env python
"""
utils

Common library

Spare parts
"""

# Borrowed from: www.garyrobinson.net
def slice_it(li, cols=10):
    start = 0
    for i in xrange(cols):
        stop = start + len(li[i::cols])
        yield li[start:stop]
        start = stop
        
def get_span(li, span=4):
    i = 0
    while i < len(li):
        yield li[i : i+span]
        i += span