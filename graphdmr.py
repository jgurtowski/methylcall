#!/usr/bin/env python

import sys
from collections import namedtuple
from operator import attrgetter
from itertools import repeat, islice

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.path.append("/sonas-hs/schatz/hpc/home/gurtowsk/workspace")
sys.path.append(".")
from ectools.nucio import fileIterator, lineRecordIterator, lineItemIterator
from ectools.misc import iterApply
from ectools.args import CLArgument, parseArgs, getHelpStr, argflag, arglist
from ectools.log import logger
from ectools.algorithm import binarySearch, expandRegion
from bsalignio import getPositions, methValue, baseCoverage


Region = namedtuple('Region', ["chr","start","end"])
RegionTypes = [str, int, int]

description = ("Usage: graphdmr.py [options] regions.txt bsalign1 bsalign2 [bsalign3 ..]\n\n"
               "Graph DMRs from a region file format chr[TAB]start[TAB]end\n\n")

arglist = [["mincov", "mincov", int ,4,
            ("Minimum coverage that must be present in all samples for a "
             "particular site to be included in the graph (Default: 4)" )],
           ["log", "log", argflag, False, "Turn on some logging"],
           ["labels", "labels", arglist, [], 
            ("Labels for methyl libraries. Should be in the saame "
             "order input files")],
           ["pout", "pout", str, "out",
            "prefix for output file (Default: out)"]]

cl_arglist = map(CLArgument._make, arglist)
(arg_map, remaining_args) = parseArgs(sys.argv[1:], cl_arglist)

if not len(sys.argv) >= 4:
    sys.exit(getHelpStr(description, cl_arglist) + "\n")

log = logger(sys.stderr if arg_map["log"] else None)

regionfile = remaining_args[0]
methylfiles = remaining_args[1:]

log("Reading Positions from Files")
positions = getPositions(methylfiles)

pp = PdfPages(arg_map["pout"] + ".pdf")

regmk = lambda arr : Region._make(map(lambda t,v : t(v), RegionTypes, arr[:3]))
regit = lambda fh: iterApply(regmk, lineItemIterator(fh))

for region in fileIterator(regionfile, regit):
    if not region.chr == positions[0].samples[0].ref:
        raise Exception, "Regions are for a different Chromosome"
    
    log("Working on %s:%d-%d" % region)
    def searchexp(item):
        if item > region.end:
            return -1
        elif item < region.start:
            return 1
        return 0

    exp = lambda pos : pos >= region.start and pos <= region.end
    fexp = lambda f : lambda p : f(p.pos)
    
    rstart,rend = expandRegion(positions,fexp(exp),
                               binarySearch(positions, fexp(searchexp)))

    filt_pos = islice(positions,rstart,rend+1)
    
    #filt_pos = filter(lambda p: exp(p.pos), positions)

    covcheck = lambda c : c >= arg_map["mincov"]
    cov_filt_pos = filter(lambda position: all(map(covcheck,map(baseCoverage,position.samples))),
                          filt_pos)
    
    pos_by_sample = zip(*map(attrgetter("samples"), cov_filt_pos))
    methyl_by_sample = map(lambda bases : map(methValue,bases),
                           pos_by_sample)
    
    
    posints = map(lambda cp: cp.pos - cov_filt_pos[0].pos, cov_filt_pos)
    pdata = zip(repeat(posints), methyl_by_sample)

    plt.ylim([0,1.2])
    map(lambda d : plt.plot(*d), pdata)

    plt.title("%s: %d-%d" % region)
    plt.legend(arg_map["labels"])

    plt.savefig(pp, format="pdf")
    plt.clf()

pp.close()


