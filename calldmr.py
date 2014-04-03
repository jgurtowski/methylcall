#!/usr/bin/env python

import sys
from operator import attrgetter, itemgetter
from itertools import izip, islice, count
from collections import namedtuple

from scipy.stats.mstats import kruskalwallis

sys.path.append("/sonas-hs/schatz/hpc/home/gurtowsk/workspace")
sys.path.append(".")
sys.path.append("../jamesdmr")

from ectools.misc import takeFirstWhile
from ectools.args import CLArgument, parseArgs, getHelpStr, argflag
from ectools.log import logger
from ectools.nucio import FileOrStream

from bsalignio import getPositions, baseCoverage, methValue


BSWindow = namedtuple("BSWindow", ["start","end", "positions"])
KWResult = namedtuple("KWResult", ["p","h", "window"])

def BSWindowGen(positions, window_size=1000, step_size=500, start=0):
    '''Takes a list of 'BSPosition's and gives subsets
    window_size and step_size are in base pairs

    start_opt is for internal optimization so we don't need to search the whole
    list every time

    positions need to be sorted 
    '''
    last_base = positions[-1].pos
    exp = lambda p : p.pos > start and p.pos <= win_end
    stop = lambda p : p.pos > win_end
    start_opt = 0
    while True:
        win_end = start + window_size
        if(win_end > last_base):
            return
        good = takeFirstWhile(lambda ps: exp(ps[1]), lambda ps: stop(ps[1]),
                              izip(count(start_opt), islice(positions,start_opt,None)))
        if bool(good):
            start_opt = good[0][0]
        
        yield BSWindow(start, win_end, map(itemgetter(1),good))
        start += step_size

description = ("Usage: calldmr.py [options] bsalign1 bsalign2 [bsalign3...]\n\n" 
               "Call DMRs from BsmoothAlign output using Kruskal-Wallis test\n\n")

arglist = [["windowsize", "windowsize", int, 500,
            ("Size of windows for evaluating DMRs (Default: 500)")],
           ["stepsize","stepsize", int, 250,
            ("Size of step between windows (Default: 250)")],
           ["mincov", "mincov", int, 4,
            ("Minimum coverage that must be present in all samples for a "
             "particular base to be included in the analysis (Default: 4)" )],
           ["minwinsites", "minwinsites",int, 20,
            ("Minimum number of sites (C's, CpG's ...) in a window (Default: 20)")],
           ["q", "q", float, 0.05,"q value for FDR. (Default: 0.05)"],
           ["log","log",argflag, False, "Turn on some logging"],
           ["out", "out", str, "stdout", "Output file"]]

cl_arglist = map(CLArgument._make, arglist)
(arg_map,remaining_args) = parseArgs(sys.argv[1:], cl_arglist)

if not len(remaining_args) >= 3:
    sys.exit(getHelpStr(description,cl_arglist) + "\n")

log  = logger(sys.stderr if arg_map["log"] else None)

ofh = sys.stdout if arg_map["out"] == "stdout" else open(arg_map["out"],"w")

log("Reading Positions from Files")
positions = getPositions(remaining_args)
chrom = positions[0].samples[0].ref


windowsize = arg_map["windowsize"]
stepsize = arg_map["stepsize"]
log("Computing Kruskal-Wallis test on windows of size %d and step %d" % (windowsize,stepsize))
testresults = []
for window in BSWindowGen(positions, windowsize, stepsize):
    
    covcheck = lambda c : c >= arg_map["mincov"]
    filt_pos = filter(lambda position: all(map(covcheck,map(baseCoverage,position.samples))),
                      window.positions)
    
    if not len(filt_pos) >= arg_map["minwinsites"]:
        continue
    
    pos_by_sample = zip(*map(attrgetter("samples"), filt_pos))
    methyl_by_sample = map(lambda bases : map(methValue,bases),
                           pos_by_sample)
    
    try:
        (h,p) = kruskalwallis(*methyl_by_sample)
    except Exception as e:
        sys.stderr.write("Error: %s \n" % e)
        continue

    testresults.append(KWResult(p, h, window))

log("Sorting Results")    
testresults.sort()

log("Writing Output")
m = float(len(testresults))
q = arg_map["q"]
k = 0
 
for result in testresults:
    if result.p > ((k+1) / m * q):
        break
    ofh.write("\t".join(map(str,[chrom,result.window.start,
                                 result.window.end,result.h,result.p])))
    ofh.write("\n")
    k+=1
    
ofh.close()






