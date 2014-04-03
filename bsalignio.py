from collections import namedtuple
from itertools import imap, izip
from operator import attrgetter

import gzip

from ectools.nucio import lineRecordIterator, fileIterator

BSBaseSample = namedtuple('BSBaseSample', ["ref","off",
                               "strand", "mstr",
                               "mcy","ustr",
                               "ucy", "filt_cycle",
                               "filt_readline",
                               "filt_allele","filt_mapq",
                               "filt_baseq"])
BSBaseSampleTypes = [str, int, str, str, int ,str, int , int ,int ,int ,int ,int]


BSPosition = namedtuple("BSPosition", ["pos","samples"])


def baseCoverage(basesample):
    return basesample.mcy + basesample.ucy


def methValue(basesample):
    '''Calculate the methvalue from a BSBaseSample'''
    total = basesample.mcy + basesample.ucy
    if not total > 0:
        return 0
    return float(basesample.mcy) / total


def getPositions(files):
    '''Takes a list of input files and returns
    a list of BSPosition's
    
    '''
    #func to filter header out
    head_filt = lambda line : False if line.startswith("ref") else True
    
    record_it = lambda fh : lineRecordIterator(fh, BSBaseSample, BSBaseSampleTypes,
                                               head_filt, "\t")

    opener = lambda fn : gzip.open if fn.endswith(".gz") else open
    
    #read all data into a matrix
    observations = map(lambda fn : list(fileIterator(fn, record_it, opener(fn))), files)

    return map(BSPosition._make,
               izip(imap( attrgetter("off"), observations[0]), izip(*observations)))
