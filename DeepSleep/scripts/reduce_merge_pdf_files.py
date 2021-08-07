import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')

from glob import glob
import os
import re
import uproot
import awkward
import numpy as np
from modules.AnaDict import AnaDict

'''
Script to get only necessary branches for
pdfweight and event matching
'''

#YEARS   = ['2016','2017','2018']
YEARS = ['2018']
PDF_DIR = '/cms/data/store/user/bcaraway/NanoAODv7/PDF/'
KEEP_VARS = ['run','luminosityBlock','event','LHEPdfWeight', 'LHE_HT']
OUT_DIR = '/cms/data/store/user/ttxeft/NanoAODv7/PDF/'

def main():
    for y in YEARS:
        sample_dirs = glob(f'{PDF_DIR}/{y}/*')
        for sample_dir in sample_dirs:
            process_dir(sample_dir)

def process_dir(s_dir):
    roofiles = glob(f'{s_dir}/**/*.root', recursive=True)
    # open root files one by one
    metaData = {'lhapdf':''} # store LHA ID here
    pdfData = {} # store vars here
    for i, roofile in enumerate(roofiles):
        with uproot.open(roofile) as roo:
            tree = roo.get('Events')
            # retieve info and concatenate into a few arrays
            if i != 0:
                pdfData = {var: aj_concatenate([pdfData[var],tree.array(var)]) for var in KEEP_VARS}
            else:
                metaData['lhapdf'] = re.search(r'LHA IDs \d+ - \d+', tree['LHEPdfWeight'].title.decode()).group()
                pdfData = {var: tree.array(var) for var in KEEP_VARS}
    # store arrays into 1 file per sample
    out_dict = AnaDict({**pdfData,**metaData})
    # save in '/cms/data/store/user/ttxeft/NanoAODv7/PDF/'
    out_dict.to_pickle(f'{OUT_DIR}/'+s_dir.split('/PDF/')[-1]+'.pkl')


def aj_concatenate(arrays):
    try: # assume we are handling awkward arrays
        contents = np.concatenate([j.flatten() for j in arrays])
        counts = np.concatenate([j.counts for j in arrays])
        return awkward.JaggedArray.fromcounts(counts, contents)
    except AttributeError: # catch ndarray objects
        return np.concatenate([arr for arr in arrays])
        

if __name__ == '__main__':
    main()
