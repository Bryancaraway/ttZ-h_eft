#from ROOT import TFile                                                                                         
import uproot
import sys
import os
import re
import concurrent.futures
import multiprocessing
import numpy as np
#executor = concurrent.futures.ThreadPoolExecutor()

def worker(roo):
    if 'kodiak' in roo:
        roo = roo.replace('root://kodiak-se.baylor.edu//','/cms/data/')
    t=uproot.open(roo)['Events']

    gw    = t.array('genWeight'  )
    ttbb  = t.array('genTtbarId')
    ttbb  = ttbb % 100
    scale = t.array('LHEScaleWeight').pad(9).fillna(1)
    ps    = t.array('PSWeight'      ).pad(4).fillna(1)
    # need pdf as well
    pdf_up   = t.array('pdfWeight_Up') * np.sign(gw)
    pdf_down = t.array('pdfWeight_Down') * np.sign(gw)
    tot_count= sum(gw>=0.0) - sum(gw<0.0)
    #
    ttbb_count = sum((gw>=0.0) & (ttbb>=51)) - sum((gw<0.0) & (ttbb>=51))
    #
    scale = np.sign(gw) * scale
    sc_tot = [ sum( scale[:,i] ) for i in range(9)]
    mur_up_tot, mur_down_tot   = sc_tot[7], sc_tot[1]
    muf_up_tot, muf_down_tot   = sc_tot[5], sc_tot[3]
    murf_up_tot, murf_down_tot = sc_tot[8], sc_tot[0]
    sc_ttbb = [ sum( scale[(ttbb>=51)][:,i] ) for i in range(9)]
    mur_up_ttbb, mur_down_ttbb   = sc_ttbb[7], sc_ttbb[1]
    muf_up_ttbb, muf_down_ttbb   = sc_ttbb[5], sc_ttbb[3]
    murf_up_ttbb, murf_down_ttbb = sc_ttbb[8], sc_ttbb[0]
    #
    ps = np.sign(gw) * ps
    ps_tot = [ sum( ps[:,i] ) for i in range(4) ]
    isr_up_tot, isr_down_tot = ps_tot[2], ps_tot[0]
    fsr_up_tot, fsr_down_tot = ps_tot[3], ps_tot[1]
    ps_ttbb = [ sum( ps[(ttbb>=51)][:,i] ) for i in range(4) ]
    isr_up_ttbb, isr_down_ttbb = ps_ttbb[2], ps_ttbb[0]
    fsr_up_ttbb, fsr_down_ttbb = ps_ttbb[3], ps_ttbb[1]
    #
    pdf_up_tot, pdf_down_tot = sum(pdf_up), sum(pdf_down)
    pdf_up_ttbb, pdf_down_ttbb = sum(pdf_up[(ttbb>=51)]), sum(pdf_down[(ttbb>=51)])
    #
    return [tot_count,  # 0
            mur_up_tot, mur_down_tot, muf_up_tot, muf_down_tot, murf_up_tot, murf_down_tot, # 1-6
            isr_up_tot, isr_down_tot, fsr_up_tot, fsr_down_tot, # 7-10
            pdf_up_tot, pdf_down_tot, # 11-12
            ttbb_count,  # 13
            mur_up_ttbb, mur_down_ttbb, muf_up_ttbb, muf_down_ttbb, murf_up_ttbb, murf_down_ttbb, # 14-19
            isr_up_ttbb, isr_down_ttbb, fsr_up_ttbb, fsr_down_ttbb,# 20-23
            pdf_up_ttbb, pdf_down_ttbb] # 24-25
    
def sys_worker(roo):
    if 'kodiak' in roo:
        roo = roo.replace('root://kodiak-se.baylor.edu//','/cms/data/')
    t=uproot.open(roo)['Events']

    gw    = t.array('genWeight'  )
    ttbb  = t.array('genTtbarId')
    ttbb  = ttbb % 100
    #
    tot_count= sum(gw>=0.0) - sum(gw<0.0)
    ttbb_count = sum((gw>=0.0) & (ttbb>=51)) - sum((gw<0.0) & (ttbb>=51))
    #

    return [tot_count,  # 0
            ttbb_count]  # 2
            
 
def process(sample, target_file, xsec, pool):
    with open(target_file) as roo_files:
        roo_list  = [line.strip('\n') for line in roo_files.readlines()]
        #
        #print(roo_list)
        if 'erdOn' in sample or 'hdamp' in sample or 'UE' in sample:
            p_worker = sys_worker
        else:
            p_worker = worker
        print('\nFor sample: {0} with xsec {1}'.format(sample, xsec))
        print(p_worker)
        result_list = np.array(pool.map(p_worker, ( roo for roo in roo_list)))
        #
        if 'erdOn' in sample or 'hdamp' in sample or 'UE' in sample:
            results = [ sum(result_list[:,i]) for i in range(2)]
            print(f'\n For sample: with xsec : Tot Count : tt_bb Yield')
            print('{0}, {1}, {2}, {3}'.format(sample, xsec, *results[:]))
        else:
            results = [ sum(result_list[:,i]) for i in range(26)]
            print(f'\n For sample: with xsec : Tot Count : Total Yield mu r up/down, mu f up/down, mu rf up/down, ISR up/down, FSR up/down : pdf up/down')
            print('{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, {12}, {13}, {14}'.format(sample, xsec, *results[:13]))
            #
            print(f'\n tt_bb Count : tt_bb Yield mu r up/down, mu f up/down, mu rf up/down, ISR up/down, FSR up/down : pdf up/down')
            print('{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, {12}'.format(*results[13:]))
        #


def main(cfg_file,get_s=None):
    with open(sys.argv[1]) as cfg_file:
        pool = multiprocessing.Pool() # have to define this here
        for sample in cfg_file.readlines():
            if '#' is sample[0]: continue
            if sample.replace('\n','') is '' : continue
            if get_s is not None and get_s not in sample: continue
            if 'owheg' not in sample and 'erdOn' not in sample and 'UE' not in sample and 'hdamp' not in sample : continue
            s_name, s_dir, s_file, _, xsec = [ _.strip() for _ in sample.split(',')[:5] ]
            #if not os.path.exists(s_dir+'/'+s_file) : raise NameError(s_dir+'/'+s_file)
            s_dir = 'xrdcp -f root://cmseos.fnal.gov//' + re.search(r'store/[a-zA-Z0-9_/]*', s_dir).group()
            os.system(f'{s_dir}/{s_file} .')
            process(s_name, s_file, xsec, pool)
            os.system(f'rm {s_file}')

if __name__ == '__main__':

    cfg_file = sys.argv[1]
    sample=None
    if len(sys.argv) == 3: sample=sys.argv[2]
    if not os.path.exists(sys.argv[1]) : raise NameError(sys.argv[1])
    #
    main(cfg_file,sample)
