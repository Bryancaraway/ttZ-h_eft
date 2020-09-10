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
    t=uproot.open(roo)['Events']

    gw    = t.array('genWeight'  )
    ttbb  = t.array('genTtbarId')
    ttbb  = ttbb % 100
    scale = t.array('LHEScaleWeight').pad(9).fillna(1)
    ps    = t.array('PSWeight'      ).pad(4).fillna(1)
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
    


    return [tot_count,  # 0
            mur_up_tot, mur_down_tot, muf_up_tot, muf_down_tot, murf_up_tot, murf_down_tot, # 1-6
            isr_up_tot, isr_down_tot, fsr_up_tot, fsr_down_tot, # 7-10
            ttbb_count,  # 11
            mur_up_ttbb, mur_down_ttbb, muf_up_ttbb, muf_down_ttbb, murf_up_ttbb, murf_down_ttbb, # 12-17
            isr_up_ttbb, isr_down_ttbb, fsr_up_ttbb, fsr_down_ttbb] # 18-21
 
def process(sample, target_file, xsec, pool):
    with open(target_file) as roo_files:
        roo_list  = [line.strip('\n') for line in roo_files.readlines()]
        #
        #print(roo_list)
        print('\nFor sample: {0} with xsec {1}'.format(sample, xsec))
        result_list = np.array(pool.map(worker, ( roo for roo in roo_list)))
        #p_count, n_count = sum(result_list[:,0]), sum(result_list[:,1])
        results = [ sum(result_list[:,i]) for i in range(22)]
        #
        print(f'Tot Count: {results[0]}')
        print('Total Yield Mu r up/down {0}/{1}, mu f up/down {2}/{3}, mu rf up/down {4}/{5}'.format(*results[1:7]))
        print('Total Yield ISR up/down {0}/{1}, FSR up/down {2}/{3}'.format(*results[7:11]))
        #
        print(f'Tot ttbb Count: {results[11]}')
        print('Tot ttbb Yield Mu r up/down {0}/{1}, mu f up/down {2}/{3}, mu rf up/down {4}/{5}'.format(*results[12:18]))
        print('Tot ttbb Yield ISR up/down {0}/{1}, FSR up/down {2}/{3}'.format(*results[18:]))


def main(cfg_file,get_s=None):
    with open(sys.argv[1]) as cfg_file:
        pool = multiprocessing.Pool() # have to define this here
        for sample in cfg_file.readlines():
            if '#' is sample[0]: continue
            if sample.replace('\n','') is '' : continue
            if get_s is not None and get_s not in sample: continue
            if 'owheg' not in sample and 'erdOn' not in sample and 'UE' not in sample and 'hdamp' : continue
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
