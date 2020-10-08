#from ROOT import TFile                                                                                         
import uproot
import sys
import os
import concurrent.futures
import multiprocessing
import re
import numpy as np
executor = concurrent.futures.ThreadPoolExecutor()

def worker(roo):
    if 'kodiak' in roo:
        roo = roo.replace('root://kodiak-se.baylor.edu//','/cms/data/')
    t=uproot.open(roo)['Events']
    gw = t.array('genWeight')#, executor=executor, blocking=True)
    scale = t.array('LHEScaleWeight').pad(9).fillna(1) * np.sign(gw)
    pdf_up   = t.array('pdfWeight_Up') * np.sign(gw)
    pdf_down = t.array('pdfWeight_Down') * np.sign(gw)
    #
    sc_tot = [ sum( scale[:,i] ) for i in range(9)]
    mur_up_tot, mur_down_tot   = sc_tot[7], sc_tot[1]
    muf_up_tot, muf_down_tot   = sc_tot[5], sc_tot[3]
    murf_up_tot, murf_down_tot = sc_tot[8], sc_tot[0]
    #
    pdf_up_tot, pdf_down_tot = sum(pdf_up), sum(pdf_down)
    #
    p_count= sum(gw>=0.0)
    n_count= sum(gw<0.0)
    tot_count = p_count - n_count
    return [p_count, n_count, # 0-1
            tot_count,        # 2
            mur_up_tot, mur_down_tot, muf_up_tot, muf_down_tot, murf_up_tot, murf_down_tot, # 3-8
            pdf_up_tot, pdf_down_tot # 9-10
    ]
 
def process(sample, target_file, pool):
    with open(target_file) as roo_files:
        roo_list  = [line.strip('\n') for line in roo_files.readlines()]
        #
        print('\nFor sample: {}'.format(sample))
        result_list = np.array(pool.map(worker, ( roo for roo in roo_list)))
        #results = sum(result_list[:,0]), sum(result_list[:,1])
        results  = [ sum(result_list[:,i]) for i in range(11) ]
        #
        #print(results)
        #
        print('p_events, n_events: {}, {}'.format(results[0], results[1]))
        print('Sample : Tot events : u r up/down, mu f up/down, mu rf up/down : pdf up/down')
        print('{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}'.format(sample, *results[2:]))


def main(cfg_file,get_s=None):
    with open(sys.argv[1]) as cfg_file:
        pool = multiprocessing.Pool() # have to define this here
        for sample in cfg_file.readlines():
            if '#' is sample[0]: continue
            if sample.replace('\n','') is '' : continue
            if get_s is not None and get_s not in sample: continue
            s_name, s_dir, s_file, = [ _.strip() for _ in sample.split(',')[:3] ]
            #if not os.path.exists(s_dir+'/'+s_file) : raise NameError(s_dir+'/'+s_file)
            print('uscms' in s_dir)
            if 'uscms' in s_dir:
                s_dir = 'xrdcp -f root://cmseos.fnal.gov//' + re.search(r'store/[a-zA-Z0-9_/]*', s_dir).group()            
                print(f'{s_dir}/{s_file} .')
                os.system(f'{s_dir}/{s_file} .')
            else: 
                s_file = f'{s_dir}/{s_file}'
            #print(f'{s_dir}/{s_file} .')
            #os.system(f'{s_dir}/{s_file} .')
            process(s_name, s_file, pool)
            os.system(f'rm {s_file}')

if __name__ == '__main__':

    cfg_file = sys.argv[1]
    sample=None
    if len(sys.argv) == 3: sample=sys.argv[2]
    if not os.path.exists(sys.argv[1]) : raise NameError(sys.argv[1])
    #
    main(cfg_file,sample)
    

