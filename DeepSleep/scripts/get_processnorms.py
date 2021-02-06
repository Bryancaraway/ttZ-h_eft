import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import uproot
import numpy as np
import pandas as pd
import json
import re
from glob import glob
#import functools
import config.ana_cff as cfg
from config.sample_cff import sample_cfg, process_cfg
from modules.AnaDict import AnaDict

def main():
    json_dict = {}
    for y in cfg.Years:
        json_dict[y] = __worker(y)
    out_json_file = f'{cfg.dataDir}/process_norms/process_norms_ttbbw_run2.json'
    with open(out_json_file,'w') as jsf:
        json.dump( json_dict,jsf, indent=4)

def __worker(y):
    _worker_dict = {}
    sample_pdir = cfg.postSkim_dir+f'{y}/'
    for mc in cfg.All_MC+['ttbb_sys']:
        if 'Data' in mc: continue # do this later
        _sub_worker_dict = {}
        for sf in process_cfg[mc]:
            sample_file = sample_pdir+f'{mc}/{sf}.pkl'
            meta = AnaDict.read_pickle(sample_file)['metaData']
            if 'ttbb' in mc:
                _sub_worker_dict[sf] = handle_ttbb(meta, sample_file,y)
            else:
                _sub_worker_dict[sf] = handle_nom(meta)
            #
        _worker_dict[mc] = handle_process(_sub_worker_dict, mc)
    return _worker_dict
        
        
def handle_process(meta, mc): 
    # need to normalize by xsec
    dummy_dict = {} # for ttbb
    tot_xsec = 0.
    if 'ttbb' in mc:
        for sample in process_cfg[mc]:
            tot_xsec += meta[sample]['xs']
            dummy_dict[sample] = {'weight': meta[sample]['weight'], 'xs': meta[sample]['xs']}
            del meta[sample]['xs'], meta[sample]['weight'] # clean up
    else:
        tot_xsec = sum([sample_cfg[sample]['xs']*sample_cfg[sample]['kf'] for sample in process_cfg[mc] if sample != 'TTZToBB']) 
    #
    _mc_dict_out = {k:0 for k in meta[process_cfg[mc][0]]} # meta[process_cfg[mc][0] is dummy
    for sample in meta:
        if sample == 'TTZToBB': continue
        norm_f = (sample_cfg[sample]['xs']*sample_cfg[sample]['kf'])/tot_xsec if 'ttbb' not in mc else dummy_dict[sample]['xs']/tot_xsec 

        for k in meta[sample]:
            _mc_dict_out[k] += meta[sample][k]*norm_f
    _mc_dict_out.update(dummy_dict)
    print(mc,'\n',_mc_dict_out)
    return _mc_dict_out
    

def handle_nom(meta):
    #
    _out_dict = {
        'mu_r_Up'       :meta['tot_events']/meta['mur_up_tot'] ,
        'mu_r_Down'     :meta['tot_events']/meta['mur_down_tot'] ,
        'mu_f_Up'       :meta['tot_events']/meta['muf_up_tot'] ,
        'mu_f_Down'     :meta['tot_events']/meta['muf_down_tot'] ,
        'mu_rf_Up'      :meta['tot_events']/meta['murf_up_tot'] ,
        'mu_rf_Down'    :meta['tot_events']/meta['murf_down_tot'] ,
        'pdfWeight_Up'  :meta['tot_events']/meta['pdf_up_tot'] ,
        'pdfWeight_Down':meta['tot_events']/meta['pdf_down_tot'] ,
    }
    return _out_dict
    
def handle_ttbb(meta,sample_file,y):
    ttbar_name = sample_file.split('/')[-1].replace('TTbb_','TTTo').replace('.pkl','')
    ttbar_meta = AnaDict.read_pickle(sample_file.replace('TTbb_','TTTo').replace('ttbb','TTBar'))['metaData']
    lumi = cfg.Lumi[y]*1000
    xsec = sample_cfg[ttbar_name]['xs']
    get_weight =  (lambda : (ttbar_meta['ttbb_events']/meta['ttbb_events'])*xsec/ttbar_meta['tot_events'] * lumi)
    def get_varnorm(pre, ps=False, pdf=False):
        var_frac_5FS = ttbar_meta[f'{pre}_ttbb']/ttbar_meta[f'{pre}_tot']
        nom_frac_5FS = ttbar_meta[f'ttbb_events']/ttbar_meta[f'tot_events']
        var_norm_4FS = meta['ttbb_events']/meta[f'{pre}_ttbb']
        #renorm_f     = ttbar_meta[f'tot_events']/ttbar_meta[f'{pre}_tot'] if renorm else 1
        #renorm_f     = ttbar_meta[f'ttbb_events']/ttbar_meta[f'{pre}_ttbb'] if renorm else 1
        renorm_f     = 1 
        if ps:
            return  (ttbar_meta[f'{pre}_ttbb']/ttbar_meta[f'ttbb_events'])/(meta[f'{pre}_ttbb']/meta[f'ttbb_events'])
        if pdf:
            return meta['ttbb_events']/meta[f'{pre}_ttbb']
        return var_norm_4FS*(var_frac_5FS/nom_frac_5FS)*renorm_f
        
        
    _out_dict = {
        'weight'         : get_weight(), # works for hdamp as well
        'xs'             : get_weight()*meta['ttbb_events']/lumi,
        #
        'mu_r_Up'        : get_varnorm('mur_up',  ),
        'mu_r_Down'      : get_varnorm('mur_down',),
        'mu_f_Up'        : get_varnorm('muf_up',  ),
        'mu_f_Down'      : get_varnorm('muf_down',),
        #
        'pdfWeight_Up'   : get_varnorm('pdf_up',pdf=True),
        'pdfWeight_Down' : get_varnorm('pdf_down',pdf=True),
        #
        'ISR_Up'         : get_varnorm('isr_up',   ps=True),
        'ISR_Down'       : get_varnorm('isr_down', ps=True),
        'FSR_Up'         : get_varnorm('fsr_up',   ps=True),
        'FSR_Down'       : get_varnorm('fsr_down', ps=True),
    }
    return _out_dict

if __name__ == '__main__':
    main()
