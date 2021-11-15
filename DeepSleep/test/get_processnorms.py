import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import numpy as np
import pandas as pd
import json
import re
from glob import glob
#import functools
import config.ana_cff as cfg
from config.sample_cff import sample_cfg, process_cfg
from modules.AnaDict import AnaDict

'''
main-> __worker --> handle_ttbb/nom/sig --> handle_process
'''

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
            elif 'ttH' in mc or 'ttZ' in mc:
                #print(_sub_worker_dict)
                handle_sig(meta, _sub_worker_dict, sf)
                #print(_sub_worker_dict)
            else:
                _sub_worker_dict[sf] = handle_nom(meta, sample_file)
            #
        if 'ttH' in mc or 'ttZ' in mc:
            handle_sig_process(_sub_worker_dict, mc, _worker_dict)
        else:
            _worker_dict[mc] = handle_process(_sub_worker_dict, mc)
    return _worker_dict
        

def handle_sig_process(meta, mc, in_dict):
    for suf in ['','_0','_1','_2','_3']: 
        sub_meta = {sample+suf: meta[sample+suf] for sample in process_cfg[mc]}
        in_dict[mc+suf.replace('_','')] = handle_process(sub_meta,mc,suf=suf)
        
def handle_process(meta, mc, suf=''): 
    # need to normalize by xsec
    dummy_dict = {} # for ttbb
    tot_xsec = 0.
    if 'ttbb' in mc:
        for sample in process_cfg[mc]:
            tot_xsec += meta[sample]['xs']
            dummy_dict[sample] = {'weight': meta[sample]['weight'], 'xs': meta[sample]['xs']}
            del meta[sample]['xs'], meta[sample]['weight'] # clean up
    else:
        tot_xsec = sum([(sample_cfg['TTZToQQ']['xs']-sample_cfg['TTZToBB']['xs'])*sample_cfg[sample]['kf'] if sample == 'TTZToQQ' else
                        sample_cfg[sample]['xs']*sample_cfg[sample]['kf'] for sample in process_cfg[mc]])
        #tot_xsec = sum([sample_cfg[sample]['xs']*sample_cfg[sample]['kf'] for sample in process_cfg[mc]]) 
    #
    _mc_dict_out = {k:0 for k in meta[process_cfg[mc][0]+suf]} # meta[process_cfg[mc][0] is dummy
    for sample in meta:
        #if sample == 'TTZToBB': continue
        q_sample = re.sub(r'_\d','',sample) if suf != '' else sample# need to get rid of trailing number
        xs = ((sample_cfg['TTZToQQ']['xs']-sample_cfg['TTZToBB']['xs'])*sample_cfg[q_sample]['kf']) if q_sample == 'TTZToQQ' else (sample_cfg[q_sample]['xs']*sample_cfg[q_sample]['kf'])
        norm_f = (xs)/tot_xsec if 'ttbb' not in mc else dummy_dict[sample]['xs']/tot_xsec 
        #
        
        for k in meta[sample]:
            if 'yield' not in k:
                _mc_dict_out[k] += meta[sample][k]*norm_f
            else:
                _mc_dict_out[k] += meta[sample][k]*meta[sample]['weight']
    _mc_dict_out.update(dummy_dict)
    if 'weight' in _mc_dict_out:
        del _mc_dict_out['weight']
    print(mc,'\n',_mc_dict_out)
    return _mc_dict_out
    

def handle_nom(meta, sample_file):
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
    isTTBar = 'sys' not in sample_file and 'TTBar' in sample_file
    if isTTBar:
        _out_dict[f'pdfweight_Up']   = meta['tot_events']/calc_pdfas_unc(meta[f'pdfweights_tot'],'hess','Up') 
        _out_dict[f'pdfweight_Down'] = meta['tot_events']/calc_pdfas_unc(meta[f'pdfweights_tot'],'hess','Down')
        _out_dict[f'alphas_Up']   = meta['tot_events']/calc_pdfas_unc(meta[f'pdfweights_tot'],'alphas','Up') 
        _out_dict[f'alphas_Down'] = meta['tot_events']/calc_pdfas_unc(meta[f'pdfweights_tot'],'alphas','Down')

    return _out_dict
def handle_sig(meta, _in_dict, _sf):
    #
    genpt_bins = [0,200,300,450,np.inf]
    for i,ptbin in enumerate(genpt_bins):
        _f  = (lambda x: x.replace('_tot',f'_{i-1}')) if i > 0 else (lambda x:x)
        _sub_out_dict = {
            'weight'        :meta['weight'],
            'yield'         :meta[_f('events_tot')],   
            'stxs_yield'    :meta[_f('events_rap_tot')],
            'mu_r_Up'       :meta[_f('events_rap_tot')]/meta[_f('mur_up_tot')] ,
            'mu_r_Down'     :meta[_f('events_rap_tot')]/meta[_f('mur_down_tot')] ,
            'mu_f_Up'       :meta[_f('events_rap_tot')]/meta[_f('muf_up_tot')] ,
            'mu_f_Down'     :meta[_f('events_rap_tot')]/meta[_f('muf_down_tot')] ,
            'mu_rf_Up'      :meta[_f('events_rap_tot')]/meta[_f('murf_up_tot')] ,
            'mu_rf_Down'    :meta[_f('events_rap_tot')]/meta[_f('murf_down_tot')] ,
            'ISR_Up'        :meta[_f('events_rap_tot')]/meta[_f('isr_up_tot')] ,
            'ISR_Down'      :meta[_f('events_rap_tot')]/meta[_f('isr_down_tot')] ,
            'FSR_Up'        :meta[_f('events_rap_tot')]/meta[_f('fsr_up_tot')] ,
            'FSR_Down'      :meta[_f('events_rap_tot')]/meta[_f('fsr_down_tot')] ,
            'pdfWeight_Up'  :meta[_f('events_rap_tot')]/meta[_f('pdf_up_tot')] ,
            'pdfWeight_Down':meta[_f('events_rap_tot')]/meta[_f('pdf_down_tot')] ,
            'pdfweight_Up'  :meta[_f('events_rap_tot')]/calc_pdfas_unc(meta[_f('pdfweights_tot')], 'hess', 'Up'),
            'pdfweight_Down':meta[_f('events_rap_tot')]/calc_pdfas_unc(meta[_f('pdfweights_tot')], 'hess', 'Down'),
            'alphas_Up'     :meta[_f('events_rap_tot')]/calc_pdfas_unc(meta[_f('pdfweights_tot')], 'alphas', 'Up'),
            'alphas_Down'   :meta[_f('events_rap_tot')]/calc_pdfas_unc(meta[_f('pdfweights_tot')], 'alphas', 'Down'),
        }
        alt_sf = _sf + f'_{i-1}' if i > 0 else _sf
        _in_dict[alt_sf] = _sub_out_dict
        
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
        'mu_rf_Up'        : get_varnorm('murf_up',  ),
        'mu_rf_Down'      : get_varnorm('murf_down',),
        #
        'pdfWeight_Up'   : get_varnorm('pdf_up',pdf=True),
        'pdfWeight_Down' : get_varnorm('pdf_down',pdf=True),
        #x
        'ISR_Up'         : get_varnorm('isr_up',   ps=True),
        'ISR_Down'       : get_varnorm('isr_down', ps=True),
        'FSR_Up'         : get_varnorm('fsr_up',   ps=True),
        'FSR_Down'       : get_varnorm('fsr_down', ps=True),
    }
    isSys = 'sys' in sample_file

    _out_dict[f'pdfweight_Up']   = (meta['ttbb_events']/calc_pdfas_unc(meta[f'pdfweights_ttbb'],'replica','Up')) \
                                   if not isSys else 0
    _out_dict[f'pdfweight_Down'] = (meta['ttbb_events']/calc_pdfas_unc(meta[f'pdfweights_ttbb'],'replica','Down')) \
                                   if not isSys else 0
    return _out_dict

def calc_pdfas_unc(_pdf, _type='hess', _ud='Up'):
    var_pdf = {
        'hess'   : (lambda _x: np.sqrt(np.sum(np.power(_x[1:-2]-_x[0],2)))),
        'replica': (lambda _x: (np.quantile(_x[1:],.84) - np.quantile(_x[1:],.16)) / 2),
        'alphas' : (lambda _x: (_x[-2] - _x[-1])/2),
    }
    nom_pdf = {
        'hess'   :(lambda _x: _x[0]),
        'replica':(lambda _x: (np.quantile(_x[1:],.84) + np.quantile(_x[1:],.16)) / 2),
        'alphas' :(lambda _x: _x[0]),
    }
    err = var_pdf[_type](_pdf)
    nom = nom_pdf[_type](_pdf)
    unc = np.nan_to_num(err/nom, nan=1)
    _out = _pdf[0]*(1.+unc) if _ud == 'Up' else _pdf[0]*(1.-unc)
    return _out # returns up, down

if __name__ == '__main__':
    main()
