import uproot
import os
import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)',
        shell=True).decode().strip('\n')+'DeepSleep/')
#
import re
from multiprocessing import Pool
#from pathos.multiprocessing import ProcessingPool as Pool
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
import matplotlib.backends.backend_pdf as matpdf
from matplotlib import rc
rc("savefig",dpi=250)
rc("figure", figsize=(6, 6*(6./8.)), dpi=200)                                                            
rc("figure", max_open_warning=600)
#
import functools
import numpy as np
import pandas as pd
import config.ana_cff as cfg
from lib.fun_library import t2Run, save_pdf, getZhbbWeight, getLaLabel
from makeDatacard import MakeDataCard

nn = 'withbbvl_NN'

target = ['Zh_M','Zh_pt']
doNNcuts = True
tbins_map = {
    'Zh_M' :[50,75,90,105,120,140,200], # new format
    'Zh_pt':[200,300,450,np.inf],
    'Zh_score':[0,.25,.50,1.0],
    'nBottoms_drLeptonCleaned':[1.5,2.5,3.5,4.5,10],
    'n_ak4jets':[3.5,4.5,5.5,6.5,7.5,14],
    'n_b_inZh':[-0.5,0.5,1.5,2.5]}


def get_sumw_sumw2(df, weights, year):
    #ret_target = (lambda df: (df[nn].to_numpy(), df['Zh_M'].to_numpy()))
    #sumw = np.array([
    #    hist2d[year][i_bin]( *ret_x_y(df[df['pt_bin']==i_bin]), 
    #                         weights=weights[df['pt_bin']==i_bin])[0] 
    #    for i_bin in range(1,len(pt_bins))
    #])
    #sumw2 = np.array([
    #    hist2d[year][i_bin]( *ret_x_y(df[df['pt_bin']==i_bin]), 
    #                         weights=np.power(weights[df['pt_bin']==i_bin],2))[0]
    #    for i_bin in range(1,len(pt_bins))
    #])
    sumw  = np.histogram( df[target].to_numpy(),
                          bins = tbins_map[target],
                          weights=weights)[0]
    sumw2 = np.histogram( df[target].to_numpy(),
                          bins = tbins_map[target],
                          weights=np.power(weights,2))[0]

    return sumw, sumw2
#


class PlotsFromDC(MakeDataCard):
    '''
    borrow methods used to get samples with variables
    '''
    pt_bins = [0,200,300,450]
    def __init__(self, 
                 sig = cfg.Sig_MC+cfg.sig_sys_samples, 
                 bkg = cfg.Bkg_MC+cfg.bkg_sys_samples,  # does not include QCD
                 years = cfg.Years, 
                 isblind=False):
        # start getting signal and bkg 
        super().__init__(sig,bkg,years,isblind,get_sumw_sumw2)
        self.sig = sig
        self.bkg = bkg
        self.data = cfg.Data_samples
        self.years =  years
        #self.years = ['2018']
        self.isblind = isblind
        self.dc_bins = 1#len(tbins_map[target])#len(pt_bins[1:])
        self.bkg_v  = super().weights + super().weight_sys + ['process',nn,'sample']+target
        self.sig_v  = self.bkg_v + ['genZHpt']
        self.data_v = ['process',nn]+target
        #
        #self.getdata()
    def getData(self):
        super().getdatav2()
        return self.data_dict

    #def updatedict(self, p, v, y=''):
    #    if doNNcuts:
    #        v = v +[nn]
    #    sub_f_dir = 'data_files' if 'Data' in p else 'mc_files'
    #    if not os.path.exists(f'{self.file_dir}{y}/{sub_f_dir}/{p}_val.pkl'): return 
    #    df = pd.read_pickle(f'{self.file_dir}{y}/{sub_f_dir}/{p}_val.pkl').filter(items=v)
    #    if 'TTbb' in p : df = pd.concat([df,pd.read_pickle(f'{self.file_dir}{y}/{sub_f_dir}/{p}_val.pkl').filter(regex=r'\w*_weight')],axis='columns') # to get special ttbb normalizations
    #    if 'Data' not in p: # clip mu_rf, isr/fsr, pdf at 3sigma percentile, 99.7% (.15%,99.85%)
    #        func = np.nanpercentile # because this is a long name
    #        [df[v_str].clip(func(df[v_str].values,.15), func(df[v_str].values,99.85), inplace=True ) 
    #         for v_str in [ f'{s}_{ud}' for s in ['mu_r','mu_f','mu_rf','ISR','FSR','pdfWeight'] for ud in ['Up','Down']] ]
    #    #df['Zh_pt'].clip(pt_bins[0]+1,pt_bins[-1]+1, inplace=True)
    #    #df['pt_bin'] = pd.cut(df['Zh_pt'], bins=pt_bins+[500],
    #    #                      labels=[i_bin for i_bin in range(len(pt_bins))])
    #    if doNNcuts:
    #        group = df[df[nn] >= 0.0].groupby(by='process')
    #    else:
    #        group = df.groupby(by='process')
    #    # extract sys type (if any)
    #    sys = '' 
    #    if p in cfg.all_sys_samples: 
    #        sys =  '_'+p.split('_')[-1] # format [process_name]_[systype]
    #        if 'hdamp' in sys and 'TTbb' in p: sys = sys.replace('hdamp','hdamp_ttbb') 
    #        #if 'JES' in sys or 'JER' in sys:   sys = sys.replace('Up',f'_{y}Up').replace('Down',f'_{y}Down')
    #        if 'JES' in sys or 'JER' in sys or 'UE' in sys or 'hdamp' in sys:   sys = sys.replace('Up',f'_{y}Up').replace('Down',f'_{y}Down')
    #    data_dict = {f"{n.replace('Data','data_obs')}_{y}{sys}": g for n,g in group} # iterate over name and content
    #    return data_dict


def initDF():
    data = PlotsFromDC().getData()
    #processes = ['ttZ','ttH','TTBar','tt_bb','tt_2b','ttX','VJets','other']
    processes = ['ttZ','ttH','TTBar','tt_bb','tt_2b','single_t','ttX','VJets']
    #rare_p = {'ttX','VJets','other'}    
    k_list = ['Zh_M','Zh_pt',nn,'process']
    df = pd.DataFrame()
    for p in processes:
        for y in ['2016','2017','2018']:
            data[f'{p}_{y}']['event_weight'] = getZhbbWeight(data[f'{p}_{y}'],y)
            data[f'{p}_{y}']['event_weight2'] = data[f'{p}_{y}']['event_weight']**2
            #if p in rare_p:
            #    data[f'{p}_{y}']['process'] = 'rare'
            df = pd.concat([df,data[f'{p}_{y}'].filter(items=k_list+['event_weight','event_weight2'])], 
                           axis='rows', ignore_index=True)
    #
    # organize by process
    cuts = (lambda df_: df_[(df_[nn] > 0.80) & (df_['Zh_pt'] > 450)])
    df = cuts(df)
    ordered_p = ['ttZ','ttH','tt_bb','tt_2b','TTBar','single_t','ttX','VJets']
    getv = (lambda k: [df[k][df['process'] == p].values for p in ordered_p])
    #
    h = getv('Zh_M')
    w = getv('event_weight')
    w2 = getv('event_weight2')
    i     = [np.sum(w_i) for w_i in w]
    i_err     = [np.sqrt(np.sum(w2_i)) for w2_i in w2]
    #i_err = np.sqrt(np.sum(w2,axis=0))
    l,c   = np.array([np.array(getLaLabel(p)) for p in ordered_p]).T
    c = [_.replace('gold','k') for _ in c]
    l   = np.array([f'{x} ({y:3.1f}+/-{z:3.1f})' for x,y,z in zip(l,i,i_err)])
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.11,
            right=0.88,
            hspace=0.2,
            wspace=0.2
        )
    n, bins, _ = ax.hist(
        h,
        bins=tbins_map['Zh_M'],
        histtype='step',
        weights = w,
        color   = c,
        label    = l
    )
    addMCStat(ax,h,n,w,i,c,bins)
    #
    from matplotlib.ticker import AutoMinorLocator
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both', direction='in', top=True, right=True)
    #self.fig.text(0.105,0.89, r"$\bf{CMS}$ $Simulation$", fontsize = self.fontsize)
    fig.text(0.105,0.89, r"$\bf{CMS}$ $Preliminary $", fontsize = 12)
    fig.text(0.635,0.89, f'137'+r' fb$^{-1}$ (13 TeV)',  fontsize = 12)
    ax.set_xlabel('Z/H $m_{sd}$ ($NN > 0.80$; Z/H $p_{t}>450$) (GeV)', fontsize = 12)
    #self.ax.set_ylabel(f"{'%' if self.doNorm else 'Events'} / {(self.bin_w[0].round(2) if len(np.unique(self.bin_w.round(4))) == 1 else 'bin')}")#fontsize = self.fontsize)
    ax.set_ylabel('Events / bin')
    #print(self.ax.get_ylim()[1], self.ax.get_ylim()[1] * 1.10 )
        
    #plt.xlim(self.bin_range)
    ax.set_yscale('log')
    ax.set_xlim(bins[0],bins[-1])
    ax.set_ylim(ymin=ax.get_ylim()[0],ymax=ax.get_ylim()[1]*7.5)
    ax.legend(framealpha = 0, ncol=2, fontsize='xx-small')
    #plt.grid(True)
    #plt.setp(patches_, linewidth=0)
    #handles, labels = self.ax.get_legend_handles_labels()
    #hatch_patch = Patch(hatch=10*'X', label='Stat Unc.',  fc='w')
    #handles = handles + [hatch_patch]
    #labels  = labels + ['Stat Unc.']
    #self.ax.legend(handles,labels, framealpha = 0, ncol=2, fontsize='xx-small')
    #self.ax.set_ylim(ymin=self.ax.get_ylim()[0],ymax=self.ax.get_ylim()[1]*(10 if self.doLog else 1.50))
    #if self.doSave: plt.savefig(f'{self.saveDir}{self.xlabel}_.pdf', dpi = 450)
    #if self.doShow: plt.show()
    plt.savefig('pdf/zhm_80nn_450pt_run2.pdf')

def addMCStat(ax,h,n,w,i,c,bins):
    bin_c = (bins[1:]+bins[:-1])/2
    yerr = np.array([np.histogram(h[s], bins=bins, weights=np.power(w[s],2))[0] for s in range(len(h))])
    for j in range(len(i)):
        ax.errorbar(x=bin_c, y=n[j], 
                    yerr= (np.sqrt(yerr[j])),
                    fmt='.', ms=3, color=c[j])
    #
#




def main():
    df = initDF()


if __name__ == '__main__':
    main()
