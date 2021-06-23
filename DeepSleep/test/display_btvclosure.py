import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator, FixedLocator, FormatStrFormatter
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)

import json
import numpy as np
import config.ana_cff as cfg
from modules.AnaDict import AnaDict
from config.sample_cff import sample_cfg, process_cfg
from lib.fun_library import save_pdf, getLaLabel

@save_pdf("btvclosure_study.pdf")
def main():
    for year in cfg.Years:
        r_ratio = json.load(open(cfg.dataDir+'/btagw_r_ratio/'+f"btagw_r_ratio_{year}.json", 'r'))
        for p in ['ttZ','ttH','ttbb','TTBar']:
            p_hists = get_hists(p,year)
            apply_rratio(p_hists,np.array(r_ratio[p]['r_ratio']))
            makeplots(p_hists,p,year)

def makeplots(p_hists,p,year):
    k_dict = {
        'nj': 'jet multiplicity',
        'ht':'ht [GeV]',
        'b_score':'deepCSV score (leading jet)',
    }
    bin_dict = {
        'nj'     : np.arange(-0.5,10.5,1),
        'ht'     : np.arange(0,1550,50),
        'b_score': np.arange(0,1.05,0.05),
    }
    for kinem in k_dict:
        fig, ax = initPlt()
        #
        ax[0].step(x=bin_dict[kinem],y=np.append(p_hists[kinem],0),            where='post', color='k',label='no SF')
        ax[0].step(x=bin_dict[kinem],y=np.append(p_hists[kinem+'_sf'],0),      where='post', color='blue',label='SF')
        ax[0].step(x=bin_dict[kinem],y=np.append(p_hists[kinem+'_sf_corr'],0), where='post', color='red',label='SF+corr')
        ax[0].set_xlim(bin_dict[kinem][0],bin_dict[kinem][-1])
        #
        ax[1].axhline(1, color='k', linestyle='--', dashes=(4,8), snap=True)
        ax[1].step(x=bin_dict[kinem],y=np.append(p_hists[kinem+'_sf']/p_hists[kinem],0),      where='post', color='blue',)
        ax[1].step(x=bin_dict[kinem],y=np.append(p_hists[kinem+'_sf_corr']/p_hists[kinem],0), where='post', color='red', )
        #ax[1].set_xlim(bin_dict[kinem][0],bin_dict[kinem][-1])
        #
        endPlt(fig,ax,p,year,kinem=k_dict[kinem])
        #plt.show()

def apply_rratio(hists, r_):
    hists['nj'], hists['nj_sf'] = hists['nj_yield'], hists['btw_yield']
    del hists['nj_yield'], hists['btw_yield'] 
    hists['nj_sf_corr'] = hists['nj_sf']*r_
    #
    hists['htvsn_2dsf_corr'] = hists['htvsn_2dsf']*r_[:,np.newaxis]
    hists['bvsn_2dsf_corr']  = hists['bvsn_2dsf']*r_[:,np.newaxis]
    #
    for t1 in ['','sf','sf_corr']:
        for t2 in ['htvsn','bvsn']: 
            hists[t2.replace('vsn','').replace('b','b_score')+t1.replace('sf','_sf')] = np.sum(hists[f'{t2}_2d{t1}'],axis=0)
            del hists[f'{t2}_2d{t1}']

def get_hists(p,year):
    out_dict = {
        'nj_yield':None,
        'htvsn_2d':None,
        'bvsn_2d':None,
        'btw_yield':None,
        'htvsn_2dsf':None,
        'bvsn_2dsf':None,
    }
    for mc in process_cfg[p]:
        mD = AnaDict.read_pickle(f'{cfg.postSkim_dir}/{year}/{p}/{mc}.pkl')['metaData']
        out_dict = update_dict(mD,out_dict)
    return out_dict

def update_dict(md_,out_):
    out_ = {k: md_[k]*md_['weight'] if out_[k] is None else md_[k]*md_['weight']+out_[k] for k in out_}
    return out_


def initPlt():
    fig, ax = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
    fig.subplots_adjust(top=0.88,bottom=0.11,left=0.11,right=0.90,hspace=0.0, wspace=0.0)
    return fig, ax

def endPlt(fig,ax,title,year,kinem=''):
    fig.text(0.110,0.89, r"$\bf{CMS}$ $Preliminary $", fontsize = 12)
    fig.text(0.695,0.89, f'{round(cfg.Lumi[year],1)}'+r' fb$^{-1}$ (13 TeV)',  fontsize = 12)
    fig.text(0.215,0.835, rf"Process: {title.replace('TTBar','ttlf')}", fontsize=12)
    #
    ax[0].legend(framealpha = 0)
    ax[0].set_ylabel('Events / bin')
    ax[1].set_ylabel('Ratio')
    plt.xlabel(kinem)
    ax[0].grid(color='k', linestyle=':', alpha=0.25)
    #
    ax[0].tick_params(which='both', direction='in', top=True, right=True)
    ax[0].xaxis.set_minor_locator(AutoMinorLocator())
    ax[0].yaxis.set_minor_locator(AutoMinorLocator())
    #
    ax[1].xaxis.set_minor_locator(AutoMinorLocator())
    ax[1].yaxis.set_major_formatter(FormatStrFormatter('%g'))
    ax[1].yaxis.set_major_locator(FixedLocator([.50,1,1.50]))
    ax[1].set_ylim(0.2,1.8)
    ax[1].yaxis.set_minor_locator(AutoMinorLocator())
    ax[1].tick_params(which='both', direction='in', top=True, right=True)
    ax[1].grid(color='k', linestyle=':', alpha=0.25)
    


if __name__ == '__main__':
    main()
