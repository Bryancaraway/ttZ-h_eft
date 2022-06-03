import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
import matplotlib.pyplot as plt
from matplotlib import rc
rc("figure", max_open_warning=600)
#rc("figure", figsize=(8, 6*(6./8.)), dpi=200)

import json
import numpy as np
import pandas as pd
import config.ana_cff as cfg
from lib.fun_library import t2Run, save_pdf, getZhbbBaseCuts, getZhbbWeight, getLaLabel, import_mpl_settings, upperlefttext, CMSlabel

mfp = cfg.master_file_path
ana_cuts = getZhbbBaseCuts

#pt_bins = [0,200,300,450,np.inf]
pt_bins = [200,300,450,np.inf]


def make_heatmap(df,x_bins,y_bins,p_type, opt_title=None, c_label='Exptected Yield'):
    fig, ax = plt.subplots()
    fig.subplots_adjust(
        top=0.88,
        bottom=0.11,
        left=0.11,
        right=0.88,
        #hspace=0.0 
        wspace=0.2
    )
    print(df)
    #if opt_title is not None:
    #    fig.suptitle(opt_title.format(p_type))
    #CMSlabel(fig, ax, opt='')
    CMSlabel(fig,ax,altloc=False,opt='Simulation', lumi='nl')
    #CMSlabel(fig,ax,altloc=False,opt='Simulation Preliminary', lumi='nl')
    c = ax.pcolor(x_bins,y_bins,df.to_numpy(), vmin=0., vmax=df.to_numpy().max()*1.20)
    #c = ax.pcolor(eff_dict['pt_eta_sf'])
    for i in range(1,len(x_bins)):
        for j in range(1,len(y_bins)):
            text = ax.text(x_bins[i] - (x_bins[i]-x_bins[i-1])/2, 
                           y_bins[j] - (y_bins[j]-y_bins[j-1])/2, 
                           f'{round(df.iloc[j-1,i-1],2):.2f}' if df.to_numpy().max() < 1 else f'{round(df.iloc[j-1,i-1],1):.1f}',
                           ha="center", va="center", color="w")
    #
    #ax.set_xscale("Log")
    ax.tick_params('both', direction='in', which='both', color='k', top=True, right=True)
    ax.set_xticks(x_bins)
    ax.set_xticklabels(['200','300','450',r'$\infty$'] if len(x_bins) == 4 else [str(i) for i in x_bins])
    ax.tick_params(axis='x', pad=5)
    ax.set_yticks(y_bins)
    ax.set_yticklabels(['200','300','450',r'$\infty$'] if len(y_bins) == 4 else [str(i) for i in y_bins])
    plt.minorticks_off()
    #
    #ax.set_xlabel(r"Reconstructed ${p}_{\mathrm{T}}^{\mathrm{Z/H\; cand.}}$ [GeV]", fontsize=10)
    #ax.set_ylabel(rf"Simulated ${{p}}_{{\mathrm{{T}}}}^{{\mathrm{{{p_type}}}}}$ [GeV]", usetex=True, fontsize=10)
    ax.set_xlabel(r"Reconstructed $\mathsf{p}_{\text{T}}^{\text{Z/H\;cand.}}$ \raisebox{0.25ex}{[}$\text{GeV}$\raisebox{0.25ex}{]}", fontsize=10, usetex=True)
    ax.set_ylabel(rf"Simulated $\mathsf{{p}}_{{\text{{T}}}}^{{\text{{{p_type}}}}}$ "+r"\raisebox{0.25ex}{[}$\text{GeV}$\raisebox{0.25ex}{]}", usetex=True, fontsize=10)
    #fig.suptitle(f'{lep}_{year}')
    cbar = fig.colorbar(c, pad=.01, format='%.2f' if df.to_numpy().max() < 1 else '%.1f')
    cbar.ax.tick_params('y', direction='in', which='both', right=True)
    cbar.ax.set_ylabel(c_label, fontsize=10)
    plt.tight_layout()
    #plt.show()


@save_pdf("stxs_sigsens_yields.pdf")
def make_stxs_sigsens_yields(sel):
    x_bins = [200,300,450,600]
    y_bins = pt_bins[:-1] + [600]
    #y_bins = [200,300,450,600] #[0,200,300,450,600]
    for p in ['ttZ','ttH']:
        r_df = sel[p]
        make_heatmap(r_df.clip(0,np.inf), x_bins, y_bins, p.replace('tt',''), opt_title='NN > 0.8, {} mass bin')        

#@save_pdf("stxs_sigsens_response_paper.pdf")
@save_pdf("stxs_sigsens_response_paper_final.pdf")
def make_stxs_sigsens_response(inc,sel):
    x_bins = [200,300,450,600]
    y_bins = pt_bins[:-1] + [600]
    #y_bins = [200,300,450,600] #[0,200,300,450,600]
    for p in ['ttZ','ttH']:
        inc_y = np.array([inc[p][p+str(i)] for i in range(5-len(y_bins),4)])
        r_df = sel[p].divide(inc_y,axis='rows') * 100
        make_heatmap(r_df.clip(0,np.inf), x_bins, y_bins, p.replace('tt',''), c_label='Folding matrix ${M}_{ij}$ [\%]')        


@save_pdf("stxs_yields.pdf")
def make_stxs_yields(inc,sel):
    x_bins = [200,300,450,600]
    y_bins = pt_bins[:-1] + [600]
    #y_bins = [200,300,450,600]#[0,200,300,450,600]
    for p in ['ttZ','ttH']:
        r_df = sel[p]
        make_heatmap(r_df.clip(0,np.inf), x_bins, y_bins, p.replace('tt',''))        
    
@save_pdf("stxs_response_paper.pdf")
def make_stxs_response(inc,sel):
    x_bins = [200,300,450,600]
    #y_bins = [200,300,450,600]#[0,200,300,450,600]
    y_bins = pt_bins[:-1] + [600]
    for p in ['ttZ','ttH']:
        #inc_y = np.array([inc[p][p+str(i)] for i in range(1,4)])#range(4)])
        inc_y = np.array([inc[p][p+str(i)] for i in range(5-len(y_bins),4)])
        r_df = sel[p].divide(inc_y,axis='rows') * 100
        make_heatmap(r_df.clip(0,np.inf), x_bins, y_bins, p.replace('tt',''), c_label='Folding matrix ${M}_{ij}$ [\%]')        


@save_pdf("inc_response.pdf")
def make_inc_response(inc,sel):
    x_bins = [200,300,450,600]
    y_bins = [0,200,600]
    for p in ['ttZ','ttH']:
        inc_y = np.array( [ inc[p][p+str(0)] ] + [ sum([inc[p][p+str(i)] for i in range(1,4)]) ]  )
        temp_df = pd.concat([sel[p].iloc[0,:],sel[p].iloc[1:,:].sum(axis='rows')], axis='columns').T
        r_df = temp_df.divide(inc_y,axis='rows')

        make_heatmap(r_df.clip(0,np.inf), x_bins, y_bins, p.replace('tt','')) 

def get_sel_yields(sel_cuts=None):
    sel_cuts = {'ttZ': (lambda _:_), 'ttH':(lambda _:_)} if sel_cuts is None else sel_cuts
    _out_dict = {}
    for p in ['ttZ','ttH']:
        o_df = pd.DataFrame(np.zeros(shape=(len(pt_bins)-1,3)),
                            columns=[(pt_bins[i-1],pt_bins[i]) for i in range(len(pt_bins)-3,len(pt_bins))],#range(2,len(pt_bins))], # reco pt
                            index=[(pt_bins[i-1],pt_bins[i]) for i in range(1,len(pt_bins))],   # gen  pt
        )
        for y in cfg.Years:
            p_file = f'{mfp}/{y}/mc_files/{p}_val.pkl'
            p_df = pd.read_pickle(p_file)
            p_df['tot_weight'] = getZhbbWeight(p_df,y)#.filter(items=['Zh_pt','genZHpt','genZHstxs'])
            p_df = sel_cuts[p](p_df)
            p_df = p_df[((p_df['genZHstxs'] == True) & (p_df['process']==p) & (ana_cuts(p_df)==True))].filter(items=['Zh_pt','genZHpt','tot_weight'])
            #p_df = p_df[(ana_cuts(p_df)==True)].filter(items=['Zh_pt','genZHpt','tot_weight'])
            #for i in range(1,len(pt_bins)): # gen pt 
            #    for j in range(1,len(pt_bins)):#2,len(pt_bins)): # reco pt 
            #        print(i,j)
            #        pt_cut = ((p_df['Zh_pt'] >= pt_bins[j-1]) & (p_df['Zh_pt'] < pt_bins[j]) &
            #                  (p_df['genZHpt'] >= pt_bins[i-1]) & (p_df['genZHpt'] < pt_bins[i]))
            #        #
            #        o_df.iloc[i-1, j-1] = o_df.iloc[i-1, j-1] + sum(p_df.loc[(pt_cut==True), 'tot_weight'])
            for i_index, i in enumerate(o_df.index): # gen pt 
                for j_column, j in enumerate(o_df.columns):# reco pt 
                    pt_cut = ((p_df['Zh_pt'] >= j[0]) & (p_df['Zh_pt'] < j[1]) &
                              (p_df['genZHpt'] >= i[0]) & (p_df['genZHpt'] < i[1]))
                    #
                    o_df.iloc[i_index, j_column] = o_df.iloc[i_index, j_column] + sum(p_df.loc[(pt_cut==True), 'tot_weight'])
            #
        #
        _out_dict[p] = o_df
    return _out_dict

def get_tot_yields():
    _out_dict = {}
    p_df = json.load(open(cfg.dataDir+'/process_norms/process_norms_ttbbw_run2.json','r'))
    for p in ['ttZ','ttH']:
        o_df = {**{p:0},**{p+str(i):0 for i in range(4)}}
        for y in cfg.Years:
            o_df[p] = o_df[p]+p_df[y][p]['stxs_yield']
            for i in range(4):
                o_df[p+str(i)] = o_df[p+str(i)]+p_df[y][p+str(i)]['stxs_yield']
            #
        #
        _out_dict[p] = o_df
    return _out_dict
            
def main():
    nn = cfg.nn
    import_mpl_settings(i=1)
    sel_y      = get_sel_yields()
    sigsens_y  = get_sel_yields(sel_cuts = {'ttZ':(lambda df : df[((df[nn]>0.8) & (df['Zh_M']>80)  & (df['Zh_M']<115))]),
                                            'ttH':(lambda df : df[((df[nn]>0.8) & (df['Zh_M']>115) & (df['Zh_M']<155))])})
    tot_y      = get_tot_yields()
    #
    #make_stxs_yields(tot_y,sel_y)
    #make_stxs_sigsens_yields(sigsens_y)
    #
    make_stxs_sigsens_response(tot_y,sigsens_y)
    #make_stxs_response(tot_y, sel_y)
    #make_inc_response(tot_y, sel_y)
    

if __name__ == '__main__':
    main()

