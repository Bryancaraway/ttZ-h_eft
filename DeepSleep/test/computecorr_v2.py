import uproot
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.gridspec as gridspec

plt.rc("font", size=10, family="sans-serif", **{"sans-serif" : [u'TeX Gyre Heros', u'Helvetica', u'Arial']})
#plt.rc("xaxis", labellocation='right')
#plt.rc("yaxis", labellocation='top')
plt.rc("legend", fontsize=10, scatterpoints=1, numpoints=1, borderpad=0.15, labelspacing=0.2,
                 handlelength=0.7, handletextpad=0.25, handleheight=0.7, columnspacing=0.6,
                 fancybox=False, edgecolor='none', borderaxespad=0.15)
plt.rc("savefig", dpi=200)
#plt.rc("figure", figsize=(3.375, 3.375*(6./8.)), dpi=200)
plt.rc("text", usetex=True)
plt.rc("text.latex", preamble='\n'.join([r"\usepackage{amsmath}",
                                         r"\usepackage{helvet}",
                                         r"\usepackage{sansmath}",
                                         r"\sansmath"]))
#plt.rc("hatch", linewidth=0.3)
import pandas as pd
import numpy as np
import re

c_dict = {
    'stxs': [p+str(i) for p in ['r_ttZ','r_ttH'] for i in range(4)],
    'inc' : ['r_ttZ','r_ttH'],
}

l_dict = {
    'r_ttZ' : r'${\mu}_{\mathrm{ttZ}}$',
    'r_ttZ0': '[0,200]',
    'r_ttZ1': '[200,300]',
    'r_ttZ2': '[300,450]',
    'r_ttZ3': r'$[450,\infty]$',
    'r_ttH' : r'${\mu}_{\mathrm{ttH}}$',
    'r_ttH0': '[0,200]',
    'r_ttH1': '[200,300]',
    'r_ttH2': '[300,450]',
    'r_ttH3': r'$[450,\infty]$',
}



def main():
    #for _file,_list in c_dict.items():
    #dodiffCorr('robustHesse_stxs.root', c_dict['stxs'])
    fdir = 'fitdiag_roots/corr/'
    dodiffCorr(fdir+'robustHesse_stxs.root', c_dict['stxs'])
    doincCorr( fdir+'robustHesse_inc.root',  c_dict['inc'])
    doincCorr( fdir+'robustHesse_inc_blind.root',  c_dict['inc'])

class TH2D:

    def __init__(self, th2d_from_root):
        self.th2d = th2d_from_root
    
    def getbylabels(self, str_1, str_2):
        i = self.th2d.xlabels.index(str_1)
        j = self.th2d.ylabels.index(str_2)
        return self.th2d.values[i,j]
    

def dodiffCorr(i_file, c_list):
    roo =  uproot.open(i_file)
    h_corr = TH2D(roo["h_correlation"])
    #print(h_corr.getbylabels('r_ttZ0','r_ttZ2'))
    print("Post-fit correlation between Signal POIs")
    df = pd.DataFrame()
    for c_i in c_list:
        for c_j in c_list:
            df.loc[c_i,c_j] = h_corr.getbylabels(c_i,c_j)  
    df=df.reindex(index=df.index[::-1])
    print(df)

    fig, (ax) = plt.subplots(2,2, figsize=(6.75,6.75))#, sharex=True)#, sharey=True)#,  sharex=True)
    fig.subplots_adjust(
        top=0.85,
        bottom=0.14,
        left=0.14,
        right=0.85,
        wspace=0.0,
        hspace=0.0
    )
    ax[0,0].set_aspect("equal")
    ax[0,1].set_aspect("equal")
    ax[1,0].set_aspect("equal")
    ax[1,1].set_aspect("equal")
    CMSlabel(fig,ax)
    _  = ax[0,1].matshow(df.iloc[:4,4:], cmap='bwr', vmin=-1, vmax=1) # upper-left 
    _  = ax[1,0].matshow(df.iloc[4:,:4], cmap='bwr', vmin=-1, vmax=1) # lower-left
    _  = ax[0,0].matshow(df.iloc[:4,:4], cmap='bwr', vmin=-1, vmax=1) # upper-right
    im = ax[1,1].matshow(df.iloc[4:,4:], cmap='bwr', vmin=-1, vmax=1) # lower-left
    # annotate, y, x, axes
    annotate_ax(df, c_list[4:][::-1], c_list[:4], ax[0,0])
    annotate_ax(df, c_list[4:][::-1], c_list[4:], ax[0,1])
    annotate_ax(df, c_list[:4][::-1], c_list[:4], ax[1,0])
    annotate_ax(df, c_list[:4][::-1], c_list[4:], ax[1,1])
    #im = ax.matshow(df, cmap='bwr', vmin=-1, vmax=1)
    ## We want to show all ticks...
    ax[0,0].tick_params('x', direction='in', which='both',top=True, bottom=False, labeltop=False)
    ax[0,0].tick_params('y', direction='in', which='both',left=True, right=False)
    ax[0,1].tick_params('x', direction='in', which='both',top=True, bottom=False, labeltop=False)
    ax[0,1].tick_params('y', direction='in', which='both', left=False, right=True, labelleft=False)
    #ax[0,1].yaxis.set_visible(False)
    ax[1,0].tick_params('x', direction='in', which='both',top=False, bottom=True, labeltop=False)
    ax[1,0].tick_params('y', direction='in', which='both', left=True, right=False)
    ax[1,1].tick_params('x', direction='in', which='both',top=False, bottom=True, labeltop=False)
    ax[1,1].tick_params('y', direction='in', which='both', left=False, right=True, labelleft=False)
    #ax[1,1].yaxis.set_visible(False)    
    #
    ax[1,0].set_xticks(range(len(c_list[:4])))
    ax[1,1].set_xticks(range(len(c_list[4:])))
    ax[0,0].set_yticks(range(len(c_list[:4][::-1])))
    ax[1,0].set_yticks(range(len(c_list[4:][::-1])))
    ## ... and label them with the respective list entries
    l_list = [l_dict[c] for c in c_list]
    ax[1,0].set_xticklabels(l_list[:4], fontsize=10)
    ax[1,1].set_xticklabels(l_list[4:], fontsize=10)
    ax[1,0].xaxis.set_ticks_position('bottom')
    ax[1,1].xaxis.set_ticks_position('bottom')
    ax[0,0].set_yticklabels(l_list[:4][::-1], fontsize=10)
    ax[1,0].set_yticklabels(l_list[4:][::-1], fontsize=10)    
    #
    ax[1,0].set_xlabel(r'${\sigma}_{\mathrm{ttZ}}$')
    ax[1,1].set_xlabel(r'${\sigma}_{\mathrm{ttH}}$')
    ax[0,0].set_ylabel(r'${\sigma}_{\mathrm{ttH}}$')
    ax[1,0].set_ylabel(r'${\sigma}_{\mathrm{ttZ}}$')

    ##
    #
    cbar_ax = fig.add_axes([0.86,0.14,0.04,0.71])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cticks = np.arange(-1.0,1.25,.25)
    cbar.set_ticks(cticks)
    cbar_ax.set_ylabel('Correlation')
    ##
    cbar.set_ticklabels(["{0:.2f}".format(i) for i in cticks])
    #plt.show()
    plt.savefig("../pdf/stxs_new_corr.pdf")
    plt.close('all')
    #plt.minorticks_off()


def annotate_ax(df,list1,list2, ax):
    for i,c_i in enumerate(list1):
        for j,c_r in enumerate(list2):
            text = ax.text(j, i, round(df.loc[c_i,c_r],2),
                           ha="center", va="center", color="k")

def doincCorr(i_file, c_list):
    roo =  uproot.open(i_file)
    h_corr = TH2D(roo["h_correlation"])
    #print(h_corr.getbylabels('r_ttZ0','r_ttZ2'))
    print("Post-fit correlation between Signal POIs")
    df = pd.DataFrame()
    for c_i in c_list:
        for c_j in c_list:
            df.loc[c_i,c_j] = h_corr.getbylabels(c_i,c_j)  
    print(df)
    # print out correlation matrix between signal and tt+LF, tt+bb
    print(['r_ttH','r_ttZ','tt_qsc','CMS_ttbbnorm'])
    for i in ['r_ttH','r_ttZ','tt_qsc','CMS_ttbbnorm']:
        h_line = i+': '
        for j in ['r_ttH','r_ttZ','tt_qsc','CMS_ttbbnorm']:
            h_line += str(h_corr.getbylabels(i,j))+'\t'
        print(h_line)
    #

    fig, ax = plt.subplots(figsize=(6.75,6.75))
    fig.subplots_adjust(
        top=0.85,
        bottom=0.14,
        left=0.14,
        right=0.85,
        wspace=0.0,
        hspace=0.0
    )
    CMSlabel(fig,ax)
    im = ax.matshow(df, cmap='bwr', vmin=-1, vmax=1)
    # We want to show all ticks...
    ax.set_xticks(range(len(c_list)))
    ax.set_yticks(range(len(c_list)))
    # ... and label them with the respective list entries
    l_list = [l_dict[c] for c in c_list]
    ax.set_xticklabels(l_list, fontsize=10)
    ax.xaxis.set_ticks_position('bottom')
    ax.set_yticklabels(l_list, fontsize=10)
    #
    for i,c_i in enumerate(c_list):
        for j,c_j in enumerate(c_list):
            text = ax.text(j, i, round(df.loc[c_i,c_j],2),
                           ha="center", va="center", color="k")

    cbar_ax = fig.add_axes([0.86,0.14,0.04,0.71])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cticks = np.arange(-1.0,1.25,.25)
    cbar.set_ticks(cticks)
    cbar.set_ticklabels(["{0:.2f}".format(i) for i in cticks])
    cbar_ax.set_ylabel('Correlation')
    #plt.tight_layout()
    #plt.show()
    plt.savefig(f"../pdf/inc_new_corr{'_blinded' if 'blind' in i_file else ''}.pdf")
    plt.close('all')
    #plt.minorticks_off()


def upperlefttext(s):
    trans = gca().transAxes + matplotlib.transforms.ScaledTranslation(3/72, -3/72, gcf().dpi_scale_trans)
    return text(0, 1, s, transform=trans, ha='left', va='top')

def CMSlabel(fig=None, ax=None):
    if fig is None:
        fig = plt.gcf()
    if ax is None:
        ax = plt.gca()

    if type(ax) is np.ndarray:
        ax0, ax1 = ax[0,0], ax[0,1]
    #
    else:
        ax0, ax1 = ax, ax
    trans = ax0.transAxes + matplotlib.transforms.ScaledTranslation(0/72, 3/72, fig.dpi_scale_trans)
    ax0.text(0, 1, r"\textbf{CMS} {\footnotesize \textit{Preliminary}}",
            transform=trans, ha='left', va='baseline')
    trans = ax1.transAxes + matplotlib.transforms.ScaledTranslation(0/72, 3/72, fig.dpi_scale_trans)
    ax1.text(1, 1, r"$138\,\mathrm{fb}^{\text{-1}}$ (13 TeV)",
            transform=trans, ha='right', va='baseline')

if __name__ == '__main__':
    main()

