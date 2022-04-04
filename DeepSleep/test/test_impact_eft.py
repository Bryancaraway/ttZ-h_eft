import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')
import pandas as pd
import numpy as np
import re
from lib.fun_library import save_pdf, import_mpl_settings
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import rc
from matplotlib.ticker import AutoMinorLocator
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)                                                            
import config.ana_cff as cfg

'''
Script that demonstrates how wc impacts
our expected signal in terms of pt distribution
'''

nn = cfg.nn

kbins = {
    'Zh_pt': [200,300,450,600],
    'Zh_M' : cfg.sdm_bins,
    #nn     : [0.,   0.21 , 0.58 , 0.70 , 0.82 , 0.90 , 1.00],
    nn : { # for 2018 only
        '200':[0.0,0.17,0.51,0.63,0.78,0.88,1.00],
        '300':[0.0,0.21,0.58,0.70,0.82,0.90,1.00],
        '450':[0.0,0.31,0.69,0.78,0.87,0.92,1.00],},
    'reco_pt_bin': [1,2,3],
    'reco_m_bin': [0,1,2,3],
    'gen_pt_bin': [0,1,2,3],
    'genZHpt': [200,300,450,550,650],
}
klabel = {
    'Zh_pt': (lambda ax: (ax.set_xlabel('reco. Z/H $p_{T}$ (GeV)'), 
                          ax.set_ylabel(r'(d$\sigma^{EFT}$$/$d$p_{T}$) $/$ (d$\sigma^{SM}$$/$d$p_{T}$)'))),
    'Zh_M' : (lambda ax: (ax.set_xlabel('reco. Z/H $m_{sd}$ (GeV)'), 
                          ax.set_ylabel(r'(d$\sigma^{EFT}$$/$d$m_{sd}$) $/$ (d$\sigma^{SM}$$/$d$m_{sd}$)'))),
    nn     : (lambda ax: (ax.set_xlabel('NN'), 
                          ax.set_ylabel(r'(d$\sigma^{EFT}$$/$d$NN$) $/$ (d$\sigma^{SM}$$/$d$NN$)'))),
    'gen_pt_bin' : (lambda ax: (ax.set_xlabel('gen Z/H $p_{T}$ bin'), 
                                ax.set_ylabel(r'(d$\sigma^{EFT}$$/$d$p_{T}$) $/$ (d$\sigma^{SM}$$/$d$p_{T}$)'))),
    'genZHpt' : (lambda ax: (ax.set_xlabel('gen Z/H $p_{T}$ (GeV)'), 
                             ax.set_ylabel(r'(d$\sigma^{EFT}$$/$d$p_{T}$) $/$ (d$\sigma^{SM}$$/$d$p_{T}$)'))),
}
         
class TestEFTImpact():
    fit_f   = 'EFT_Parameterization_test.npy'
    wc_ranges = {
        #old'ctW'  :[-2.44, 2.40]  ,
        #old'ctZ'  :[-2.33, 2.31]  ,
        #old'ctp'  :[-11.64, 41.58],
        #old'cpQM' :[-12.33, 16.04],
        #old'ctG'  :[-0.50, 0.48]  ,
        #old'cbW'  :[-7.40, 7.40]  ,
        #old'cpQ3' :[-6.94, 7.20]  ,
        #old'cptb' :[-16.67, 16.67],
        #old'cpt'  :[-20.35, 16.53]

        'ctW'  :[-1.03, 0.94]  ,
        'ctZ'  :[-0.99, 1.02]  ,
        'ctp'  :[0.13, 30.18],
        'cpQM' :[-4.77, 5.63],
        #'ctG'  :[-0.50, 0.48]  ,
        'cbW'  :[-4.55, 4.65]  ,
        'cpQ3' :[-3.86, 2.88]  ,
        'cptb' :[-9.44, 10.15],
        'cpt'  :[-8.50, 5.39]
        }
    samples=['ttbb','ttbbjet','ttjets','ttZ','ttH']
    df = pd.read_pickle(fit_f)

    @save_pdf('ttbb_wc_single_effect.pdf')
    def run_singles(self, years=['2018']):
        norm_change = {}
        for y in years:
            df = self.df[y]['ttH']
            for wc, r in self.wc_ranges.items():
                norm_change[wc] = [self.calc_norm(df,wc,v) for v in r]
    
        wcs   = np.array([w for w in norm_change])
        norms = np.array([norm_change[w] for w in norm_change])
        ind   = np.array([ [i,i] for i in range(len(norm_change))])
        for i, w in enumerate(norm_change):
            plt.scatter([i,i],norm_change[w],c='k')
            plt.plot([i,i],norm_change[w],c='k')
        #plt.scatter(ind.flatten(), norms.flatten())
        plt.axhline(1,c='r',ls='--')
        plt.fill_between([-10,10],1.21,0.82,color='r',alpha=0.5,label=f'Exp $1\sigma$ run2 limit')
        plt.fill_between([-10,10],1.42,0.64,color='k',alpha=0.25,label=f'Exp $2\sigma$ run2 limit')
        plt.xticks([i for i in range(len(norm_change))], wcs)
        plt.ylabel(r'$\sigma$(ttbb)$^{EFT}$/$\sigma$(ttbb)$^{SM}$')
        plt.xlim(-1,10)
        plt.title('Rate change w/resp. to WC for tt+bb')
        plt.legend()
        plt.show()
        exit()
    


    def run_singles_jon(self, kinem, samples, cut, title, year='2018', sepSig=False, sepSigreco=False, pt_int='200'):
        master_s = samples[0]
        genpt_samples   = self.handle_sepSig(samples, year, 'gen_pt_bin')
        recopt_samples  = self.handle_sepSig(samples, year, 'reco_pt_bin')
        recoptm_samples = self.handle_sepSig_pt_mass(samples, year)
        for wc, r in self.wc_ranges.items():
            #ax2=ax.twinx()
            bins = np.array(kbins[kinem][pt_int] if kinem == nn else kbins[kinem])
            for i in range(1,len(r)):
                fig,ax = self.initPlot()
                genptsc_sumw = []
                genptsc_sumw2 = []
                recoptsc_sumw = []
                recoptsc_sumw2 = []
                recoptmsc_sumw = []
                recoptmsc_sumw2 = []
                def get_totals(samples_,sumw_,sumw2_):
                    for s in samples_:
                        df = self.df[year][s]
                        df = cut(df)
                        try: 
                            nsm, nsm2 = self.get_sumw_sumw2(df[kinem].clip(bins[0],bins[-1]), df['SM']*df['weight'] , bins)
                        except ValueError:
                            nsm, nsm2 = self.get_sumw_sumw2(df[kinem], df['SM']*df['weight'] , bins)
                        eftw = self.get_eftw(df, wc, r[i])
                        #print(s,nsm,sum(eftw*df['weight'])/sum(df['SM']*df['weight']) )
                        sumw_.append(nsm   * sum(eftw*df['weight'])/sum(df['SM']*df['weight']) )
                        sumw2_.append(nsm2 * ((sum(eftw*df['weight']))/(sum(df['SM']*df['weight'])))**2)
                get_totals(genpt_samples, genptsc_sumw, genptsc_sumw2)
                get_totals(recopt_samples, recoptsc_sumw, recoptsc_sumw2)
                get_totals(recoptm_samples, recoptmsc_sumw, recoptmsc_sumw2)
                #
                df = self.df[year][master_s]
                df = cut(df)
                eftw = self.get_eftw(df, wc, r[i])
                try:
                    nsm, nsm2 = self.get_sumw_sumw2(df[kinem].clip(bins[0],bins[-1]), df['SM']*df['weight'] , bins)                
                    neft, neft2 = self.get_sumw_sumw2(df[kinem].clip(bins[0],bins[-1]), eftw*df['weight'] , bins)
                except ValueError:
                    nsm, nsm2 = self.get_sumw_sumw2(df[kinem], df['SM']*df['weight'] , bins)
                    neft, neft2 = self.get_sumw_sumw2(df[kinem], eftw*df['weight'] , bins)

                ax.step(
                    x = bins,
                    y = np.append(sum(genptsc_sumw)/nsm,0),
                    where='post', label=f'{master_s} gen pt scaled'
                )
                #ax.step(
                #    x = bins,
                #    y = np.append(sum(recoptsc_sumw)/nsm,0),
                #    where='post', label=f'{master_s} reco pt scaled'
                #)
                ax.step(
                    x = bins,
                    y = np.append(sum(recoptmsc_sumw)/nsm,0),
                    where='post', label=f'{master_s} reco pt&mass scaled'
                )
                ax.step(
                    x = bins,y = np.append(neft/nsm,0),
                    where='post', label=f'{master_s} raw EFT'
                )
                ax.step(
                    x = bins,y = np.append(nsm/nsm,0),
                    where='post', label=f'{master_s} SM'
                )
                ax.errorbar(**self.errbar_kwargs(nsm, nsm2, sum(genptsc_sumw), sum(genptsc_sumw2), bins))
                #ax.errorbar(**self.errbar_kwargs(nsm, nsm2, sum(recoptsc_sumw), sum(recoptsc_sumw2), bins))
                ax.errorbar(**self.errbar_kwargs(nsm, nsm2, sum(recoptmsc_sumw), sum(recoptmsc_sumw2), bins+.005))
                ax.errorbar(**self.errbar_kwargs(nsm, nsm2, neft, neft2, bins-.005))
                ax.errorbar(**self.errbar_kwargs(nsm, nsm2, nsm, nsm2, bins))

                #ax2.step(x=bins, y = np.append(nsm,0), where='post',linestyle='--')
                #ax2.errorbar(x=(bins[1:]+bins[:-1])/2, y = nsm, yerr=np.sqrt(nsm2), fmt='.')
                
                #for i,xy in enumerate(zip((b[1:]+b[:-1])/2,neft/nsm)):
                #ax.text(x=xy[0],y=xy[1], s=f'{slabel} {tot_sm[i]:.1f}:{tot_eft[0][i]:.1f}:{tot_eft[1][i]:.1f}', fontsize=8)
                #
                fig.suptitle(title.format(wc=wc, r=r))
                ax.legend()
                self.endPlot(ax,title,kinem,bins)

    #@save_pdf('eft_bkgsig_comparison_nn09.pdf')
    def run_singles_test(self, kinem, samples, cut, title, year='2018', sepSig=False, sepSigreco=False, pt_int='200'):
        if sepSig or sepSigreco:
            samples = self.handle_sepSig(samples, year, 'reco_pt_bin' if sepSigreco else 'gen_pt_bin')
        for wc, r in self.wc_ranges.items():
            fig,ax = self.initPlot()
            #ax2=ax.twinx()
            #dott2b = True
            bins = np.array(kbins[kinem][pt_int] if kinem == nn else kbins[kinem])
            for j,s in enumerate(samples):                
                #sam_text = f'{s}'
                df = self.df[year][s]
                slabel = s
                df = cut(df)
                if 'ttbb' in s:
                    df = df[df['process'] == 'tt_B']
                    slabel = s
                    #if dott2b:
                    #    df = df[df['process'] == 'tt_2b']
                    #    slabel = 'tt_2b' 
                    #    dott2b=False # only do this for the first ttbb
                    #else: 
                    #    df = df[df['process'] == 'tt_bb']
                elif s == 'ttjets':
                    df = df[df['process'] == 'TTBar']
                    slabel = 'tt_lf'
            
                try: 
                    nsm, nsm2 = self.get_sumw_sumw2(df[kinem].clip(bins[0],bins[-1]), df['SM']*df['weight'] , bins)
                except ValueError:
                    nsm, nsm2 = self.get_sumw_sumw2(df[kinem], df['SM']*df['weight'] , bins)
                #
                for i in range(0,len(r)):
                    eftw = self.get_eftw(df, wc, r[i])
                    try:
                        neft, neft2 = self.get_sumw_sumw2(df[kinem].clip(bins[0],bins[-1]), eftw*df['weight'] , bins)
                    except ValueError:
                        neft, neft2 = self.get_sumw_sumw2(df[kinem], eftw*df['weight'] , bins)
                    #tot_eft[i] = (neft/nsm)*tot_sm
                    ax.step(
                        x = bins,
                        y = np.append(neft/nsm,0),
                        where='post', label=f'{slabel}:{r[i]}'
                    )
                    ax.errorbar(**self.errbar_kwargs(nsm, nsm2, neft, neft2, bins+j*5))

                #ax2.step(x=bins, y = np.append(nsm,0), where='post',linestyle='--')
                #ax2.errorbar(x=(bins[1:]+bins[:-1])/2, y = nsm, yerr=np.sqrt(nsm2), fmt='.')
                
                #for i,xy in enumerate(zip((b[1:]+b[:-1])/2,neft/nsm)):
                #ax.text(x=xy[0],y=xy[1], s=f'{slabel} {tot_sm[i]:.1f}:{tot_eft[0][i]:.1f}:{tot_eft[1][i]:.1f}', fontsize=8)
            #
            fig.suptitle(title.format(wc=wc, r=r))
            self.endPlot(ax,title,kinem,bins)

    def errbar_kwargs(self,nsm, nsm2, neft, neft2, bins):
        y    = neft/nsm
        #yerr = y*np.sqrt( nsm2/np.power(nsm,2) +  neft2/np.power(neft,2) )
        yerr = np.sqrt(neft2)/nsm
        kwargs = {'x':(bins[1:]+bins[:-1])/2, 'y':y, 'yerr':yerr, 'fmt':'.'}
        return kwargs

    def get_sumw_sumw2(self,df_,w_,bins):
        sumw, _ = np.histogram(df_, bins=bins, weights=w_)
        sumw2, _ = np.histogram(df_, bins=bins, weights=np.power(w_,2))
        return sumw, sumw2
            
    def initPlot(self):
        fig, ax = plt.subplots(1,1, sharey=True)
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.08,
            right=0.90,
            hspace=0.0,
            wspace=0.0
        )
        return fig,ax

    def endPlot(self,ax,title,kinem, bins):
        ax.legend(ncol=2, fontsize='x-small')
        ax.set_xlim(bins[0],bins[-1])
        ax.set_ylim(0.4)
        klabel[kinem](ax)
        #ax2.set_ylabel('SM yield')
        #ax2.set_yscale('log')
        #ax2.set_ylim(.1)
        #fig.suptitle(f'EFT impact, NN > 0.9, {wc}:{r}')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        plt.tight_layout()

    def get_eftw(self, df, wc, v):
        p = df[f'{wc}_{wc}']*v*v
        q = df[f'{wc}']*v
        r = df['SM']
        _eft_w = p + q + r
        return _eft_w

    def calc_norm(self,df,wc,v):
        p = sum(df[f'{wc}_{wc}']*v*v)
        q = sum(df[f'{wc}']*v)
        r = sum(df['SM'])
        return p/r + q/r + r/r
    
    def getMCStat_info(self,n, b, h, w):
        #x = (self.bins[1:]+self.bins[:-1])/2
        x = (b[1:]+b[:-1])/2
        xerr = (b[1:]-b[:-1])/2
        yerr = np.histogram(h, bins=b, weights=np.power(w,2))[0]
        yerr = np.sqrt(yerr)
        y = n
        return x,xerr,y,yerr

    def handle_sepSig(self,samples, year, pt_bin):
        sigs = re.findall(r'(ttZ|ttH)', ' '.join(samples))
        r_o_g = pt_bin.split('_')[0]
        bin_range = range(kbins[pt_bin][0],kbins[pt_bin][-1]+1)
        out_  = []
        for sig in sigs:
            out_ += [f'{sig}{r_o_g}pt{i}' for i in bin_range]
            sub_df = (lambda df_, i : df_[year][sig][df_[year][sig][pt_bin] == i])
            self.df[year].update({f'{sig}{r_o_g}pt{i}': sub_df(self.df, i) for i in bin_range})
        return out_
    def handle_sepSig_pt_mass(self,samples,year):
        sigs = re.findall(r'(ttZ|ttH)', ' '.join(samples))
        pt_bin = 'reco_pt_bin'
        m_bin = 'reco_m_bin'
        pt_bin_range = range(kbins[pt_bin][0],kbins[pt_bin][-1]+1)
        m_bin_range = range(kbins[m_bin][0],kbins[m_bin][-1]+1)
        out_  = []
        for sig in sigs:
            out_ += [f'{sig}reco{i}{j}' for i in pt_bin_range for j in m_bin_range]
            sub_df = (lambda df_, i,j : df_[year][sig][(df_[year][sig][pt_bin] == i) & (df_[year][sig][m_bin] == j)])
            self.df[year].update({f'{sig}reco{i}{j}': sub_df(self.df, i, j) for i in pt_bin_range for j in m_bin_range})
        #
        return out_

@save_pdf("eft_impacts_template_ttbb2017.pdf")
def plot_template_impacts(eft_impacts):
    pt_bins = [200,300,450,np.inf]
    mass_bins = [50,75,105,145,200]
    #for sig in ['ttZ','ttH','ttbb']:
    for sig in ['ttbb']:
        for i_bin in range(2,len(pt_bins)):
            #for i_bin in range(1,len(pt_bins)):
            for j_bin in range(1,len(mass_bins)):
                cut_fun = (lambda _x : _x[(_x['Zh_pt']>pt_bins[i_bin-1]) & (_x['Zh_pt']<pt_bins[i_bin]) & (_x['Zh_M']>mass_bins[j_bin-1]) & (_x['Zh_M']<mass_bins[j_bin]) & (_x[nn]>0.0)])
                eft_impacts.run_singles_test(nn, [sig], cut_fun, f'Z/H pT {i_bin}, Z/H mass {j_bin},'+' {wc}', pt_int=str(pt_bins[i_bin-1]), year='2017')
        

    
if __name__=='__main__':
    import_mpl_settings(1)
    eft_impacts = TestEFTImpact()
    #eft_impacts.run_singles()
    cut_nntight = (lambda x : x[(x[nn] > .8)])
    cut_nnloose = (lambda x : x[(x[nn] > 0)])
    cut_nnloose_rapid =  (lambda x : x[(x[nn] > 0) & (x['genZHstxs']==1)])
    plot_template_impacts(eft_impacts)
    exit()
    #eft_impacts.run_singles_test(['ttH','ttZ','ttbb','ttbb'],cut)
    #@save_pdf('eft_bkgsig_comparison_nn09.pdf')
    #save_pdf('eft_ttZ_masscomparison_nn00.pdf')(eft_impacts.run_singles_test)('Zh_M', ['ttZ'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=True)
    #save_pdf('eft_ttH_masscomparison_nn00.pdf')(eft_impacts.run_singles_test)('Zh_M', ['ttH'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=True)
    #save_pdf('eft_ttZ_nncomparison_nn00.pdf')(eft_impacts.run_singles_test)(nn, ['ttZ'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=True)
    #save_pdf('eft_ttH_nncomparison_nn00.pdf')(eft_impacts.run_singles_test)(nn, ['ttH'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=True)
    
    # ---- Jons study
    #save_pdf('eft_ttZ_ptcomparison4Jon_nn00.pdf')(eft_impacts.run_singles_jon)('Zh_pt', ['ttZ'],cut_nnloose_rapid, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=True)
    #save_pdf('eft_ttH_ptcomparison4Jon_nn00.pdf')(eft_impacts.run_singles_jon)('Zh_pt', ['ttH'],cut_nnloose_rapid, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=True)
    #save_pdf('eft_ttZ_masscomparison4Jon_nn00.pdf')(eft_impacts.run_singles_jon)('Zh_M', ['ttZ'],cut_nnloose_rapid, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=True)
    #save_pdf('eft_ttH_masscomparison4Jon_nn00.pdf')(eft_impacts.run_singles_jon)('Zh_M', ['ttH'],cut_nnloose_rapid, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=True)
    #save_pdf('eft_ttZ_nncomparison4Jon_nn00.pdf')(eft_impacts.run_singles_jon)(nn, ['ttZ'],cut_nnloose_rapid, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=True)
    #save_pdf('eft_ttH_nncomparison4Jon_nn00.pdf')(eft_impacts.run_singles_jon)(nn, ['ttH'],cut_nnloose_rapid, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=True)
    # ---- Jons study


    #save_pdf('eft_ttZ_pt_nocomparison_nn00.pdf')(eft_impacts.run_singles_test)('Zh_pt', ['ttZ'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=False)
    #save_pdf('eft_ttH_pt_nocomparison_nn00.pdf')(eft_impacts.run_singles_test)('Zh_pt', ['ttH'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=False)
    #save_pdf('eft_ttZ_mass_nocomparison_nn00.pdf')(eft_impacts.run_singles_test)('Zh_M', ['ttZ'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=False)
    #save_pdf('eft_ttH_mass_nocomparison_nn00.pdf')(eft_impacts.run_singles_test)('Zh_M', ['ttH'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=False)
    #save_pdf('eft_ttZ_nn_nocomparison_nn00.pdf')(eft_impacts.run_singles_test)(nn, ['ttZ'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=False)
    #save_pdf('eft_ttH_nn_nocomparison_nn00.pdf')(eft_impacts.run_singles_test)(nn, ['ttH'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=False)
    #
    #
    save_pdf('eft_ttbb_masscomparison_nn00.pdf')(eft_impacts.run_singles_test)('Zh_M', ['ttbb'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}')    
    save_pdf('eft_ttbb_ptcomparison_nn00.pdf')(eft_impacts.run_singles_test)('Zh_pt',  ['ttbb'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}')    
    save_pdf('eft_ttbb_nncomparison_nn00.pdf')(eft_impacts.run_singles_test)(nn,       ['ttbb'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}')    
    #
    #cut_nnloose = (lambda x : x[((x[nn] > 0) & (x['SM'] < 1))])
    #save_pdf('eft_tt_bkg_masscomparison_nn00.pdf')(eft_impacts.run_singles_test)('Zh_M', ['ttbb', 'ttjets'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}')    
    #save_pdf('eft_tt_bkg_ptcomparison_nn00.pdf')(eft_impacts.run_singles_test)('Zh_pt',  ['ttbb', 'ttjets'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}')    
    #save_pdf('eft_tt_bkg_nncomparison_nn00.pdf')(eft_impacts.run_singles_test)(nn,       ['ttbb', 'ttjets'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}')    

    #### OLD
    #save_pdf('eft_ttZ_genptcomparison_nn00.pdf')(eft_impacts.run_singles_test)('gen_pt_bin', ['ttZ'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=True)
    #save_pdf('eft_ttH_genptcomparison_nn00.pdf')(eft_impacts.run_singles_test)('gen_pt_bin', ['ttH'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}', sepSig=True)
    #save_pdf('eft_bkgsig_masscomparison_nn08.pdf')(eft_impacts.run_singles_test)(['ttH','ttZ','ttbb','ttbb'],cut_nntight, 'EFT impact, NN > 0.8, {wc}:{r}')
    #save_pdf('eft_bkgsig_masscomparison_nn00.pdf')(eft_impacts.run_singles_test)(['ttH','ttZ','ttbb','ttbb','ttjets'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}')
    #save_pdf('eft_bkg_masscomparison_nn00.pdf')(eft_impacts.run_singles_test)('Zh_pt', ['ttbb','ttjets'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}')
    #save_pdf('eft_bkg_masscomparison_nn08.pdf')(eft_impacts.run_singles_test)(['ttbb','ttbb'],cut_nntight, 'EFT impact, NN > 0.8, {wc}:{r}')
    #
    #save_pdf('eft_ttbb_jet_masscomparison_nn00.pdf')(eft_impacts.run_singles_test)('Zh_M', ['ttbb','ttbbjet'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}')    
    #save_pdf('eft_ttbb_jet_ptcomparison_nn00.pdf')(eft_impacts.run_singles_test)('Zh_pt',  ['ttbb','ttbbjet'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}')    
    #save_pdf('eft_ttbb_jet_nncomparison_nn00.pdf')(eft_impacts.run_singles_test)(nn,       ['ttbb','ttbbjet'],cut_nnloose, 'EFT impact, NN > 0.0, {wc}:{r}')    
