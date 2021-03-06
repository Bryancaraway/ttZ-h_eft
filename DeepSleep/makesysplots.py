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
from lib.fun_library import t2Run, save_pdf
'''
Script to grab all histograms from datacard root files
and plot the effect of each systematic per process
'''
# globals
dc_dir = 'Higgs-Combine-Tool'
file_pre = 'datacard' # format file_pre_[year].txt/.root
years    = ['2016','2017','2018']
tag = '' if len(sys.argv) < 2 else sys.argv[1]+'_'

sys_mult = {'2016':
            {'ak8JERUp':
             {'Zhpt1': 2,
              'Zhpt2': 2,
              'Zhpt3': 0,}
         },
            '2017':
            {'ak8JERUp':
             {'Zhpt1': 3.5,
              'Zhpt2': 1.25,
              'Zhpt3': 0.5,}
         },
            '2018':
            {'ak8JERUp':
             {'Zhpt1': 2,
              'Zhpt2': .75,
              'Zhpt3': 0.5,}
         }
}

@t2Run
def main():
    for y in years:
        print(y)
        roofile = f'{dc_dir}/{file_pre}_{tag}{y}.root'
        with uproot.open(roofile) as roo:
            histos = {hist.decode().replace(';1',''):rHist(roo[hist]) for hist in roo} # store hist objects in dict
            #print(histos.keys()) # hist name format is ALWAYS bin_process_sys or bin_process (nominal)
            processes = set()
            uncert    = set()
            bins      = set()
            for k in histos.keys():
                process = k.replace(k.split('_')[0]+'_','')
                bins.add(k.split('_')[0])
                sys = re.search(r'(_(hdamp|pref|mu|mu_r|mu_f|trigeff|btglf|btghf|lepsf|pdf|isr|fsr))?_[a-zA-Z0-9]*(Up|Down|On)',process)
                
                if sys is not None: 
                    sys= sys.group()
                    process = process.replace(sys , '')
                    uncert.add(sys.lstrip('_').replace('Up','').replace('Down',''))
                processes.add(process)
            #
            PlotSystematics(processes,uncert,bins,histos,y)
        #

class rHist(object):
    '''
    trim down uproot TH1D objects 
    to interpretable object class
    '''
    def __init__(self, hist):
        self.edges = hist.edges
        self.values = hist.values
        self.variances = hist.variances

class PlotSystematics:
    #
    def __init__(self,process_list,uncert,bins,histos,year):
        self.process_list = sorted(list(process_list)) # this doesnt need to be a set anymore
        self.uncert = sorted(list(uncert))
        self.bins   = sorted(list(bins))
        self.histos = histos
        self.year = year
        # methods 
        # add pdf-er here
        #save_pdf(f'datacard_{tag}{self.year}_sys.pdf')(self.run_makeplots)(self.makeplots)
        #save_pdf(f'datacard_{tag}{self.year}_sys_bypt.pdf')(self.run_makeplots)(self.makeplots_bypt)
        save_pdf(f'datacard_{tag}{self.year}_sys_bychannel.pdf')(self.run_makeplots)( functools.partial(self.makeplots_bychannel, syss=['UE','hdamp','fsr','JES','JER']))
        #
        #save_pdf(f'datacard_{tag}{self.year}_sys_dataMC.pdf')(lambda x: x())( functools.partial(self.makeplots_bydataMC, sys='ak8JERUp'))
    
    def run_makeplots(self,func):
        for p in self.process_list:
            if f'{self.bins[0]}_{p}' not in self.histos: continue 
            func(p)
    #def run_makeplots(self,func):
    #    pool = Pool()
    #    run_over_processes = [p for p in self.process_list if f'{self.bins[0]}_{p}' in self.histos]
    #    _ = pool.map(func, run_over_processes)
    #    print(_[1][0].)
    #    #print(plt.gcf().number)

    def makeplots_bydataMC(self, sys=''):
        print(self.process_list)
        for b in self.bins:
            fig, ax, ax2 = self.init_fig_axes()
            #

            edges = self.histos[f'{b}_data_obs'].edges
            data_values  = self.histos[f'{b}_data_obs'].values
            mc_values    = np.nansum([self.histos[f'{b}_{k}'].values for k in self.process_list 
                                      if f'{b}_{k}' in self.histos and 'data' not in k],axis=0)
            mc_variances = np.nansum([self.histos[f'{b}_{k}'].variances for k in self.process_list 
                                      if f'{b}_{k}' in self.histos and 'data' not in k],axis=0)
            sys_values   = np.nansum([self.histos[f'{b}_{k}_{sys}'].values for k in self.process_list 
                                      if f'{b}_{k}_{sys}' in self.histos and 'data' not in k],axis=0)
            #

            ax.errorbar(x=(edges[1:]+edges[:-1])/2, y=data_values, yerr=np.sqrt(data_values), 
                        fmt='.',
                        label='Data', color='k')

            ax2.errorbar(x=(edges[1:]+edges[:-1])/2, y=data_values/np.where(mc_values == 0, np.finfo(np.float).eps, mc_values), 
                         yerr=np.sqrt(data_values)/np.where(data_values ==0, np.finfo(np.float).eps,data_values), 
                         #yerr=np.sqrt(sumw2)/nom, 
                         fmt='.',
                         label='Data', color='k')
            #
            ax2_nom = np.where(np.append(mc_values,0) == 0, np.finfo(np.float).eps, np.append(mc_values,0))            
            ax.step( x=edges, y=np.append(mc_values, 0),        where='post', label='MC', color='blue')
            ax2.step(x=edges, y=np.append(mc_values,0)/ax2_nom, where='post',             color='blue')
            #
            ax.step( x=edges, y=np.append(sys_values, 0),         where='post', label=sys,  color='red')
            ax2.step(x=edges, y=np.append(sys_values, 0)/ax2_nom, where='post',             color='red')
            #
            variation = sys_mult[self.year][sys][b] * ((np.append(sys_values, 0)/ax2_nom) - 1 ) + 1
            ax.step( x=edges, y=ax2_nom * variation,         where='post', label=sys+f" at {sys_mult[self.year][sys][b]} sigma",  color='purple')
            ax2.step(x=edges, y=  variation, where='post',             color='purple')
            #
            self.format_fig_axes(fig, ax, ax2)
            fig.suptitle(f'{b}_Data_vs_MC_{self.year}')
            ax.set_xlim(edges[0],edges[-1])
            #
            plt.xlabel('Bins from flattened NN vs. sdM')

        #return fig
        return [plt.figure(i) for i in plt.get_fignums()]

    def makeplots_bypt(self, p):
        # get hist objects per process
        hist_names = re.findall(rf'\w*{p}\w*' ,' '.join(self.histos.keys()))
        #print(p)

        #
        colors = plt.cm.gist_rainbow(np.linspace(0,1,len(self.uncert)))
        #
        for b in self.bins:
            edges = self.histos[f'{b}_{p}'].edges
            if f'{b}_{p}' not in hist_names: continue
            get_allvalues = (lambda p_str: self.histos[f'{b}_{p_str}'].values)
            nom = self.histos[f'{b}_{p}'].values 
            sumw2 = self.histos[f'{b}_{p}'].variances
            fig, ax, ax2 = self.init_fig_axes()
            ax.errorbar(x=(edges[1:]+edges[:-1])/2, y=nom, yerr=np.sqrt(sumw2), 
                        fmt='.',
                        label=p, color='k')
            ax2.errorbar(x=(edges[1:]+edges[:-1])/2, y=np.ones_like(nom), 
                         yerr=np.sqrt(sumw2)/np.where(nom ==0, np.finfo(np.float).eps,nom), 
                         #yerr=np.sqrt(sumw2)/nom, 
                         fmt='.',
                         label=p, color='k')
            for i,unc in enumerate(self.uncert):
                #print(unc)
                if not any(re.findall(rf'{unc}[Up,Down,erdOn]',' '.join(hist_names))): continue # findall generates a list of matches so not any means unc str not found

                ax2_nom = np.where(np.append(nom,0) == 0, np.finfo(np.float).eps, np.append(nom,0))
                if unc == 'erdOn':
                    var = get_allvalues(f'{p}_{unc}') 
                    if type(var) is int: continue # this basically checks if process has this systematic
                    ax.step(x=edges, y=np.append(var,0), where='post', label='erdOn', color=colors[i])
                    ax2.step(x=edges, y=np.append(var,0)/ax2_nom, where='post', color=colors[i])
                else:
                    var_up   = get_allvalues(f'{p}_{unc}Up') 
                    var_down = get_allvalues(f'{p}_{unc}Down')
                    if type(var_up) is int or type(var_down) is int: continue # this basically checks if process has this systematic
                    ax.step(x=edges, y=np.append(var_up,0),   where='post', label=f'{unc}', color=colors[i])
                    ax.step(x=edges, y=np.append(var_down,0), where='post', color=colors[i])
                    ax2.step(x=edges, y=np.append(var_up,0)  /ax2_nom,   where='post', color=colors[i])
                    ax2.step(x=edges, y=np.append(var_down,0)/ax2_nom, where='post', color=colors[i])
            #
            self.format_fig_axes(fig, ax, ax2)
            #
            fig.suptitle(f'{b}_{p}_{self.year}')
            ax.set_xlim(edges[0],edges[-1])
            #
            plt.xlabel('Bins from flattened NN vs. sdM')

        #return fig
        return [plt.figure(i) for i in plt.get_fignums()]
            

    def makeplots_bychannel(self, p, syss=None):
        # get hist objects per process
        hist_names = re.findall(rf'\w*{p}\w*' ,' '.join(self.histos.keys()))
        #print(p)

        #
        colors = plt.cm.gist_rainbow(np.linspace(0,1,len(self.uncert)))
        #
        if syss is None: uncerts = self.uncert
        else:
            uncerts = [u for u in self.uncert for sys in syss if sys in u]
        #print(uncerts)
        for b in self.bins:
            edges = self.histos[f'{b}_{p}'].edges
            if f'{b}_{p}' not in hist_names: continue
        
            get_allvalues = (lambda p_str: self.histos[f'{b}_{p_str}'].values)
            nom = self.histos[f'{b}_{p}'].values 
            sumw2 = self.histos[f'{b}_{p}'].variances
            for i,unc in enumerate(uncerts):
                #print(unc)
                if not any(re.findall(rf'{unc}[Up,Down,erdOn]',' '.join(hist_names))): continue # findall generates a list of matches so not any means unc str not found
                #if sys not in unc: continue
        
                fig, ax, ax2 = self.init_fig_axes()
                ax2_nom = np.where(np.append(nom,0) == 0, np.finfo(np.float).eps, np.append(nom,0))
                if unc == 'erdOn':
                    var = get_allvalues(f'{p}_{unc}') 
                    if type(var) is int: continue # this basically checks if process has this systematic
                    ax.step(x=edges, y=np.append(var,0), where='post', label='erdOn', color=colors[i])
                    ax2.step(x=edges, y=np.append(var,0)/ax2_nom, where='post', color=colors[i])
                else:
                    var_up   = get_allvalues(f'{p}_{unc}Up') 
                    var_down = get_allvalues(f'{p}_{unc}Down')
                    if type(var_up) is int or type(var_down) is int: continue # this basically checks if process has this systematic
                    if sys == '':
                        ax.step(x=edges, y=np.append(var_up,0),   where='post', label=f'{unc}', color=colors[i])
                        ax.step(x=edges, y=np.append(var_down,0), where='post', color=colors[i])
                        ax2.step(x=edges, y=np.append(var_up,0)  /ax2_nom,   where='post', color=colors[i])
                        ax2.step(x=edges, y=np.append(var_down,0)/ax2_nom, where='post', color=colors[i])
                    else:
                        ax.step(x=edges, y=np.append(var_up,0),   where='post', label=f'{unc}Up',   color=colors[i])
                        ax.step(x=edges, y=np.append(var_down,0), where='post', label=f'{unc}Down', color=colors[i+5])
                        ax2.step(x=edges, y=np.append(var_up,0)  /ax2_nom,   where='post', color=colors[i])
                        ax2.step(x=edges, y=np.append(var_down,0)/ax2_nom,   where='post', color=colors[i+5])
                        #
        
                #
                ax.errorbar(x=(edges[1:]+edges[:-1])/2, y=nom, yerr=np.sqrt(sumw2), 
                            fmt='.',
                            label=p, color='k')
                ax2.errorbar(x=(edges[1:]+edges[:-1])/2, y=np.ones_like(nom), 
                             yerr=np.sqrt(sumw2)/np.where(nom ==0, np.finfo(np.float).eps,nom), 
                             #yerr=np.sqrt(sumw2)/nom, 
                             fmt='.',
                             label=p, color='k')
                #
                self.format_fig_axes(fig, ax, ax2)
                #
                fig.suptitle(f'{b}_{p}_{self.year}')
                ax.set_xlim(edges[0],edges[-1])
                #
                plt.xlabel('Bins from flattened NN vs. sdM')

        #return fig
        return [plt.figure(i) for i in plt.get_fignums()]



    def makeplots(self, p):
        # get hist objects per process
        hist_names = re.findall(rf'\w*{p}\w*' ,' '.join(self.histos.keys()))
        get_allvalues = (lambda p_str : sum([self.histos[f'{b}_{p_str}'].values    for b in self.bins if f'{b}_{p_str}' in hist_names]))
        get_sumw2      = (lambda p_str : sum([self.histos[f'{b}_{p_str}'].variances for b in self.bins if f'{b}_{p_str}' in hist_names]))
        #print(p)
        nom = get_allvalues(p) 
        sumw2 = get_sumw2(p)   
        edges = self.histos[f'{self.bins[0]}_{p}'].edges
        #
        colors = plt.cm.gist_rainbow(np.linspace(0,1,len(self.uncert)))
        #
        fig, ax, ax2 = self.init_fig_axes()
        #
        ax.errorbar(x=(edges[1:]+edges[:-1])/2, y=nom, yerr=np.sqrt(sumw2), 
                     fmt='.',
                     label=p, color='k')
        ax2.errorbar(x=(edges[1:]+edges[:-1])/2, y=np.ones_like(nom), 
                     yerr=np.sqrt(sumw2)/np.where(nom ==0, np.finfo(np.float).eps,nom), 
                     #yerr=np.sqrt(sumw2)/nom, 
                    fmt='.',
                    label=p, color='k')

        for i,unc in enumerate(self.uncert):
            #print(unc)
            ax2_nom = np.where(np.append(nom,0) == 0, np.finfo(np.float).eps, np.append(nom,0))
            #ax2_nom = np.append(nom,0)
            if unc == 'erdOn':
                var = get_allvalues(f'{p}_{unc}') 
                if type(var) is int: continue # this basically checks if process has this systematic
                ax.step(x=edges, y=np.append(var,0), where='post', label='erdOn', color=colors[i])
                ax2.step(x=edges, y=np.append(var,0)/ax2_nom, where='post', color=colors[i])
            else:
                var_up   = get_allvalues(f'{p}_{unc}Up') 
                var_down = get_allvalues(f'{p}_{unc}Down')
                if type(var_up) is int or type(var_down) is int: continue # this basically checks if process has this systematic
                ax.step(x=edges, y=np.append(var_up,0),   where='post', label=f'{unc}', color=colors[i])
                ax.step(x=edges, y=np.append(var_down,0), where='post', color=colors[i])
                ax2.step(x=edges, y=np.append(var_up,0)  /ax2_nom,   where='post', color=colors[i])
                ax2.step(x=edges, y=np.append(var_down,0)/ax2_nom, where='post', color=colors[i])
        #
        self.format_fig_axes(fig, ax, ax2)
        #
        fig.suptitle(f'{p}_{self.year}')
        ax.set_xlim(edges[0],edges[-1])
        #
        plt.xlabel('Bins from flattened NN vs. sdM')
        
        return [plt.figure(i) for i in plt.get_fignums()]

    @staticmethod
    def init_fig_axes():
        fig, (ax, ax2) = plt.subplots(2,1, sharex=True, gridspec_kw={'height_ratios':[3,1]})
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.11,
            right=0.88,
            hspace=0.0,
            wspace=0.2
        )
        return fig, ax, ax2

    @staticmethod
    def format_fig_axes(fig, ax, ax2):
        # set ax format
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='both', direction='in', top=True, right=True)
        ax.set_ylabel('Events / bins')
        ax.grid(True)
        ax.legend(ncol=2, fontsize='xx-small', framealpha = 0)
        # set ax2 format
        ax2.axhline(1, color='k', linewidth='1', linestyle='--', dashes=(4,8), snap=True)
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax2.yaxis.set_major_locator(FixedLocator([.85,1,1.15]))
        ax2.yaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(which='both', direction='in', top=True, right=True)
        ax2.set_ylim(0.70,1.30)
        ax2.set_ylabel('var/nom')
        ax2.grid(True)

if __name__ == '__main__':
    # just main
    main()
    
