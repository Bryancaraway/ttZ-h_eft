import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
from modules.plotAna import Plotter, StackedHist, Hist
import operator as op
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import config.ana_cff as cfg
#
ttbar_dict = {}
w_dict     = {}
eleData_dict = {}
muData_dict = {}
lep_dict = {}
for y in cfg.Years:
    tt_dir = f'files/{y}/mc_files/TTBarLep_val.pkl'
    ele_dir = f'files/{y}/data_files/EleData_val.pkl'
    mu_dir = f'files/{y}/data_files/MuData_val.pkl'
    #tt_pd = pd.read_pickle(tt_dir)
    ele_pd = pd.read_pickle(ele_dir) 
    ele_pd = ele_pd[ele_pd['passSingleLepElec'] == 1]
    ele_pd = ele_pd[ele_pd['SAT_Pass_HEMVeto_DataOnly_drLeptonCleaned'] == True] if y == '2018' else ele_pd
    mu_pd = pd.read_pickle(mu_dir) 
    mu_pd = mu_pd[mu_pd['passSingleLepMu'] == 1]
    mu_pd = mu_pd[mu_pd['SAT_Pass_HEMVeto_DataOnly_drLeptonCleaned'] == True] if y == '2018' else mu_pd
    
    plt.hist(
        x = mu_pd['Lep_pt'][mu_pd['Pass_trigger_muon'] == 0],
        bins=50,
        range = (0,500),
        histtype='step',
        density=True,
        label = f'MuData_{y}'
    )
plt.xlabel('Mu_pt')
plt.ylabel('Normed events')
plt.title(f'Muon_Data_FAIL_trigger_muon')
plt.yscale('log')
plt.legend()
plt.show()

#    plt.hist(
#        x = ele_pd['Lep_pt'][ele_pd['Pass_trigger_electron'] == 0],
#        bins=50,
#        range = (0,500),
#        histtype='step',
#        density=True,
#        label = f'EleData_{y}'
#    )
#plt.xlabel('Ele_pt')
#plt.ylabel('Normed events')
#plt.title(f'Electron_Data_FAIL_trigger_electron')
#plt.yscale('log')
#plt.legend()
#plt.show()
exit()
        #density = True,
    #plt.hist(
    #    x = [ele_pd['Lep_pt'][ele_pd['passSingleLepElec'] == 1], ele_pd['Lep_pt'][(ele_pd['passSingleLepElec'] == 1) & (ele_pd['Pass_trigger_electron'] == 1)]],
    #    bins=50,
    #    range = (0,500),
    #    histtype='step',
    #    #density = True,
    #    label =  ['w/o Pass_trigger_electron','with Pass_trigger_electron']
    #)
    #plt.xlabel('Ele_pt')
    #plt.ylabel('Events')
    #plt.title(f'Electron_Data_{y}')
    #plt.yscale('log')
    #plt.legend()
    #plt.show()
    #plt.clf()
    #
    #plt.hist(
    #    x = [mu_pd['Lep_pt'][mu_pd['passSingleLepMu'] == 1], mu_pd['Lep_pt'][(mu_pd['passSingleLepMu'] == 1) & (mu_pd['Pass_trigger_muon'] == 1)]],
    #    bins=50,
    #    range = (0,500),
    #    histtype='step',
    #    #density = True,
    #    label =  ['w/o Pass_trigger_muon','with Pass_trigger_muon']
    #)
    #plt.xlabel('Mu_pt')
    #plt.ylabel('Events')
    #plt.title(f'Muon_Data_{y}')
    #plt.yscale('log')
    #plt.legend()
    #plt.show()
    #plt.clf()

    #x_bins= np.append(np.arange(25,200,5), [225,275,325,375,450])
    #
#    sns.regplot(
#        x=np.clip(tt_pd['Lep_pt'][tt_pd['passSingleLepElec'] == 1],0,500),
#        y=tt_pd['Stop0l_trigger_eff_Electron_pt'][tt_pd['passSingleLepElec'] == 1],
#        x_bins=x_bins,
#        x_estimator=np.mean,
#        label=f'{y}',
#        fit_reg=False,
#        scatter_kws={'s':10}
#    )
#plt.xlabel('Elec_pt')
#plt.xlim(0,500)
#plt.ylim(0,1)
#plt.legend()
#plt.ylabel('Elec_pt_eff')
#plt.title('ttbar_MC')
#plt.grid(True)
#plt.show()
#plt.clf()
#
#for y in cfg.Years:#['2018']:
#    tt_dir = f'files/{y}/mc_files/TTBarLep_val.pkl'
#    ele_dir = f'files/{y}/data_files/EleData_val.pkl'
#    mu_dir = f'files/{y}/data_files/MuData_val.pkl'
#    tt_pd = pd.read_pickle(tt_dir)
#    #ele_pd = pd.read_pickle(ele_dir) 
#    #mu_pd = pd.read_pickle(mu_dir) 
#    x_bins= np.append(np.arange(25,200,5), [225,275,325,375,450])
#    sns.regplot(
#        x=np.clip(tt_pd['Lep_pt'][tt_pd['passSingleLepMu'] == 1],0,500),
#        y=tt_pd['Stop0l_trigger_eff_Muon_pt'][tt_pd['passSingleLepMu'] == 1],
#        x_bins=x_bins,
#        x_estimator=np.mean,
#        label=f'{y}',
#        fit_reg=False,
#        scatter_kws={'s':10}
#    )
#plt.xlabel('Mu_pt')
#plt.xlim(0,500)
#plt.ylim(0,1)
#plt.legend()
#plt.ylabel('Mu_pt_eff')
#plt.title(f'ttbar_MC')
#plt.grid(True)
#plt.show()
#plt.clf()

    #ttbar_dict[y] = tt_pd
    #w_dict[y] = (tt_pd['weight']* np.sign(tt_pd['genWeight']) 
    #             * tt_pd['Stop0l_topptWeight']
    #             * (tt_pd['SAT_HEMVetoWeight'] if  y == '2018' else 1.0 )
    #             * (pd.concat([tt_pd['Stop0l_trigger_eff_Electron_pt'][tt_pd['passSingleLepElec']==1],tt_pd['Stop0l_trigger_eff_Muon_pt'][tt_pd['passSingleLepMu']==1]]).sort_index())
    #             * tt_pd['BTagWeight'] 
    #             * tt_pd['puWeight']  
    #             * (tt_pd['PrefireWeight'] if y != '2018' else 1.0)
    #             * tt_pd['ISRWeight'])
    #w_dict[y] = w_dict[y]/(w_dict[y].sum())
    ##
    #eleData_dict[y] = ele_pd[['Lep_pt','Lep_eta','Lep_phi']][((ele_pd['passSingleLepElec'] == 1) & (ele_pd['run']>=319077) & ele_pd['SAT_Pass_HEMVeto_DataOnly'])]
    #muData_dict[y] = mu_pd  [['Lep_pt','Lep_eta','Lep_phi']][((mu_pd['passSingleLepMu'] == 1)    & (mu_pd['run']>=319077)  &  mu_pd['SAT_Pass_HEMVeto_DataOnly'])]
    #lep_dict[y] = pd.concat([eleData_dict[y],muData_dict[y]])
    #print(lep_dict[y])

exit()

plt.hist2d(
    x=lep_dict['2018']['Lep_eta'],
    y=lep_dict['2018']['Lep_phi'],
    bins=[30,18],
    range=((-3,3),(-3.15,3.15)),
    cmin=1,
    norm=colors.LogNorm()
)
plt.xlabel('Lep_eta')
plt.ylabel('Lep_phi')
plt.title('Data_postHEM')
plt.colorbar()
plt.show()
plt.clf()
exit()

plt.hist(
    x = [muData_dict[y] for y in cfg.Years],
    bins=20,
    range = (0,500),
    histtype='step',
    weights = [np.ones(len(muData_dict[y]))/len(muData_dict[y]) for y in cfg.Years],
    label =  ['MuData_'+y for y in cfg.Years]
    )
plt.xlabel('Mu_pt')
plt.legend()
plt.show()
plt.clf()

plt.hist(
    x = [eleData_dict[y] for y in cfg.Years],
    bins=20,
    range = (0,500),
    histtype='step',
    weights = [np.ones(len(eleData_dict[y]))/len(eleData_dict[y]) for y in cfg.Years],
    label =  ['EleData_'+y for y in cfg.Years]
    )
plt.xlabel('Ele_pt')
plt.legend()
plt.show()
plt.clf()

plt.hist(
    x = [ttbar_dict[y] for y in cfg.Years],
    bins=20,
    range = (0,500),
    histtype='step',
    weights = [w_dict[y] for y in cfg.Years],
    label =  ['ttbar_'+y for y in cfg.Years]
    )
plt.xlabel('Lep_pt')
plt.legend()
plt.show()
plt.clf()
    


