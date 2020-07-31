import uproot
import sys
import numpy as np
import awkward
import concurrent.futures
from collections import defaultdict
#import functools
from modules.AnaDict import AnaDict
from modules.getdata import getData
from lib.fun_library import clop_pear_ci, deltaR, fillne
import cfg.deepsleepcfg as cfg

year = sys.argv[1]

def pickle_mcdata():

    roodir  = '/cms/data/store/user/ttxeft/Skim_nanoAOD/lep_study_files/' 
    
    n_events = AnaDict()

    mu_id_map = {0:'loose',1:'med',2:'tight'}
    el_id_map = {0:'fail', 1:'veto',2:'loose',3:'med',4:'tight'}

    
    def mc_ana():
        roofile = f'{roodir}MC_{year}_lep.root'
        with uproot.open(roofile) as f:
            t = {'ttzh':f.get('Training_bb/TTZH'),
                 'tt'  :f.get('Training_bb/TTBarLep')}
            for key in t.keys():
                mu_data, el_data, gen_data, weight = get_tree_data(t[key], getgen=True)
                calc_lepSel_eff(mu_data,el_data, key, gen_data)
                #n_events[key] = calc_n_wLepSel(mu_data, el_data, weight)
                #del mu_data, el_data, weight
            ##
            exit()
            s_over_sqrt_b = {k : n_events['ttzh'][k]['sum_w']/np.sqrt(n_events['tt'][k]['sum_w']) for k in n_events['tt']}
            print(f'\n\nTTZH over sqrt TTBarLep, {year} MC\n')
            [print(f'{k:100}\t s/sqrt(b): {v}') for k,v in sorted(s_over_sqrt_b.items(), key=lambda item: item[1], reverse=True)]         
    #
    def data_ana():
        roofile = f'{roodir}Data_{year}_lep.root'
        with uproot.open(roofile) as f:
            t = {'el' : f.get('Training_bb/EleData'),
                 'mu' : f.get('Training_bb/MuData')}
            for key in t.keys():
                mu_data, el_data, weight = get_tree_data(t[key])
                n_events[key] = calc_n_wLepSel(mu_data, el_data, weight, opt=key)
                del mu_data, el_data, weight
            n_data = {k : n_events['el'][k]['sum_w']+n_events['mu'][k]['sum_w'] for k in n_events['el']}
            print(f'\n\n{year} Data Yield\n')
            [print(f"{k:100}\t Yield: {v:10}\t Ele Yield: {n_events['el'][k]['sum_w']}\t Mu Yield: {n_events['mu'][k]['sum_w']}") for k,v in sorted(n_data.items(), key=lambda item: item[1], reverse=True)]

    mc_ana()
    data_ana()
    n_events.to_pickle(f'mcdata_lepsel_{year}.pkl')

def calc_lepSel_eff(mu, el, mc, gen):
    event_mask = getData.DF_Container.event_mask
    gen_mu = gen[abs(gen['GenPart_pdgId']) == 13]
    gen_el = gen[abs(gen['GenPart_pdgId']) == 11]
    def get_gen_matched(lep, gen, lep_type):
        gen_pt, gen_eta, gen_phi, gen_mother = map(fillne, [gen['GenPart_pt'], gen['GenPart_eta'], gen['GenPart_phi'], gen['GenPart_mother']])
        lep_gen_match = []
        for i in range(gen_pt.shape[1]):
            #lep_gen_match.append((deltaR(lep[f'{lep_type}_eta'], lep[f'{lep_type}_phi'], gen_eta[:,i], gen_phi[:,i]) < 0.1) & (abs((lep[f'{lep_type}_pt']/gen_pt[:,i]) - 1) < .30))
            lep_gen_match.append((deltaR(lep[f'{lep_type}_eta'], lep[f'{lep_type}_phi'], gen_eta[:,i], gen_phi[:,i]) > 0.1) | (abs((lep[f'{lep_type}_pt']/gen_pt[:,i]) - 1) > .30))
        no_match = (sum(lep_gen_match) > 0)
        gen_nom = gen[abs(gen['GenPart_eta']) < 1.2][no_match.sum()>0]
        lep_nom = lep[no_match][lep[f'{lep_type}_eta'][no_match].counts > 0]
        lep_nom = lep_nom[gen_nom['GenPart_pt'].counts>0]
        gen_nom = gen_nom[gen_nom['GenPart_pt'].counts>0]
        for i in range(10):
            print(f"Event {i+1}: GEN Info: Gen pt:{gen_nom['GenPart_pt'][i]}, Gen Eta:{gen_nom['GenPart_eta'][i]}, Gen Phi:{gen_nom['GenPart_phi'][i]}, Gen Mom:{gen_nom['GenPart_mother'][i]}")
            print(f"Event {i+1}: RECO Info:{lep_type} pt:{lep_nom[f'{lep_type}_pt'][i]}, {lep_type} Eta:{lep_nom[f'{lep_type}_eta'][i]}, {lep_type} Phi:{lep_nom[f'{lep_type}_phi'][i]}, {lep_type} Id:{lep_nom[f'{lep_type}_FlagId'][i]}, {lep_type} miniIso:{lep_nom[f'{lep_type}_miniPFRelIso_all'][i]}")
        exit()
        return lep[(sum(lep_gen_match) > 0)]
        
    mu = get_gen_matched(mu, gen_mu, 'Muon')
    print(mu[mu['Muon_pt'].counts > 0])
    exit()
    el = get_gen_matched(el, gen_el, 'Electron')
    #
    mu_num = mu[cfg.lep_sel['muon'](mu)][event_mask] # med id, .2 miniIso
    mu_den = mu[mu['Muon_pt'] > 30][event_mask]
    #
    el_num = el[((cfg.lep_sel['electron'][year](el)))][event_mask] # tight id, .1 miniIso
    el_den = el[((el['Electron_cutBasedNoIso'] >= 2) & (el['Electron_pt'] > (30 if year == '2016' else 35)))][event_mask]
    #
    bins = {'pt' : { 'Muon'    : np.append(np.linspace(30,100,5+1),[125,150,175,200,250,300,400,500]),
                     'Electron': np.append(np.linspace((30 if year == '2016' else 35),100,5+1),[125,150,175,200,250,300,400,500])
                 },
            'eta': {'Muon'     : [-2.6,-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,-0.2,
                                  0.2,0.3,0.9,1.2,1.6,2.1,2.4,2.6],
                    'Electron' : [-2.6, -2.5, -1.566, -1.444, -0.8, 0.0, 0.8, 1.444, 1.566, 2.5, 2.6] 
                }
        }
    #
    import matplotlib.pyplot as plt
    from matplotlib.ticker import AutoMinorLocator
    #
    def plot_eff(num, den, lep):
        fig, ax = plt.subplots(1,2, figsize=(16,10))

        for i,k in enumerate(['pt','eta']):
                n_num, _     = np.histogram(num[f'{lep}_{k}'].flatten(), bins=bins[k][lep])
                n_den, edges = np.histogram(den[f'{lep}_{k}'].flatten(), bins=bins[k][lep])
                bin_c = (edges[1:]+edges[:-1])/2
                bin_w = (edges[1:]-edges[:-1])
                lo, hi = clop_pear_ci(n_num,n_den,return_error=True)
                y = n_num/n_den
                yerr = [lo,hi]
                ax[i].errorbar(x=bin_c, xerr=bin_w/2,
                               y=y, yerr=yerr, fmt='.', label=f'{mc}_{year}')
                ax[i].axhline(1, color='r', linewidth='1', linestyle='--', dashes=(4,8), snap=True)
                ax[i].grid(True)
                ax[i].set_ylim(0.00,1.05)
                ax[i].set_xlabel(f"{lep} {k}{' (GeV)' if k == 'pt' else ''}")

        fig.suptitle(f'{lep} Selection Efficiency, {mc}_{year}')
        plt.show()
        plt.clf()
        plt.close()

    def plot_2d_eff(num,den, lep):
        fig, ax = plt.subplots()
        n_num, *_ = np.histogram2d(x=num[f'{lep}_eta'].flatten(), y=num[f'{lep}_pt'].flatten(),
                                   bins=[bins['eta'][lep],bins['pt'][lep]], range=((-2.6,2.6),(0,500)))
        n_den, xedges, yedges = np.histogram2d(x=den[f'{lep}_eta'].flatten(), y=den[f'{lep}_pt'].flatten(),
                                               bins=[bins['eta'][lep],bins['pt'][lep]], range=((-2.6,2.6),(0,500)))
        x,y = np.meshgrid(xedges,yedges)
        data = n_num/n_den
        plt.pcolormesh(x,y,data.T, cmap=plt.get_cmap('viridis',16), vmin=0.60, vmax=1.0)
        plt.colorbar()
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        plt.ylim(0,500)
        plt.ylabel('pT (GeV)')
        plt.xlabel('eta')
        plt.title(f'{mc}_{year} {lep} Sel Eff')
        plt.show()
        plt.clf()
        plt.close()
        #
    #
    #plot_2d_eff(el_num,el_den,'Electron')
    #plot_2d_eff(mu_num,mu_den,'Muon')
        
    plot_eff(el_num,el_den,'Electron')
    #plot_eff(mu_num,mu_den,'Muon')

def get_tree_data(tree, getgen=False):
    executor = concurrent.futures.ThreadPoolExecutor()
    getData.DF_Container.set_attr(True, year, 2, 99, 1, None,None)
    getData.DF_Container.set_current_tree_mask(tree)
    #
    mu_keys   = ["Muon_pt","Muon_eta","Muon_phi","Muon_mass",
                 "Muon_miniPFRelIso_all","Muon_pfRelIso04_all",
                 "Muon_FlagId"]  # 0,1,2 = loose, med, tight
    el_keys = ["Electron_pt","Electron_eta","Electron_phi","Electron_mass",
               "Electron_miniPFRelIso_all", 
               "Electron_cutBasedNoIso", "Electron_cutBased"] # 0,1,2,3,4 = fail, veto, loose, med, tight
    #
    mu_data = AnaDict({k:tree.array(k, executor=executor) for k in mu_keys})
    el_data = AnaDict({k:tree.array(k, executor=executor) for k in el_keys})
    mc_weight = tree.array('weight',   executor=executor)
    g_weight  = tree.array('genWeight',executor=executor) if b'genWeight' in tree.keys() else 1
    weight = mc_weight*np.sign(g_weight)
    #
    mu_data = mu_data[(mu_data['Muon_pt']     >= 30) & (abs(mu_data['Muon_eta'])     < 2.4)]
    el_data = el_data[(el_data['Electron_pt'] >= 30) & (abs(el_data['Electron_eta']) < 2.5)]
    if getgen:
        gen_info = AnaDict({k:tree.array(k, executor=executor) for k in cfg.ana_vars['genpvars']})
        gen_info['GenPart_mother'] = gen_info['GenPart_pdgId'][gen_info['GenPart_genPartIdxMother']]
        return mu_data, el_data, gen_info, weight
    #
    return mu_data, el_data, weight

def calc_n_wLepSel(mu, el, w, opt=None):
    #
    event_mask = getData.DF_Container.event_mask 
    njets_ak4  = getData.DF_Container.n_ak4jets
    mu_id_map = {0:'loose',1:'med',2:'tight'}
    el_id_map = {0:'fail', 1:'veto',2:'loose',3:'med',4:'tight'}
    lep_sel_yield = AnaDict()
    for i in [0.2,0.4]:
        for j in [0.1,0.15]:
            for k in [1,2]:
                for l in [3,4]:
                    lep_config = f'Muon miniIso: {i}, Muon Id: {mu_id_map[k]:5}, Ele miniIso: {j:5}, Ele Id(noIso): {el_id_map[l]:5}'
                    mask, lep_pt, lep_eta = test_lep_sel((mu["Muon_miniPFRelIso_all"] < i),
                                                         (mu['Muon_FlagId'] >= k ),
                                                         (el['Electron_miniPFRelIso_all'] < j), 
                                                         (el['Electron_cutBasedNoIso'] >= l),
                                                         mu,el,opt=opt)
                    lep_sel_yield[lep_config] = {'sum_w'  : w[mask & event_mask].sum(),
                                                 'weight' : w[mask & event_mask],
                                                 'Lep_pt' : lep_pt[mask & event_mask].flatten(),
                                                 'Lep_eta': lep_eta[mask & event_mask].flatten(),
                                                 'njets'  : njets_ak4[mask[event_mask]]}
                            
    for i in [0.15, 0.25]:
        for k in [1,2]:
            for l in [3,4]:
                lep_config = f'Muon Iso: {i}, Muon Id: {mu_id_map[k]:5}, Ele Id(Iso): {el_id_map[l]:5}'
                mask, lep_pt, lep_eta  = test_lep_sel((mu["Muon_pfRelIso04_all"] < i),
                                                      (mu['Muon_FlagId'] >= k ),
                                                      (el['Electron_miniPFRelIso_all'] < 99999), 
                                                      (el['Electron_cutBased'] >= l),
                                                      mu,el,opt=opt)
                lep_sel_yield[lep_config] = {'sum_w'  : w[mask & event_mask].sum(),
                                             'weight' : w[mask & event_mask],
                                             'Lep_pt' : lep_pt[mask & event_mask].flatten(),
                                             'Lep_eta': lep_eta[mask & event_mask].flatten(),
                                             'njets'  : njets_ak4[mask[event_mask]]}

    return lep_sel_yield


def test_lep_sel(mu_iso_cut, mu_id_cut, el_iso_cut, el_id_cut, mu_data, el_data, opt=None):
    mu_mask = ((mu_iso_cut)  & (mu_id_cut))
    mu_data = mu_data[mu_mask]
    el_mask = ((el_iso_cut)  & (el_id_cut))
    el_data = el_data[el_mask]
    #
    base_mask = (mu_mask[mu_mask].counts + el_mask[el_mask].counts == 1)
    #lep_pt = mu_data['Muon_pt'][base_mask] 
    lep_pt   = awkward.JaggedArray.concatenate([mu_data['Muon_pt'], el_data['Electron_pt']],axis=1)
    lep_eta  = awkward.JaggedArray.concatenate([mu_data['Muon_eta'], el_data['Electron_eta']],axis=1)

    #
    opt = '' if opt is None else opt
    mask_dict = {''  : base_mask,
                 'el': (base_mask  & (el_mask[el_mask].counts == 1)),
                 'mu': (base_mask & (mu_mask[mu_mask].counts == 1))
    }
    oneLep_mask = mask_dict[opt]
    return oneLep_mask, lep_pt, lep_eta

#
def plot_pickle():
    
    f_name = f'mcdata_lepsel_{year}.pkl'
    a_dict = AnaDict.read_pickle(f_name)
    #print(a_dict)
    miniIso_keys = [k for k in a_dict['el'] if 'mini'  in k]
    iso_keys     = [k for k in a_dict['el'] if '(Iso)' in k]
    #
    iso_dict  = processKeys(iso_keys)
    miso_dict = processKeys(miniIso_keys)
    # make plots for miniIso
    import matplotlib.pyplot as plt

    def plot_from_dict_miso(dict_, title):
        fig, ax = plt.subplots(2,2, figsize=(8,6))
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.11,
            right=0.88,
            hspace=0.4,
            wspace=0.2
        )
        for i,m_id in enumerate(miso_dict['Muon Id']):
            for j,el_id in enumerate(miso_dict['Ele Id(noIso)']):
                val = []
                for m_iso in miso_dict['Muon miniIso']:
                    val.append([dict_[f'Muon miniIso: {float(m_iso)}, Muon Id: {m_id:5}, Ele miniIso: {float(el_iso):5}, Ele Id(noIso): {el_id:5}'] for el_iso  in miso_dict['Ele miniIso']])
                #
                ax[i,j].imshow(val, cmap="YlGn")
                #ax[i,j].table(val)
                ax[i,j].set_yticks(np.arange(len(miso_dict['Muon miniIso'])))
                ax[i,j].set_yticklabels(miso_dict['Muon miniIso'])
                ax[i,j].set_ylabel('Muon miniIso')
                #
                ax[i,j].set_xticks(np.arange(len(miso_dict['Ele miniIso'])))
                ax[i,j].set_xticklabels(miso_dict['Ele miniIso'])
                ax[i,j].set_xlabel('Ele miniIso')
                ax[i,j].set_title(f'Muon Id: {m_id}, Ele Id(noIso): {el_id}')
                #
                for l in range(len(miso_dict['Muon miniIso'])):
                    for m in range(len(miso_dict['Ele miniIso'])):
                        ax[i,j].text(m,l, f'{np.array(val)[l,m]:6.2f}', ha="center", va="center", color="k")
                
        fig.suptitle(title)
        #plt.show()
    #
    def plot_from_dict_iso(dict_, title):
        fig, ax = plt.subplots(1,2)
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.11,
            right=0.88,
            hspace=0.2,
            wspace=0.4
        )
        for i,el_id in enumerate(iso_dict['Ele Id(Iso)']):
            val = []
            for m_id in iso_dict['Muon Id']:
                val.append([dict_[f'Muon Iso: {float(m_iso)}, Muon Id: {m_id:5}, Ele Id(Iso): {el_id:5}'] for m_iso  in iso_dict['Muon Iso']])
            #
            ax[i].imshow(val, cmap="YlGn")
            ax[i].set_yticks(np.arange(len(iso_dict['Muon Id'])))
            ax[i].set_yticklabels(iso_dict['Muon Id'])
            ax[i].set_ylabel('Muon Id')
            #
            ax[i].set_xticks(np.arange(len(iso_dict['Muon Iso'])))
            ax[i].set_xticklabels(iso_dict['Muon Iso'])
            ax[i].set_xlabel('Muon Iso')
            ax[i].set_title(f'Ele Id(Iso): {el_id}')
            #
            for l in range(len(iso_dict['Muon Id'])):
                for m in range(len(iso_dict['Muon Iso'])):
                    ax[i].text(m,l, f'{np.array(val)[l,m]:6.2f}', ha="center", va="center", color="k")
                
        fig.suptitle(title)
        #plt.show()
                
    #
    mifunc = plot_from_dict_miso
    ifunc  = plot_from_dict_iso

    s_over_sqrt_b = {k : a_dict['ttzh'][k]/np.sqrt(a_dict['tt'][k]) for k in a_dict['tt']}
    n_data = {k : a_dict['el'][k]+a_dict['mu'][k] for k in a_dict['el']}

    to_do = [[a_dict['tt'], 'ttbar Yield'],
             [a_dict['ttzh'], 'ttZ/h Yield'], 
             [s_over_sqrt_b, 'Sig/sqrt(BkG)'],
             [a_dict['el'], 'Ele Data Yield'],
             [a_dict['mu'], 'Muon Data Yield'],
             [n_data, 'Total Data Yield']]
    for a,b in to_do:
        mifunc(a,b)
        ifunc(a,b)
    #
    import matplotlib.backends.backend_pdf as matpdf
    pdf = matpdf.PdfPages(f'pdf/lep_sel_criteria.pdf')
    for fig_ in range(1, plt.gcf().number+1):
        pdf.savefig( fig_ )
    pdf.close()
    plt.close('all')
    
    
def processKeys(keys):
    dict_ = defaultdict(list)
    for key in keys:
        i_k = key.split(',')
        for i in i_k:
            k,v = i.split(':')
            dict_[k.strip()].append(v.strip())
    dict_ = {k:list(set(v)) for k,v in dict_.items()}
    return dict_
            

def plot_sig_over_bkg():
    f_name = f'mcdata_lepsel_{year}.pkl'
    a_dict = AnaDict.read_pickle(f_name)
    miniIso_keys = [k for k in a_dict['el'] if 'mini'  in k]
    iso_keys     = [k for k in a_dict['el'] if '(Iso)' in k]
    # tt is bkg, ttzh is sig
    import matplotlib.pyplot as plt
    pt_bins = [30,40,50,60,120,200,300,500]

    def compare_iso(n_jets=None):
        if n_jets is None: 
            cut = (lambda x: x>=2)
            add_to_title = ''
        else: 
            cut = (lambda x: x==n_jets)
            add_to_title = f'_njets=={n_jets}'
        
        fig, ax = plt.subplots(1,2, figsize=(14,8))
        fig.subplots_adjust(
            top=0.88,
            bottom=0.11,
            left=0.11,
            right=0.88,
            hspace=0.0,
            wspace=0.4
        )
        for i,k in enumerate(miniIso_keys+iso_keys):
            #if (('0.1,' not in k or '0.2,' not in k) and '(noIso)' in k) or (('(Iso): tight' not in k or '0.15' not in k) and '(Iso)' in k) : continue
            if (('0.1,' not in k or '0.2,' not in k) and '(noIso)' in k) or (('0.15' not in k) and '(Iso)' in k) : continue

            njet_cut = (lambda y: cut(y[k]['njets']))
            #sig = a_dict['ttzh'][k]['Lep_pt'][njet_cut(a_dict['ttzh'])]
            #bkg = a_dict['tt'][k]['Lep_pt'][njet_cut(a_dict['tt'])]
            sig = a_dict['ttzh'][k]['njets'][njet_cut(a_dict['ttzh'])]
            bkg = a_dict['tt'][k][  'njets'][njet_cut(a_dict['tt'])]
            sig_w = a_dict['ttzh'][k]['weight'][njet_cut(a_dict['ttzh'])]
            bkg_w = a_dict['tt'][k]['weight'][njet_cut(a_dict['tt'])]

            #n_sig, edges = np.histogram(np.clip(sig,pt_bins[0], pt_bins[-1]), bins=pt_bins, range=(pt_bins[0], pt_bins[-1]), weights=  sig_w)
            #n_bkg, _     = np.histogram(np.clip(bkg,pt_bins[0], pt_bins[-1]), bins=pt_bins,   range=(pt_bins[0], pt_bins[-1]), weights=bkg_w)
            n_sig, edges = np.histogram(np.clip(sig, 2, 9), bins=[2,3,4,5,6,7,8,9],   range=(2,9), weights=  sig_w)
            n_bkg, _     = np.histogram(np.clip(bkg, 2, 9), bins=[2,3,4,5,6,7,8,9],   range=(2,9), weights=bkg_w)
            bin_c = (edges[1:]+edges[:-1])/2
            bin_w = (edges[1:]-edges[:-1])
            y    = n_sig/np.sqrt(n_bkg)
            #yerr = abs(y) * np.sqrt( np.power( np.sqrt() ,2) ) 
            ax[0].errorbar(x=bin_c, xerr=bin_w/2,
                           y=y, #yerr=yerr
                           fmt='.', label=k)
            #
        #ax[0].set_xlabel('Lep pt (GeV)')
        ax[0].set_xlabel('n_jets')
        #ax[0].set_xscale('log')
        ax[0].set_ylabel('Sig/sqrt(Bkg)')
        ax[0].grid(True)
        ax[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax[0].set_title(f'Lep Selection Significance, {year}{add_to_title}')
        ax[1].axis('off')
        plt.show()


    compare_iso()
    
    #import matplotlib.backends.backend_pdf as matpdf
    #pdf = matpdf.PdfPages(f'pdf/lep_sel_criteria.pdf')
    #for fig_ in range(1, plt.gcf().number+1):
    #    pdf.savefig( fig_ )
    #pdf.close()
    #plt.close('all')
        

if __name__ == '__main__':
    #
    pickle_mcdata()
    #
    #plot_pickle()
    #
    #plot_sig_over_bkg()
    

