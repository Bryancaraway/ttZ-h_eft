# Class to organize in dictionary form, 
# the relevent variables for run case
# i.e, handle variables sensitive to year,
# data vs mc, and jes, jer variations
import config.ana_cff as cfg
import re

class AnaVars:
    
    year   = None
    isData = False
    jec_sys= None # this should be a choice of JES up/down or JER up/down or left empty

    def __init__(self, year, isData, jec_sys=None):

        self.year     = year
        self.isData   = isData
        self.jec_sys  = ('' if jec_sys is None else jec_sys) # one of the values in jec_dict
        self.jec_type = ''  # ak4/8jer, both, ak8jms, ak8jmr
        if   'jes' in self.jec_sys:
            self.jec_type = 'both'
        elif 'jer' in self.jec_sys:
            self.jec_type = re.search(r'ak\djer' ,self.jec_sys).group()
        elif 'jms' in self.jec_sys or 'jmr' in self.jec_sys:
            #self.jec_type = self.jec_sys
            self.jec_type = 'ak8jmsr'
        #
        ud  = re.search(r'(Up|Down)', self.jec_sys)
        self.ud = ud.group() if ud else ''
        #
        self.jec_dict = {
            f'jesRelativeSample{self.year}{self.ud}' :f'_jesRelativeSample_{self.year}{self.ud}',
            f'jesHF{self.year}{self.ud}'             :f'_jesHF_{self.year}{self.ud}',
            f'jesAbsolute{self.year}{self.ud}'       :f'_jesAbsolute_{self.year}{self.ud}',
            f'jesEC2{self.year}{self.ud}'            :f'_jesEC2_{self.year}{self.ud}',
            f'jesBBEC1{self.year}{self.ud}'          :f'_jesBBEC1_{self.year}{self.ud}',
            f'jesAbsolute{self.ud}'                  :f'_jesAbsolute{self.ud}',
            f'jesEC2{self.ud}'                       :f'_jesEC2{self.ud}',
            f'jesBBEC1{self.ud}'                     :f'_jesBBEC1{self.ud}',
            f'jesHF{self.ud}'                        :f'_jesHF{self.ud}',
            f'jesRelativeBal{self.ud}'               :f'_jesRelativeBal{self.ud}',
            f'jesFlavorQCD{self.ud}'                 :f'_jesFlavorQCD{self.ud}',
            f'jesHEMIssue{self.ud}'                  :f'_jesHEMIssue{self.ud}', # only for 2018
            #
            f'ak4jer{self.ud}'                    :f'_jer{self.ud}',
            f'ak8jer{self.ud}'                    :f'_jer{self.ud}',
            #
            f'jms{self.ud}'                       :f'_jms{self.ud}',
            f'jmr{self.ud}'                       :f'_jmr{self.ud}',
            #
            '':'',
        }
        #
        self.var_fact = {
            'ak4jer': {
                'Jet_pt'           : f"Jet_pt{self.jec_dict[self.jec_sys]}" ,
                'Jet_mass'         : f"Jet_mass{self.jec_dict[self.jec_sys]}" ,
                'MET_pt'           : f"MET_T1_pt{self.jec_dict[self.jec_sys]}",
                'MET_phi'          : f"MET_T1_phi{self.jec_dict[self.jec_sys]}",
            },
            'ak8jer': {
                'FatJet_pt'        : f"FatJet_pt{self.jec_dict[self.jec_sys]}" ,
                'FatJet_msoftdrop' : f"FatJet_msoftdrop{self.jec_dict[self.jec_sys]}" , 
                'MET_pt'           : f"MET_T1_pt{self.jec_dict[self.jec_sys]}",
                'MET_phi'          : f"MET_T1_phi{self.jec_dict[self.jec_sys]}",
            },
            'both' : {
                'Jet_pt'           : f"Jet_pt{self.jec_dict[self.jec_sys]}" ,
                'Jet_mass'         : f"Jet_mass{self.jec_dict[self.jec_sys]}" ,
                'FatJet_pt'        : f"FatJet_pt{self.jec_dict[self.jec_sys]}" ,
                'FatJet_msoftdrop' : f"FatJet_msoftdrop{self.jec_dict[self.jec_sys]}" , 
                'MET_pt'           : f"MET_T1_pt{self.jec_dict[self.jec_sys]}",
                'MET_phi'          : f"MET_T1_phi{self.jec_dict[self.jec_sys]}",
            },
            'ak8jmsr' : {
                'FatJet_msoftdrop' : f"FatJet_msoftdrop{self.jec_dict[self.jec_sys]}" , 
            },
            '': {
                'MET_pt' :'MET_pt'  if not self.isData else 'MET_pt',
                'MET_phi':'MET_phi' if not self.isData else 'MET_phi',
            }
        } # need this for defualt ['']['']
        
        #
        self.q_fact = (lambda k: self.var_fact[self.jec_type].get(k,k))
        #
        #self.ak4   = cfg.ana_vars['ak4lvec']['TLVarsLC']+cfg.ana_vars['ak4vars']
        #self.ak8   = cfg.ana_vars['ak8lvec']['TLVarsLC']+cfg.ana_vars['ak8vars']+cfg.ana_vars['ak8sj']
        #self.event = cfg.ana_vars['valvars']+(['run']+cfg.ana_vars['dataHLT_all']+cfg.ana_vars[f'dataHLT_{self.year}'] if self.isData else
        #                                      cfg.ana_vars['sysvars_mc']+cfg.ana_vars[f'sysvars_{self.year}'])+(cfg.ana_vars['HEM_veto'] if self.year == '2018' else [])
        #self.gen   = cfg.ana_vars['genpvars']
        #self.RC    = cfg.ana_vars['valRCvars']

        
    def get(self, var_list):
        return map(self.q_fact, var_list)

    # magic methods #
    def __getitem__(self,key):
        #if key in self.allowed_keys: return getattr(self,key)
        #else:
        try:
            return self.var_fact[self.jec_type].get(key,key)
        except:
            print(f"Can't get key: {key}, instead setting as dummy branch")
            return 'weight'
    
 
