# Class to organize in dictionary form, 
# the relevent variables for run case
# i.e, handle variables sensitive to year,
# data vs mc, and jes, jer variations
import config.ana_cff as cfg

class AnaVars:
    
    year   = None
    isData = False
    jec_sys= None # this should be a choice of JES up/down or JER up/down or left empty
    jec_dict = { 'norm': { '': '', # to query the right string with jec variation
                           'JESUp': '_jesTotalUp',
                           'JESDown': '_jesTotalDown',
                           'JERUp': '_jerUp',
                           'JERDown': '_jerDown'},
                 'cap': { '': '',
                          'JESUp': '_JESUp',
                          'JESDown': '_JESDown',
                          'JERUp': '_JERUp',
                          'JERDown': '_JERDown'}
                 }
    #rtc_dict = { 'ResolvedTopCandidate_JERUp_discriminator': 'ResolvedTopCandidate_JESUp_discriminator',
    #             'ResolvedTopCandidate_JERDown_discriminator': 'ResolvedTopCandidate_JESDown_discriminator'}
    
    LC = '_drLeptonCleaned'
    #
    allowed_keys = ['ak4','ak8','event','gen','RC']

    def __init__(self, year, isData, jec_sys=None, jec_type=None):

        self.year    = year
        self.isData  = isData
        self.jec_sys = ('' if jec_sys is None else jec_sys)
        self.jec_type = ('' if jec_type is None else jec_type) # ak4 or ak8 // or sj for sdM (not implemented yet)
        #
        #q_rtc = (lambda k: self.rtc_dict.get(k,k)) # i dont have jer variations on the rtc candidate score so..
        q_rtc  = (lambda k: k.replace('JER','JES') if 'JER' in k else k) 
        
        self.var_fact = {
            'ak4': {
                'Jet_pt'             : f"Jet_pt{self.jec_dict['norm'][self.jec_sys]}" ,
                'Jet_mass'           : f"Jet_mass{self.jec_dict['norm'][self.jec_sys]}" ,
                f'Jet_pt{self.LC}'   : f"Jet_pt{self.jec_dict['norm'][self.jec_sys]}{self.LC}" ,
                f'Jet_mass{self.LC}' : f"Jet_mass{self.jec_dict['norm'][self.jec_sys]}{self.LC}" ,
                f'Jet_lepcleaned_idx{self.LC}' : f"Jet_lepcleaned_idx{self.jec_dict['norm'][self.jec_sys]}{self.LC}",
                
                'MET_pt'             : f"MET_pt{self.jec_dict['norm'][self.jec_sys]}",
                f'nJets30{self.LC}'  : f"nJets30{self.jec_dict['cap'][self.jec_sys]}{self.LC}",
                f'nBottoms{self.LC}' : f"nBottoms{self.jec_dict['cap'][self.jec_sys]}{self.LC}",
                'ResolvedTopCandidate_discriminator': (q_rtc)(f"ResolvedTopCandidate{self.jec_dict['cap'][self.jec_sys]}_discriminator"),
                'ResolvedTopCandidate_j1Idx'        : (q_rtc)(f"ResolvedTopCandidate{self.jec_dict['cap'][self.jec_sys]}_j1Idx"),
                'ResolvedTopCandidate_j2Idx'        : (q_rtc)(f"ResolvedTopCandidate{self.jec_dict['cap'][self.jec_sys]}_j2Idx"),
                'ResolvedTopCandidate_j3Idx'        : (q_rtc)(f"ResolvedTopCandidate{self.jec_dict['cap'][self.jec_sys]}_j3Idx")
            },
            'ak8': {
                f'FatJet_pt{self.LC}'       : f"FatJet_pt{self.jec_dict['norm'][self.jec_sys]}{self.LC}" ,
                f'FatJet_msoftdrop{self.LC}': f"FatJet_msoftdrop{self.jec_dict['norm'][self.jec_sys]}{self.LC}" , # might change to sj
            },
            '':{'':''}} # need this for defualt ['']['']
        
        #
        self.q_fact = (lambda k: self.var_fact[self.jec_type].get(k,k))
        #
        self.ak4   = cfg.ana_vars['ak4lvec']['TLVarsLC']+cfg.ana_vars['ak4vars']
        self.ak8   = cfg.ana_vars['ak8lvec']['TLVarsLC']+cfg.ana_vars['ak8vars']+cfg.ana_vars['ak8sj']
        self.event = cfg.ana_vars['valvars']+(['run']+cfg.ana_vars['dataHLT_all']+cfg.ana_vars[f'dataHLT_{self.year}'] if self.isData else
                                              cfg.ana_vars['sysvars_mc']+cfg.ana_vars[f'sysvars_{self.year}'])+(cfg.ana_vars['HEM_veto'] if self.year == '2018' else [])
        self.gen   = cfg.ana_vars['genpvars']
        self.RC    = cfg.ana_vars['valRCvars']

        
    def get(self, var_list):
        return map(self.q_fact, var_list)

    # magic methods #
    def __getitem__(self,key):
        if key in self.allowed_keys: return getattr(self,key)
        else:
            try:
                return self.var_fact[self.jec_type].get(key,key)
            except:
                print(f"Can't get key: {key}, instead setting as dummy branch")
                return 'weight'
    
 
