### ######## ###
## Config for ##
# jme handling #
## ########## ##
#import subprocess as sb
import numpy as np
from config.ana_cff import cdir
# get current working directory according to git
#_wdir = sb.check_output('echo $(git rev-parse --show-toplevel)', shell=True).decode().strip('\n')+'/DeepSleep/'
_wdir, _cdir = cdir()
'''
For now this serves as just
a file referencing map
for ak8 jme correction module
'''

jme_dir  = _wdir+'data'

jec_vc   = {'2016':'Summer16_07Aug2017_V11_MC', # good 
            '2017':'Fall17_17Nov2017_V32_MC',   # good 
            '2018':'Autumn18_V19_MC'            # good 
        }
jer_vc   = {'2016':'Summer16_25nsV1b_MC',       # good 
            '2017':'Fall17_V3b_MC',             # V3 (not V3b...) 
            '2018':'Autumn18_V7b_MC'            # V8 (or V7, not V7b)
        }


ak8_jme_files = {'jec'  : (lambda year : 
                           [f'{jme_dir}/jec/{jec_vc[year]}/{jec_vc[year]}_L1FastJet_AK8PFPuppi.jec.txt',
                            f'{jme_dir}/jec/{jec_vc[year]}/{jec_vc[year]}_L2Relative_AK8PFPuppi.jec.txt',
                            f'{jme_dir}/jec/{jec_vc[year]}/{jec_vc[year]}_L3Absolute_AK8PFPuppi.jec.txt',
                            f'{jme_dir}/jec/{jec_vc[year]}/{jec_vc[year]}_L1FastJet_AK4PFPuppi.jec.txt',
                            f'{jme_dir}/jec/{jec_vc[year]}/{jec_vc[year]}_L2Relative_AK4PFPuppi.jec.txt',
                            f'{jme_dir}/jec/{jec_vc[year]}/{jec_vc[year]}_L3Absolute_AK4PFPuppi.jec.txt']),
                 'junc' : (lambda year : 
                           [f'{jme_dir}/jec/{jec_vc[year]}/{jec_vc[year]}_Uncertainty_AK8PFPuppi.junc.txt',
                            f'{jme_dir}/jec/{jec_vc[year]}/{jec_vc[year]}_Uncertainty_AK4PFPuppi.junc.txt']),
                 'jer'  : (lambda year : 
                           [f'{jme_dir}/jer/{jer_vc[year]}/{jer_vc[year]}_PtResolution_AK8PFPuppi.jr.txt',
                            f'{jme_dir}/jer/{jer_vc[year]}/{jer_vc[year]}_PtResolution_AK4PFPuppi.jr.txt']),
                 'jersf': (lambda year : 
                           [f'{jme_dir}/jer/{jer_vc[year]}/{jer_vc[year]}_SF_AK8PFPuppi.jersf.txt',
                            f'{jme_dir}/jer/{jer_vc[year]}/{jer_vc[year]}_SF_AK4PFPuppi.jersf.txt']) 
             }


puppicorr_gen = { '2016':(lambda x : 1.0062610-1.0616051*np.power(x*0.079990008,-1.2045377)),
                  '2017':(lambda x : 1+0.0090283*np.power(x,(-2*(0.0099852)))-7.30123*np.power(x,-1)),
                  '2018':(lambda x : 1+0.0231069*np.power(x,(-2*(0.0704761)))-8.42917*np.power(x,-1))  
              }

puppicorr_reco0eta1p3 = {'2016':(lambda x: 1.0930198-0.00015006789*x+(3.4486612e-07)*np.power(x,2)+(-2.6810031e-10)*np.power(x,3)+(8.6744023e-14)*np.power(x,4)+(-1.0011359e-17)*np.power(x,5)),
                         '2017':(lambda x: 1.04323417805+(8.20581677106e-05)*x+(-2.23790959145e-08)*np.power(x,2)+(-5.56816212196e-12)*np.power(x,3)+(-2.42702058503e-17)*np.power(x,4)+(5.23731618031e-19)*np.power(x,5)),
                         '2018':(lambda x: 1.06263222851+(2.97332221436e-05)*x+(-7.31785851974e-09)*np.power(x,2)+(2.53798731754e-13)*np.power(x,3)+(1.68571767997e-16)*np.power(x,4)+(-6.77123709437e-20)*np.power(x,5))
                     }
puppicorr_reco1p3eta2p5 = {'2016':(lambda x: 1.2721152-0.00057164036*x+(8.3728941e-07)*np.power(x,2)+(-5.2043320e-10)*np.power(x,3)+(1.4537521e-13)*np.power(x,4)+(-1.5038869e-17)*np.power(x,5)),
                           '2017':(lambda x: 1.11549406241+(-2.01425972518e-05)*x+(8.36181961894e-09)*np.power(x,2)+(4.39451437171e-11)*np.power(x,3)+(1.04302756829e-14)*np.power(x,4)+(-2.10404344784e-17)*np.power(x,5)),
                           '2018':(lambda x: 1.11889161475+(2.68579882197e-05)*x+(-4.30234840782e-09)*np.power(x,2)+(8.27000377942e-12)*np.power(x,3)+(1.45823446695e-15)*np.power(x,4)+(-1.65979484436e-18)*np.power(x,5))
                     }

puppicorr_massReso_0eta1p3   = (lambda x: 1.0927351+(4.1426227e-05)*x+(-1.3736806e-07)*np.power(x,2)+(1.2295818e-10)*np.power(x,3)+(-4.1970754e-14)*np.power(x,4)+(4.9237927e-18)*np.power(x,5))
puppicorr_massReso_1p3eta2p5 = (lambda x: 1.1649278+(-0.00012678903)*x+(1.0594037e-07)*np.power(x,2)+(6.0720870e-12)*np.power(x,3)+(-1.9924275e-14)*np.power(x,4)+(3.6440065e-18)*np.power(x,5))

#puppi_sdm_jms_jmr = {'jms':{'value':0.999,
#                             'up'   :0.999+0.004,
#                             'down' :0.999-0.004},
#                     'jmr':{'value':1.079,
#                            'up'   :1.079+0.105,
#                            'down' :1.079-0.105}
#                  }
puppi_sdm_jms_jmr = {'jms':{'2016':1.014,
                            '2017':0.982,
                            '2018':0.982},
                     'jmr':{'2016':1.086,
                            '2017':1.092,
                            '2018':1.092}
                  }
