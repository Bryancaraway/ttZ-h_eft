# run full analysis on kodiak's batch system
# use case: python submitjobs.py
#
import sys
import os
import argparse
import config.ana_cff as cfg
#
# add parsargs at some point for year, rundata, minibatch
parser = argparse.ArgumentParser(description='Run analysis over many samples and years using batch system')
parser.add_argument('-s', dest='samples', type=str, 
                    choices=cfg.All_MC+cfg.Data_samples+cfg.tt_sys_samples+['all','mc','data','tt','ttsys','ttbb','ttbbsys']+cfg.Sig_EFT_MC, 
                    required=True, help='analyze all, mc only, or data only')
parser.add_argument('-y', dest='year', type=str, choices=cfg.Years+['all'],
                    required=True, help='year or all years')
parser.add_argument('-j', dest='jec',     type=str, required=False, help='Run with specified jec variation', choices=['JESUp','JESDown','JERUp','JERDown','all'], default=None)
parser.add_argument('--jjec', dest='jetjec', type=str, required=False, help='Which jet to compute jec', choices=['ak4','ak8','all'], default=None)
#parser.add_argument('-t', dest='tag',     type=str, required=False, help='Optional tag to add to output file', default='')
args = parser.parse_args()

sample_dict = {'all':cfg.All_MC+cfg.Data_samples,
               'mc':cfg.All_MC,
               'data':cfg.Data_samples,
               'tt'   : cfg.ttbar_samples,
               'ttsys': cfg.tt_sys_samples,
               'ttbb' : cfg.tt_bb,
               'ttbbsys': cfg.tt_bb_sys}

# submit a job for each background / signal / data sample
job_script = 'scripts/runAna.sh' if args.jec is None else 'scripts/runAna_jec.sh'
samples = sample_dict.get(args.samples,[args.samples])
years = [args.year] if args.year != 'all' else cfg.Years
#jecs  = [args.jec]  if args.jec != 'all' else ['JESUp','JESDown','JERUp','JERDown']
jecs    = ['JESUp','JESDown','JERUp','JERDown'] if args.jec == 'all' else [args.jec]
jetjecs = ['ak4','ak8'] if args.jetjec == 'all' else [args.jetjec]

for year in years:
    log_dir = f'log/{year}/'
    os.system(f'rm {log_dir}*{year}*')
    for sample in samples:
        for jetjec in jetjecs:
            for jec in jecs:
                add_args  = ''
                add_out_name = ''
                if jec is not None:
                    add_args = f',jec={jec},jjec={jetjec}'
                    add_out_name = '_'+jec+'_'+jetjec
                #
                out_name  = sample+'_'+year+add_out_name
                if sample == 'TTBarSemi_pow' and year == '2018':
                    ppn = 4
                else:
                    ppn = 1
                pass_args = f'-v sample={sample},year={year}{add_args}'
                command   = f'qsub -l nodes=1:ppn={ppn} -o {log_dir}{out_name}.out -e {log_dir}{out_name}.err {pass_args} {job_script}'
                print(command)
                os.system(command)

