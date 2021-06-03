import sys
import os
import subprocess as sb
from multiprocessing import Pool
from glob import glob
import re

f_dir = '/cms/data/store/user/bcaraway/NanoAODv7/PreProcessed/2017/TTbb_EFT_2017.txt'


def main():
    with open(f_dir,'r') as _f :
        f_list = [l.replace('root://kodiak-se.baylor.edu//','/cms/data/').rstrip('\n') for l in  _f.readlines()]
        c_list = [[f_list[i],f_list[i+1]] for i in range(0,len(f_list),2)]
        pool = Pool()
        out_list = pool.map(__worker, c_list)
        os.system(f"mv {f_dir} {f_dir}old")
        with open(f_dir,'w') as o_f:
            o_f.writelines(out_list)
        

        
def __worker(files_):
    out_file = re.sub(r'(prod201\dMC)','merged\g<1>',files_[0])
    os.system("hadd {0} {1} {2}".format(out_file, files_[0],files_[1]))
    return out_file.replace('/cms/data/','root://kodiak-se.baylor.edu//')+'\n'

if __name__ == '__main__':
    main()
