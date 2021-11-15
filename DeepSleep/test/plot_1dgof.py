import sys
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output(
        'echo $(git rev-parse --show-cdup)', 
        shell=True).decode().strip('\n')+'DeepSleep/')
#
from matplotlib import rc
rc("savefig",dpi=250)
rc("figure", max_open_warning=600)
rc("figure", figsize=(8, 6*(6./8.)), dpi=200)

from lib.fun_library import save_pdf
from post_fit import PostFit

def main():
    #f_list = ['fitDiagnostics_unblind_run2.root']
    f_list = ['fitDiagnostics_unconstrained_run2.root','fitDiagnostics_constrained_run2.root','fitDiagnostics_unblind_run2.root']
    for f in f_list:
        goffit = GOFonly(f)
        save_pdf(f.replace('root','pdf'))(goffit.makeplots)()

class GOFonly(PostFit):
    
    def __init__(self,fitroo):
        super().__init__(fitroo)

    def load_info(self):
        pass # hack

    def makeplots(self, doPull=True, do1Pull=True):
        self.plot_1dgof()
    

if __name__ == '__main__':
    main()
