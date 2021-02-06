import json
import subprocess as sb
import sys
import os
import time
import argparse
if __name__ == '__main__':
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')

import config.ana_cff as cfg
from lib.fun_library import t2Run

parser = argparse.ArgumentParser(description='Run analysis over specified sample and era')
parser.add_argument('--mode', dest='mode', type=str, 
                    choices=['sub','exe'],
                    required=True, help='mode to run', default='sub')
parser.add_argument('--inputjson', dest='inputjson', type=str, required=False, help="input json file", default=None)
parser.add_argument('--jnum', dest='job_number', type=int, required=False, help="job number", default=None)
args = parser.parse_args()
json_dir = f'{sys.path[1]}/log/nn/'

from modules.dnn_model import * # just grab the whole thing

def main():
    if args.mode == 'exe':
        run()
        exit()
    #
    # define sequences to test out
    def iterate_hpsettings():
        for d in [.5]:
            sequences = [
                [['Dense', 128],['Dense', 64],['Dropout',d]],
                #[['Dense', 256],['Dense', 128],['Dropout',d]],
            ]
            for seq in sequences:
                for n_epoch in [100]:
                    for batch in [5128*2]:
                        for lr in [0.0003]:
                            for g in [.25,.5,.75,1]:
                                for tt_a in [.25,0.5,.75,1,1.25,1.5,2][::-1]:
                                    for ttbb_a in [.25,.5,0.75,1,1.25,1.5,2][::-1]:
                                        for ttzh_a in [.25,.5,.75,1,1.25,1.5,2]:
                                            yield {
                                                'sequence':seq,
                                                'other_settings':{
                                                    'fl_a':[tt_a,ttbb_a,ttzh_a],
                                                    'fl_g':g,
                                                    'lr_alpha':lr
                                                },
                                                'n_epochs':n_epoch,
                                                'batch_size':batch,
                                            }
    #
    for i,m_info in enumerate(iterate_hpsettings()):
        json_name = f'testdnn_{i}.json'
        with open(json_dir+json_name,'w') as jf:
            json.dump(m_info,jf)
        command = f"qsub -l nodes=1:ppn=1 -N testDNN_{i} -o {json_dir}/job_{i}.stdout -e {json_dir}/job_{i}.stderr "
        command += f" -v json={json_name},jobnum={i} {sys.path[1]}/scripts/testDNN.sh"
        num_jobs_running = lambda: int(sb.check_output('qstat -u $USER | grep testDNN | wc -l', shell=True).decode())
        #num_jobs_running = lambda: int(sb.check_output('jobs | wc -l', shell=True).decode())
        while num_jobs_running() >= 35:
            time.sleep(30)
        print(command)
        os.system(command)
                                        
@t2Run
def run():
    from sklearn.metrics import confusion_matrix
    from sklearn.metrics import roc_auc_score
    import matplotlib.pyplot as plt
    print(args.inputjson, args.job_number)
    m_info = json.load(open(json_dir+args.inputjson,'r'))
    print(m_info)
    #
    model, testX, testY = train_model(m_info)
    y_pred = model.predict(testX)
    #
    weight = np.ones_like(y_pred[:,2])*.001
    weight = np.where(testY[:,0]==1,.10,weight)
    weight = np.where(testY[:,1]==1,.15,weight)
    sig_roc_auc = roc_auc_score(testY[:,2],y_pred[:,2],sample_weight=weight)
    print('AUC score', sig_roc_auc)
    if sig_roc_auc < .8:
        exit()
    #
    loss ,acc, auc = model.evaluate(testX, testY)#, sample_weight = testW['DNNweight'].values)
    cm = confusion_matrix(np.argmax(testY,axis=1), np.argmax(y_pred,axis=1))
    print ("\nThe confusion matrix of the test set on the trained nerual network:\n" , cm)
    # contruct score
    #
    zh_score = cm[:,2]
    print(cm[:,2])
    s_b_val = zh_score[2]*.001/((zh_score[0]*0.10+zh_score[1]*.15)**(1/2))
    s_ttbb_val = zh_score[2]*.001/((zh_score[1]*.15)**(1/2))
    cm = confusion_matrix(np.argmax(testY[y_pred[:,2]>0.6],axis=1), np.argmax(y_pred[y_pred[:,2]>0.6],axis=1)) 
    print ("\nThe confusion matrix of the test set (high confidence in zh):\n" , cm) 
    
    out_name = f"{sig_roc_auc:.3f}_{m_info['other_settings']['fl_a']}_{m_info['other_settings']['fl_g']}_{s_ttbb_val:.2f}_{s_b_val:.2f}_{args.job_number}"
    print(out_name)
    model.save_weights(cfg.dnn_ZH_dir+'/test_archs/'+out_name+'.h5')

    plt.hist(y_pred[testY[:,2] == 1][:,2],
             bins=20, range=(0,1), histtype='step', 
             weights=np.ones(len(y_pred[testY[:,2] == 1][:,2]))*.001, label='SIG')
    plt.hist(y_pred[testY[:,0] == 1][:,2],
             bins=20, range=(0,1), histtype='step', 
             weights=np.ones(len(y_pred[testY[:,0] == 1][:,2]))*.10, label='tt')
    plt.hist(y_pred[testY[:,1] == 1][:,2],
             bins=20, range=(0,1), histtype='step', 
             weights=np.ones(len(y_pred[testY[:,1] == 1][:,2]))*.15, label='ttbb')
    plt.legend()
    plt.yscale('log')
    plt.xlim(0,1)
    plt.title(out_name)
    #plt.show()
    plt.savefig(cfg.dnn_ZH_dir+'/test_archs/'+out_name+'.pdf')

if __name__ == '__main__':
    main()

