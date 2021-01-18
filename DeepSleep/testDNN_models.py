import json
if __name__ == '__main__':
    import subprocess as sb
    sys.path.insert(1, sb.check_output('echo $(git rev-parse --show-cdup)', shell=True).decode().strip('\n')+'DeepSleep/')
from modules.dnn_model import *

parser = argparse.ArgumentParser(description='Run analysis over specified sample and era')
parser.add_argument('--mode', dest='mode', type=str, 
                    choices=['sub','exe']
                    required=True, help='mode to run')
parser.add_argument('--inputjson', dest='inputjson', type=str, required=False, help="input json file", default=None)
parser.add_argument('--jnum', dest='job_number', type=int, required=False, help="job number", default=None)
args = parser.parse_args()
json_dir = f'sys.path[1]/log/nn/'


def main():
    if args.mode == 'exe':
        run()
    #
    command = f"qsub -l nodes=1:ppn=1 -N testDNN -o {json_dir}/name.stdout -e {json_dir}/name.stderr"
    command += "-v json=json,jobnum=jobnum testDNN.sh"
    m_info = {
        'sequence':[['Dense', 128],
                    ['Dense', 64],
                    ['Dropout',0.5]],
        'other_settings': {'fl_a':[0.6,2.0,1.0],
                           'fl_g':2,
                           'lr_alpha':0.0001,},
        'n_epochs':5,
        'batch_size':5128,
    }

def run():
    m_info = json.load(open(json_dir+args.inputjson,'r'))
    #
    model, testX, testY = train_model(m_info)
    from sklearn.metrics import confusion_matrix
    import matplotlib.pyplot as plt
    y_pred = model.predict(testX)
    cm = confusion_matrix(np.argmax(testY,axis=1), np.argmax(y_pred,axis=1))
    print ("\nThe confusion matrix of the test set on the trained nerual network:\n" , cm)
    # contruct score
    zh_score = cm[:,2]
    print(cm[:,2])
    s_b_val = zh_score[2]*.003/((zh_score[0]*0.14+zh_score[1]*.015)**(1/2))
    s_ttbb_val = zh_score[2]*.003/((zh_score[1]*.015)**(1/2))
    cm = confusion_matrix(np.argmax(testY[y_pred[:,2]>0.6],axis=1), np.argmax(y_pred[y_pred[:,2]>0.6],axis=1)) 
    print ("\nThe confusion matrix of the test set (high confidence in zh):\n" , cm) 
    out_name = f'{s_b_val:.2f}_{s_ttbb_val:.2f}_{args.job_number}'
    print(out_name)
    model.save_weights(cfg.dnn_ZH_dir+'/test_archs/'+out_name+'.h5')

    plt.hist(y_pred[testY[:,2] == 1][:,2],
             bins=10, range=(0,1), histtype='step', 
             weights=np.ones(len(y_pred[testY[:,2] == 1][:,2]))*.003, label='SIG')
    plt.hist(y_pred[testY[:,0] == 1][:,2],
             bins=10, range=(0,1), histtype='step', 
             weights=np.ones(len(y_pred[testY[:,0] == 1][:,2]))*.14, label='tt')
    plt.hist(y_pred[testY[:,1] == 1][:,2],
             bins=10, range=(0,1), histtype='step', 
             weights=np.ones(len(y_pred[testY[:,1] == 1][:,2]))*.015, label='ttbb')
    plt.legend()
    plt.yscale('log')
    plt.xlim(0,1)
    plt.title(out_name)
    plt.show()
    #plt.save_fig(cfg.dnn_ZH_dir+'/test_archs/'+out_name+.pdf)

if __name__ == '__main__':
    main()

