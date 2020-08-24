import pandas as pd
import re
import subprocess as sub
import matplotlib.pyplot as plt
from lib.fun_library import getZhbbWeight, getZhbbBaseCuts

x = sub.check_output('ls', shell=True).decode()
t = 'ak8'
y = '2017'
do_cuts = True
fig, ax = plt.subplots(figsize=(16,10))
#
kinem = 'NN'
kinem = 'Zh_M'

file_list = re.findall(rf'TTZH_{y}_{t}\w+_val.pkl', x)
#file_list = []

print(file_list)
for i in file_list+[f'TTZH_{y}_val.pkl']:
    df = pd.read_pickle(i)
    cuts = getZhbbBaseCuts(df)
    if do_cuts: df=df[cuts & ((df['Hbb'] == True) | (df['Zbb'] == True)) & (df['NN'] > 0.80) & (df['Zh_pt'] > 300)]
    print(i)
    print(df[kinem])
    weight = getZhbbWeight(df, y)

    ax.set_title(kinem)
    ax.hist(df[kinem], 
            #range=(0,1), 
            bins = [50,80,105,145,200],
            weights=weight, histtype='step', 
            label=f"{i.split('_')[2]}: ({sum(weight):5.2f})")

ax.legend()
#plt.yscale('log')
plt.show()
plt.close(fig)
 
