#Setup

Setup Intructions for DeepSleep on Kodiak ONLY
 
=======================================

You will need some python packages available through Anaconda3

```
mkdir tmp
cd /tmp
curl -O https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh
```

Your bashrc should have the appended section after install:

```
__conda_setup="$('/home/bcaraway/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/bcaraway/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/bcaraway/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/bcaraway/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
```

Now, it will be necessary to source bashrc

```
source ~/.bashrc
```

Then create a new conda env with python pre-installed python3

```
conda create --name deepsleep python=3
conda activate deepsleep
```

To add python packages to anaconda, simply:

```
conda install [package-name-here]
```

or alternatively

```
pip install package-name-here)
```

My installed packages can be found in packages/condaList.txt and packages/pipList.txt

=======================================

Clone to your working area and move root files to your area

```
git clone https://github.com/Bryancaraway/ZInvisible.git
cd ZInvisible
git checkout TTX
cd Tools/DeepSleep
bash setupMl.sh
```

=======================================

At the start of every session:

```
conda activate mlenv
```
If on kodiak gpu node:
```
module load cuda91
```

======================================

To process for training

```
python processData.py
```
Then train

```
python train.py
```
then validate the training

```
python validation.py
```

To process for kinematic fit run

```
python kinematicFit.py
```

=====================================