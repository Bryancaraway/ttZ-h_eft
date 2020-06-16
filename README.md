# Setup

## Setup Intructions for DeepSleep on Kodiak ONLY
 
### Setup Virtual Environment

It is necessary to setup a virtual environment to manage the analysis' python package requirements
One way of accomplishing this is through Anaconda3

For Anaconda3 installation (if not already present):

```
mkdir ~/tmp
cd ~/tmp
curl -O https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh
bash Anaconda3-2019.10-Linux-x86_64.sh
```

In .bashrc or .bash_profile (etc...) there should be an appened section after install:

```
__conda_setup="$('/home/$USER/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/$USER/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/$USER/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/$USER/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
```

To start the process of creating a virtual env through Anaconda3 (if not done so prior):

```
source ~/.bashrc
```

Create a new conda env with python pre-installed python3:


```
conda create --name ttxenv python=3.7
conda activate ttxenv
```

Note: it is necessary to execute 'conda activate ttxenv' at the start of logging in
To deactivate the virtual environment, simply:

```
conda deactivate
```

## Download Repository

Make a new directory for ttZ-h analysis workspace:

```
mkdir ttZh_analysis
cd ttZh_analysis
```
Then clone and checkout a new working branch from github

```
git clone https://github.com/Bryancaraway/ttZ-h_eft.git

git checkout -b working_branch
cd DeepSleep
```

## Setup file structure and python packages

After installing the repository, it is required to setup the necessary file structures and dependencies
for the analysis framework:

```
make init
```

Please be patient as this might take a while if you do not have python packages already cached with Anaconda3

## Run the analysis and generate .pkl files 

With everything setup, the analysis may run over skim ntuple root files (currently, only 2017 for MC and Data available)
These and other analysis relevent files can be found here although it won't be necessary to do anything with them:
(In the future, I may place my personal .pkl output files here.)

```
ls -l /cms/data/store/user/ttxeft/
```
To run through the analysis framework and populate the files directory:
Note: This will submit jobs using Kodiak's pbs batch job system. 
(This process will take 15-20 minutes to complete.) 

```
python submit_jobs.py
```

To check on the status of your jobs:
(.out and .err files will begin to appear log/ after jobs complete)

```
qstat -u $USER
```

Alternatively, the option exists to run over one MC sample or data sample at a time:

```
python runAna.py TTZH
```

These python scripts rely on custom python modules and a configureation file:
```
ls -l cfg/
```

