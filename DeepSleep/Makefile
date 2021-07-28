# Makefile for Deepsleep

# Initialize working area my installing required pip packages
hi:	
	echo "hi"
	echo "hello"
init:
	@bash setup.sh
	@pip install -r packages/requirements.txt
	@conda env create -f packages/ttxenv2.0.yml
#@conda install --yes numpy
#@conda install --yes pandas
#@conda install --yes -c conda-forge xrootd
#@conda install --yes tensorflow

.PHONY: hi init clean

clean:
	@rm -rf  *~
	@rm -rf */*~
	@rm -rf *#
	@rm -rf */*#

