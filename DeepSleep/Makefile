# Makefile for Deepsleep

# Initialize working area my installing required pip packages
hi:	
	echo "hi"
	echo "hello"
init:
	@bash setup.sh
	@pip install -r packages/pip_requirements.txt
	@codna install --file packages/conda_requirements.txt

.PHONY: hi init clean

clean:
	@rm -rf  *~
	@rm -rf */*~
	@rm -rf *#
	@rm -rf */*#

