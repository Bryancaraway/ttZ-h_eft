# Makefile for Deepsleep

# Initialize working area my installing required pip packages
hi:	
	echo "hi"
	echo "hello"
init:
	@bash setup.sh
	@pip install -r packages/requirements.txt
	@conda install numpy
	@conda install pandas

.PHONY: hi init clean

clean:
	@rm -rf  *~
	@rm -rf */*~
	@rm -rf *#
	@rm -rf */*#

