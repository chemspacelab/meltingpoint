
PYTHON=./env/bin/python
CONDA=conda

all: env pip chemhelp

setup: env pip chemhelp

env:
	${CONDA} env create -f requirements.yml -p env

pip: env
	${PYTHON} -m pip install numpy
	${PYTHON} -m pip install -r requirements.txt --no-cache-dir

chemhelp:
	git clone https://github.com/charnley/chemhelp

statsig:
	git clone https://github.com/jensengroup/statsig

data:
	mkdir -p data

#

train:
	env/bin/python training.py

#

clean:
	rm *.pyc __pycache__ _tmp_*

super-clean:
	rm -fr data env __pycache__

