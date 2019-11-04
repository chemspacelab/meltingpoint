
#PY=./env/bin/python
PY=python
CONDA=conda

all: env pip chemhelp

setup: env pip chemhelp

env:
	${CONDA} env create -f requirements.yml -p env

pip: env
	${PY} -m pip install numpy
	${PY} -m pip install -r requirements.txt --no-cache-dir

chemhelp:
	git clone https://github.com/charnley/chemhelp

statsig:
	git clone https://github.com/jensengroup/statsig

data:
	mkdir -p data

#

structures:
	${PY} prepare_structures.py

overview:
	${PY} prepare_subset.py

representations:
	${PY} prepare_representations.py --sdf data/sdf/subset_structures.sdf

kernels:
	${PY} training.py

score:
	${PY} training.py --plot --scratch _tmp_subset_

# Bradley

bradley_parse_data:
	${PY} parse_bradley.py

bradley_prepare_structures:
	${PY} prepare_structures.py --data data/melting_bradley_clean -j 30

bradley_prepare_conformers:
	${PY} prepare_structures.py --sdf data/sdf/structures.sdf.gz -j 30

bradley_prepare_subset:
	${PY} prepare_subset.py --sdf data/sdf/structures.sdf.gz -j 30

bradley_prepare_representations:
	${PY} prepare_representations.py --sdf data/sdf/subset_structures.sdf -j 30 --scratch _tmp_subset_

bradley_prepare_representations_conformers:
	${PY} prepare_representations.py --conformers -j 30 --scratch _tmp_subset_

bradley_prepare_kernels:
	${PY} training.py --get-kernels --scratch _tmp_subset_

bradley_prepare_curves:
	${PY} training.py --get-learning-curves --scratch _tmp_subset_

bradley_print_score:
	${PY} plot.py --scratch _tmp_subset_

#

clean:
	rm *.pyc __pycache__ _tmp_*

super-clean:
	rm -fr data env __pycache__

