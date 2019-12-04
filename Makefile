
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

qml-build:
	git clone https://github.com/qmlcode/qml.git qml-build -b develop
	cd qml-build; ${PY} setup.py build --compiler=intelem --fcompiler=intelem

qml: qml-build
	ln -s qml-build/lib.linux-x86_64-3.6/qml qml

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

# Bing dataset

bing_parse_data:
	${PY} parse_bing.py --data data/bing/bp_scifinder.txt
	${PY} parse_bing.py --data data/bing/mp_scifinder.txt

bing_overview:
	${PY} plot_overview.py --dict data/bing/bp_scifinder
	${PY} plot_overview.py --dict data/bing/mp_scifinder

bing_set_structures: bing_set_structures_bp bing_set_structures_mp

bing_set_structures_bp:
	mkdir -p _tmp_bing_bp_
	${PY} prepare_structures.py --scratch _tmp_bing_bp_ --datadict data/bing/bp_scifinder -j 24

bing_set_structures_mp:
	mkdir -p _tmp_bing_mp_
	${PY} prepare_structures.py --scratch _tmp_bing_mp_ --datadict data/bing/mp_scifinder -j 24
	
bing_set_representations_bp:
	mkdir -p _tmp_bing_bp_
	touch _tmp_bing_bp_/slatm.mbtypes
	rm _tmp_bing_bp_/slatm.mbtypes
	${PY} prepare_representations.py --sdf _tmp_bing_bp_/structures.sdf.gz -j 24 --scratch _tmp_bing_bp_

bing_set_representations_mp:
	mkdir -p _tmp_bing_mp_
	touch _tmp_bing_mp_/slatm.mbtypes
	rm _tmp_bing_mp_/slatm.mbtypes
	${PY} prepare_representations.py --sdf _tmp_bing_mp_/structures.sdf.gz -j 24 --scratch _tmp_bing_mp_

bing_set_kernels_bp:
	${PY} training.py --get-kernels --scratch _tmp_bing_bp_

bing_set_kernels_mp:
	${PY} training.py --get-kernels --scratch _tmp_bing_mp_

bing_set_scores_bp:
	${PY} training.py --get-learning-curves --scratch _tmp_bing_bp_

bing_set_scores_mp:
	${PY} training.py --get-learning-curves --scratch _tmp_bing_mp_

# Bradley

bradley_parse_data:
	${PY} parse_bradley.py

bradley_overview:
	${PY} plot_overview.py --dict data/melting_bradley

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

bradley_prepare_scores:
	${PY} training.py --get-learning-curves --scratch _tmp_subset_

bradley_print_score:
	${PY} plot.py --scratch _tmp_subset_

# ochem

ochem_parse_mp:
	${PY} parse_ochem.py

ochem_overview:
	${PY} plot_overview.py --sdf 



# ALL

bradall_prepare_representations:
	${PY} prepare_representations.py --sdf data/sdf/structures.sdf.gz -j 30 --scratch _tmp_all_

bradall_prepare_kernels:
	${PY} training.py --get-kernels --scratch _tmp_all_

bradall_prepare_scores:
	${PY} training.py --get-learning-curves --scratch _tmp_all_


# Print

print_score_bradley_subset:
	${PY} plot.py --scratch _tmp_bradley_sub_

print_score_bradley_all:
	${PY} plot.py --scratch _tmp_bradley_all_

print_score_bing_bp:
	${PY} plot.py --scratch _tmp_bing_bp_

print_score_bing_mp:
	${PY} plot.py --scratch _tmp_bing_mp_


#

clean:
	rm *.pyc __pycache__ _tmp_*

super-clean:
	rm -fr data env __pycache__

