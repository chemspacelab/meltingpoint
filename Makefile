
PY=python
CONDA=conda

all:
	@echo "Read the Makefile"


## SETUP

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


## DATASETS

# Bing dataset

BINGMP=_tmp_bing_mp_
BINGBP=_tmp_bing_mp_

bing_parse_data:
	${PY} parse_bing.py --data data/bing/bp_scifinder.txt --scratch ${BINGBP}
	${PY} parse_bing.py --data data/bing/mp_scifinder.txt --scratch ${BINGMP}

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

BRAD=_tmp_bradley_all_

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


# ALL

bradall_prepare_representations:
	${PY} prepare_representations.py --sdf data/sdf/structures.sdf.gz -j 30 --scratch _tmp_all_

bradall_prepare_kernels:
	${PY} training.py --get-kernels --scratch _tmp_all_

bradall_prepare_scores:
	${PY} training.py --get-learning-curves --scratch _tmp_all_



# OCHEM

OCHEMDAT=_tmp_ochem_data_
OCHEMBP=_tmp_ochem_bp_
OCHEMMP=_tmp_ochem_mp_

ochem_mp_parse:
	mkdir -p ${OCHEMMP}
	${PY} parse_ochem.py --scratch ${OCHEMMP} -j 24 --sdf \
	${OCHEMDAT}/meltingpoints_0_100.sdf.gz \
	${OCHEMDAT}/meltingpoints_100_200.sdf.gz \
	${OCHEMDAT}/meltingpoints_200_250.sdf.gz \
	${OCHEMDAT}/meltingpoints_250_300.sdf.gz \
	${OCHEMDAT}/meltingpoints_300_350.sdf.gz \
	${OCHEMDAT}/meltingpoints_350_450.sdf.gz \
	${OCHEMDAT}/meltingpoints_450_x.sdf.gz

ochem_bp_parse:
	mkdir -p ${OCHEMBP}
	${PY} parse_ochem.py --scratch ${OCHEMBP} -j 24 --sdf \
		${OCHEMDAT}/boilingpoints_all.sdf.gz

ochem_overview:
	# ${PY} plot_overview.py --dict ${OCHEMMP}/molecule_data
	${PY} plot_overview.py --dict ${OCHEMBP}/molecule_data

ochem_bp_set_xyz:
	${PY} prepare_structures.py -j 24 \
		--datadict ${OCHEMBP}/molecule_data \
		--scratch ${OCHEMBP}

ochem_bp_set_rep:
	touch ${OCHEMBP}/slatm.mbtypes
	rm ${OCHEMBP}/slatm.mbtypes
	${PY} prepare_representations.py -j 24 \
		--sdf ${OCHEMBP}/structures.sdf.gz \
		--scratch ${OCHEMBP}

ochem_bp_set_kernel:
	${PY} training.py --get-kernels --scratch ${OCHEMBP}

ochem_bp_set_score:
	${PY} training.py --get-learning-curves --scratch ${OCHEMBP}



## PRINT RESULTS

print_score_bradley_subset:
	${PY} plot.py --scratch _tmp_bradley_sub_

print_score_bradley_all:
	${PY} plot.py --scratch ${BRAD}

print_score_bing_bp:
	${PY} plot.py --scratch ${BINGBP}

print_score_bing_mp:
	${PY} plot.py --scratch ${BINGMP}

print_score_ochem_bp:
	${PY} plot.py --scratch ${OCHEMBP}

print_score_ochem_mp:
	${PY} plot.py --scratch ${OCHEMMP}


## MERGE

merge_bp:
	${PY} merge.py --sdf \
	_tmp_ochem_bp_/structures.sdf.gz \
	_tmp_bing_bp_/structures.sdf.gz \
	--name OCHEM SCIFINDER \
	--filename _fig_overlap_bp

merge_mp:
	${PY} merge.py --sdf \
	_tmp_bradley_all_/structures.sdf.gz \
	_tmp_bing_mp_/structures.sdf.gz \
	--dict \
	_tmp_ochem_mp_/molecule_data \
	--name BRADLEY SCIFINDER OCHEM \
	--filename _fig_overlap_mp


## MISC

cleancache:
	rm -r __pycache__

clean:
	rm -r *.pyc __pycache__ .pycache

super-clean:
	rm -fr data env __pycache__

