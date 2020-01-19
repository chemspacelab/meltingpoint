
PY=python
CONDA=conda
BIN=src

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
	ln -s qml-build/lib.linux-x86_64-3.6/qml ${BIN}/qml

FFLAGS=-xHost -O4 -qopenmp
fbitmap:
	cd ${BIN}; f2py -c --verbose --opt="${FFLAGS}" --compiler=intelem --fcompiler=intelem -m bitmap_kernels fbitmap.f90

chemhelp:
	git clone https://github.com/charnley/chemhelp ${BIN}/chemhelp

statsig:
	git clone https://github.com/jensengroup/statsig ${BIN}statsig

data:
	mkdir -p data


## DATASETS

# Bing dataset

BINGMP=_tmp_bing_mp_
BINGBP=_tmp_bing_bp_

bing_parse_data:
	${PY} ${BIN}/parse_bing.py --data data/bing/bp_scifinder.txt --scratch ${BINGBP}
	${PY} ${BIN}/parse_bing.py --data data/bing/mp_scifinder.txt --scratch ${BINGMP}

bing_overview:
	${PY} ${BIN}/plot_overview.py --dict data/bing/bp_scifinder
	${PY} ${BIN}/plot_overview.py --dict data/bing/mp_scifinder

bing_set_structures: bing_set_structures_bp bing_set_structures_mp

bing_set_structures_bp:
	mkdir -p _tmp_bing_bp_
	${PY} ${BIN}/prepare_structures.py --scratch _tmp_bing_bp_ --datadict data/bing/bp_scifinder -j 24

bing_set_structures_mp:
	mkdir -p _tmp_bing_mp_
	${PY} ${BIN}/prepare_structures.py --scratch _tmp_bing_mp_ --datadict data/bing/mp_scifinder -j 24
	
bing_set_representations_bp:
	mkdir -p _tmp_bing_bp_
	touch _tmp_bing_bp_/slatm.mbtypes
	rm _tmp_bing_bp_/slatm.mbtypes
	${PY} ${BIN}/prepare_representations.py --sdf _tmp_bing_bp_/structures.sdf.gz -j 24 --scratch _tmp_bing_bp_

bing_set_representations_mp:
	mkdir -p _tmp_bing_mp_
	touch _tmp_bing_mp_/slatm.mbtypes
	rm _tmp_bing_mp_/slatm.mbtypes
	${PY} ${BIN}/prepare_representations.py --sdf _tmp_bing_mp_/structures.sdf.gz -j 24 --scratch _tmp_bing_mp_

bing_set_kernels_bp:
	${PY} ${BIN}/training.py --get-kernels --scratch _tmp_bing_bp_

bing_set_kernels_mp:
	${PY} ${BIN}/training.py --get-kernels --scratch _tmp_bing_mp_

bing_set_scores_bp:
	${PY} ${BIN}/training.py --get-learning-curves --scratch _tmp_bing_bp_

bing_set_scores_mp:
	${PY} ${BIN}/training.py --get-learning-curves --scratch _tmp_bing_mp_

# Bradley

BRADLEYMP=_tmp_bradley_all_

bradley_parse_data:
	${PY} ${BIN}/parse_bradley.py --scratch ${BRAD}

bradley_overview:
	${PY} ${BIN}/plot_overview.py --dict data/melting_bradley

bradley_prepare_structures:
	${PY} ${BIN}/prepare_structures.py --data data/melting_bradley_clean -j 30

bradley_prepare_conformers:
	${PY} ${BIN}/prepare_structures.py --sdf data/sdf/structures.sdf.gz -j 30

bradley_prepare_subset:
	${PY} ${BIN}/prepare_subset.py --sdf data/sdf/structures.sdf.gz -j 30

bradley_prepare_representations:
	${PY} ${BIN}/prepare_representations.py --sdf data/sdf/subset_structures.sdf -j 30 --scratch _tmp_subset_

bradley_prepare_representations_conformers:
	${PY} ${BIN}/prepare_representations.py --conformers -j 30 --scratch _tmp_subset_

bradley_prepare_kernels:
	${PY} ${BIN}/training.py --get-kernels --scratch _tmp_subset_

bradley_prepare_scores:
	${PY} ${BIN}/training.py --get-learning-curves --scratch _tmp_subset_


# ALL

bradall_prepare_representations:
	${PY} ${BIN}/prepare_representations.py --sdf data/sdf/structures.sdf.gz -j 30 --scratch _tmp_all_

bradall_prepare_kernels:
	${PY} ${BIN}/training.py --get-kernels --scratch _tmp_all_

bradall_prepare_scores:
	${PY} ${BIN}/training.py --get-learning-curves --scratch _tmp_all_



# OCHEM

OCHEMDAT=_tmp_ochem_data_
OCHEMBP=_tmp_ochem_bp_
OCHEMMP=_tmp_ochem_mp_

ochem_mp_parse:
	mkdir -p ${OCHEMMP}
	${PY} ${BIN}/parse_ochem.py --scratch ${OCHEMMP} -j 24 --sdf \
	${OCHEMDAT}/meltingpoints_0_100.sdf.gz \
	${OCHEMDAT}/meltingpoints_100_200.sdf.gz \
	${OCHEMDAT}/meltingpoints_200_250.sdf.gz \
	${OCHEMDAT}/meltingpoints_250_300.sdf.gz \
	${OCHEMDAT}/meltingpoints_300_350.sdf.gz \
	${OCHEMDAT}/meltingpoints_350_450.sdf.gz \
	${OCHEMDAT}/meltingpoints_450_x.sdf.gz

ochem_bp_parse:
	mkdir -p ${OCHEMBP}
	${PY} ${BIN}/parse_ochem.py --scratch ${OCHEMBP} -j -1 --sdf \
		${OCHEMDAT}/boilingpoints_all.sdf.gz

ochem_overview:
	${PY} ${BIN}/plot_overview.py --dict ${OCHEMBP}/molecule_data
	${PY} ${BIN}/plot_overview.py --dict ${OCHEMMP}/molecule_data

ochem_bp_set_xyz:
	${PY} ${BIN}/prepare_structures.py -j 24 \
		--datadict ${OCHEMBP}/molecule_data \
		--scratch ${OCHEMBP}

ochem_bp_set_rep:
	touch ${OCHEMBP}/slatm.mbtypes
	rm ${OCHEMBP}/slatm.mbtypes
	${PY} ${BIN}/prepare_representations.py -j -1 \
		--sdf ${OCHEMBP}/structures.sdf.gz \
		--scratch ${OCHEMBP} \
		--representations "rdkitfp" "morgan"

ochem_bp_set_kernel:
	${PY} ${BIN}/prepare_kernels.py \
		-j 24 \
		--scratch ${OCHEMBP} \
		--representations "rdkitfp" "morgan"

ochem_bp_set_score:
	${PY} ${BIN}/training.py --get-learning-curves --scratch ${OCHEMBP}


## pubchem

PUBCHEMDATA=_tmp_pubchem_data_
PUBCHEMMP=_tmp_pubchem_mp_
PUBCHEMBP=_tmp_pubchem_bp_

pubchem_mp_parse:
	mkdir -p ${PUBCHEMMP}
	${PY} ${BIN}/parse_pubchem.py --scratch ${PUBCHEMMP} -j 24 --json ${PUBCHEMDATA}/meltingpoints.json

pubchem_bp_parse:
	mkdir -p ${PUBCHEMBP}
	${PY} ${BIN}/parse_pubchem.py --scratch ${PUBCHEMBP} -j 24 --json ${PUBCHEMDATA}/boilingpoints.json

pubchem_overview:
	${PY} ${BIN}/plot_overview.py --dict ${PUBCHEMBP}/molecule_data
	${PY} ${BIN}/plot_overview.py --dict ${PUBCHEMMP}/molecule_data

## MERGE

MERGEBP=_tmp_merge_bp_
MERGEMP=_tmp_merge_mp_
MERGEMPFIL=_tmp_mergefilter_mp_

merge_bp:
	mkdir -p ${MERGEBP}
	${PY} ${BIN}/merge.py --dict \
	${BINGBP}/molecule_data \
	${OCHEMBP}/molecule_data \
	${PUBCHEMBP}/molecule_data \
	--scratch ${MERGEBP}

merge_mp:
	mkdir -p ${MERGEMP}
	${PY} ${BIN}/merge.py --dict \
	${BINGMP}/molecule_data \
	${OCHEMMP}/molecule_data \
	${BRADLEYMP}/molecule_data \
	${PUBCHEMMP}/molecule_data \
	--scratch ${MERGEMP}

mergefiltermp:
	mkdir -p ${MERGEMPFIL}
	${PY} ${BIN}/merge.py --dict \
	${BINGMP}/molecule_data \
	${OCHEMMP}/molecule_data \
	${BRADLEYMP}/molecule_data \
	${PUBCHEMMP}/molecule_data \
	--scratch ${MERGEMPFIL} \
	--filter

merge_overview:
	${PY} ${BIN}/plot_overview.py --dict ${MERGEBP}/molecule_data
	# ${PY} ${BIN}/plot_overview.py --dict ${MERGEMP}/molecule_data

merge_bp_set_xyz:
	${PY} ${BIN}/prepare_structures.py -j 24 \
		--datadict ${MERGEBP}/molecule_data \
		--scratch ${MERGEBP}

merge_mp_set_xyz:
	${PY} ${BIN}/prepare_structures.py -j 24 \
		--datadict ${MERGEMP}/molecule_data \
		--scratch ${MERGEMP}

mergefiltermp_set_xyz:
	${PY} ${BIN}/prepare_structures.py -j 24 \
		--datadict ${MERGEMPFIL}/molecule_data \
		--scratch ${MERGEMPFIL}

merge_bp_set_rep:
	touch ${MERGEBP}/slatm.mbtypes
	rm ${MERGEBP}/slatm.mbtypes
	time ${PY} ${BIN}/prepare_representations.py -j 24 \
		--sdf ${MERGEBP}/structures.sdf.gz \
		--scratch ${MERGEBP} \
		--representations "rdkitfp" "morgan"

merge_mp_set_rep:
	touch ${MERGEMP}/slatm.mbtypes
	rm ${MERGEMP}/slatm.mbtypes
	time ${PY} ${BIN}/prepare_representations.py -j -1 \
		--sdf ${MERGEMP}/structures.sdf.gz \
		--scratch ${MERGEMP} \
		--representations "rdkitfp" "morgan"

mergefiltermp_set_rep:
	touch ${MERGEMPFIL}/slatm.mbtypes
	rm ${MERGEMPFIL}/slatm.mbtypes
	time ${PY} ${BIN}/prepare_representations.py -j 24 \
		--sdf ${MERGEMPFIL}/structures.sdf.gz \
		--scratch ${MERGEMPFIL} \
		--representations "fchl19" #"rdkitfp"



merge_bp_set_kernel:
	time ${PY} ${BIN}/prepare_kernels.py \
		-j 40 \
		--scratch ${MERGEBP} \
		--representations "rdkitfp" "morgan"

merge_mp_set_kernel:
	time ${PY} ${BIN}/prepare_kernels.py \
		-j 40 \
		--scratch ${MERGEMP} \
		--representations "rdkitfp"

mergefiltermp_set_kernel:
	time ${PY} ${BIN}/prepare_kernels.py \
		-j -1 \
		--scratch ${MERGEMPFIL} \
		--representations "rdkitfp"

merge_bp_set_scores:
	${PY} ${BIN}/training.py --get-learning-curves --scratch ${MERGEBP}

merge_mp_set_scores:
	${PY} ${BIN}/training.py --get-learning-curves --scratch ${MERGEMP}

mergefiltermp_set_scores:
	${PY} ${BIN}/training.py --get-learning-curves --scratch ${MERGEMPFIL} --names "rdkitfp"

## Small subset of druglike molecule
# 1-2 aromatic rings
FILTERMP=_tmp_filtermp_
subset_mp_set_xyz:
	mkdir -p ${FILTERMP}
	${PY} ${BIN}/parse_filter.py --sdf ${MERGEMP}/structures.sdf.gz --properties ${MERGEMP}/properties.csv --scratch ${FILTERMP}

subset_mp_set_rep:
	touch ${FILTERMP}/slatm.mbtypes
	rm ${FILTERMP}/slatm.mbtypes
	time ${PY} ${BIN}/prepare_representations.py -j 24 \
		--sdf ${FILTERMP}/structures.sdf.gz \
		--scratch ${FILTERMP} \
		--representations "slatm" "rdkitfp"
	@# time ${PY} ${BIN}/prepare_representations.py -j 24 \
	@# 	--sdf ${FILTERMP}/structures.sdf.gz \
	@# 	--scratch ${FILTERMP} \
	@# 	--representations "rdkitfp"

subset_mp_set_kernel:
	time ${PY} ${BIN}/prepare_kernels.py \
		-j -1 \
		--scratch ${FILTERMP} \
		--representations "slatm" "rdkitfp"

subset_mp_set_scores:
	touch ${FILTERMP}/properties.npy
	rm ${FILTERMP}/properties.npy
	${PY} ${BIN}/training.py --get-learning-curves --scratch ${FILTERMP} --names "rdkitfp"
	${PY} ${BIN}/training.py --get-learning-curves --scratch ${FILTERMP} --names "slatm"

subset_mp_ols:
	rm ${FILTERMP}/repr.ols.npy
	${PY} ${BIN}/training_ols.py --scratch ${FILTERMP} -j24

subset_overview:
	${PY} ${BIN}/plot_overview.py --dict ${FILTERMP}/molecules

subset_view_kernel:
	${PY} ${BIN}/plot_kernel.py --dist ${FILTERMP}/dist.slatm.npy --scratch ${FILTERMP}
	${PY} ${BIN}/plot_kernel.py --kernel ${FILTERMP}/kernel.rdkitfp.npy --scratch ${FILTERMP}

subset_mp_null:
	${PY} ${BIN}/training_null.py --scratch ${FILTERMP} -j24

## PRINT RESULTS

print_score_bradley_subset:
	${PY} ${BIN}/plot.py --scratch _tmp_bradley_sub_

print_score_bradley_all:
	${PY} ${BIN}/plot.py --scratch ${BRAD}

print_score_bing_bp:
	${PY} ${BIN}/plot.py --scratch ${BINGBP}

print_score_bing_mp:
	${PY} ${BIN}/plot.py --scratch ${BINGMP}

print_score_ochem_bp:
	${PY} ${BIN}/plot.py --scratch ${OCHEMBP}

print_score_ochem_mp:
	${PY} ${BIN}/plot.py --scratch ${OCHEMMP}

print_score_merge_mp:
	${PY} ${BIN}/plot.py --scratch ${MERGEMP}

print_score_subset_mp:
	${PY} ${BIN}/plot.py --scratch ${FILTERMP}

## MISC

cleancache:
	rm -r __pycache__

clean:
	rm -r *.pyc __pycache__ .pycache

super-clean:
	rm -fr data env __pycache__

