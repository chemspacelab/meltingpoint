#!/bin/bash
#
# Note:
# Download first page to find page ranges, then do for-loop and download all
# pages for melting and boiling phases
#

phase=$1
page=$2

echo $phase $page

wget "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/annotations/heading/$phase%20Point/JSON?heading_type=Compound&page=$page" -O data/${phase}_$page.json 2> /dev/null

