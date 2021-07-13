#!/bin/bash

CV=$1
ION=$2

INPUT_LIST="/home/bach/Documents/doctoral/projects/massbank2db_FILES/tool_output/cfm-id//candidates__cv=${CV}__ion=${ION}.txt"
MODEL_DIR="/data/bach/software/cfm-id-code-r33/supplementary_material/trained_models/esi_msms_models/${ION}_metab_se_cfm/"
PARAM_FN="${MODEL_DIR}/param_output${CV}.log"
CONFIG_FN="${MODEL_DIR}/param_config.txt"
OUTPUT_DIR="/home/bach/Documents/doctoral/projects/massbank2db_FILES/tool_output/cfm-id/predicted_spectra/cv=${CV}__ion=${ION}"

parallel \
  --pipepart \
  --cat \
  --jobs 12 \
  --use-cores-instead-of-threads \
  --eta \
  --block-size -1000 \
  --arg-file "${INPUT_LIST}" \
  "mv {} {}.txt; nice -19 ./cfm-predict {}.txt 0.001 ${PARAM_FN} ${CONFIG_FN} 0 ${OUTPUT_DIR}; rm -f {}.txt"
