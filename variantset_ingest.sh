#!/bin/bash
set -euo pipefail

readonly CANDIG_INGEST_PATH=ingest
readonly CANDIG_REPO_PATH=candig_repo
readonly REGISTRY=registry.db
readonly REFERENCENAME=GRCh37
readonly DATASET=cancogen
readonly VCFPATH="./vcfs"

# initialize the registry and dataset
"${CANDIG_INGEST_PATH}" "${REGISTRY}" "${DATASET}" "${VCFPATH}/ingest.json"

for file in "${VCFPATH}/"*.gz
do
    # add tabix index to the bgzip'ed file
    tabix -p vcf "${file}"

    # get the patient id from the filename
    patient_id=$( basename -s ".vcf.gz" "${file}" )
    
    # register the variantset
    "${CANDIG_REPO_PATH}" add-variantset \
        -I "${file}.tbi" \
        -R "${REFERENCENAME}" \
        "${REGISTRY}" \
        "${DATASET}" \
        "${patient_id}" \
        "${patient_id}_sample_1" \
        "${file}"
done