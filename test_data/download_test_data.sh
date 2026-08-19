#!/usr/bin/env bash
# download_test_data.sh
#
# Download the 5 test assemblies used to validate AssayVal.
# Requires: NCBI datasets CLI (conda install -c conda-forge ncbi-datasets-cli)
#
# Usage: bash test_data/download_test_data.sh
#
# Downloads assemblies to test_data/assemblies/ and builds metadata.csv.
# Then update config/config.yaml:
#   assembly_dir: "test_data/assemblies"
#   metadata: "test_data/assemblies/metadata.csv"
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
OUTDIR="${SCRIPT_DIR}/assemblies"
ACCESSIONS="${SCRIPT_DIR}/accessions.txt"
PARSER="${REPO_ROOT}/scripts/download/parse_metadata.py"

mkdir -p "${OUTDIR}/.tmp"

echo "Downloading 5 test assemblies to ${OUTDIR} ..."

JSONL="${OUTDIR}/.tmp/summary_raw.jsonl"
> "${JSONL}"

while IFS= read -r ACC; do
    TARGET="${OUTDIR}/${ACC}.fna"
    if [[ -s "${TARGET}" ]]; then
        echo "  already exists: ${ACC}"
        continue
    fi

    ACC_TMP="${OUTDIR}/.tmp/${ACC}"
    mkdir -p "${ACC_TMP}"

    echo "  downloading ${ACC} ..."
    datasets download genome accession "${ACC}" \
        --include genome \
        --filename "${ACC_TMP}/download.zip" \
        --no-progressbar

    unzip -q -o "${ACC_TMP}/download.zip" -d "${ACC_TMP}/unzipped"
    # Portable array fill (macOS ships bash 3.2, which lacks `mapfile`)
    FNA_FILES=()
    while IFS= read -r _fna; do FNA_FILES+=("${_fna}"); done \
        < <(find "${ACC_TMP}/unzipped" -name "*.fna" | sort)
    cat "${FNA_FILES[@]}" > "${TARGET}"

    # Collect summary JSON for metadata
    datasets summary genome accession "${ACC}" --as-json-lines >> "${JSONL}"

    rm -rf "${ACC_TMP}"
    sleep 0.3
done < "${ACCESSIONS}"

echo "Building metadata CSV ..."
python3 "${PARSER}" --jsonl "${JSONL}" --out "${OUTDIR}/metadata.csv"

echo ""
echo "Done. Test assemblies are in ${OUTDIR}/"
echo ""
echo "To run AssayVal on the test dataset, edit config/config.yaml:"
echo "  assembly_dir: \"test_data/assemblies\""
echo "  metadata:     \"test_data/assemblies/metadata.csv\""
echo ""
echo "Then run: snakemake --profile workflow/profiles/local"
