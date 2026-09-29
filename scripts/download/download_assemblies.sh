#!/usr/bin/env bash
# download_assemblies.sh
#
# Download RefSeq genome assemblies for a given taxon from NCBI and build a
# metadata CSV suitable for AssayVal.
#
# Requirements:
#   - NCBI datasets CLI  (conda install -c conda-forge ncbi-datasets-cli)
#   - Python 3 (standard library only)
#   - parse_metadata.py in the same directory as this script
#
# Usage:
#   download-assemblies [OPTIONS]        (installed command; runs from any directory)
#   bash download_assemblies.sh [OPTIONS]  (equivalent, if running the script directly)
#
# Options:
#   -t TAXON       Taxon name or NCBI taxon ID (required). Repeatable: pass -t
#                  multiple times to download the de-duplicated union of several
#                  taxa, e.g. -t "Vibrio jasicida" -t "Vibrio owensii"
#                  Examples: "Vibrionaceae", "Vibrio cholerae", 641
#   -o OUTDIR      Output directory for assemblies and metadata (default: assemblies)
#   -l LEVELS      Assembly levels, comma-separated
#                  (default: complete,chromosome,scaffold,contig)
#   -s SOURCE      Assembly-source: refseq, genbank, or all (default: refseq)
#   -n N           Randomly keep at most N assemblies per -t taxon (default: all)
#   -r SEED        Integer random seed for -n selection (default: 0)
#   -e EMAIL       NCBI e-mail address (optional but polite; or set NCBI_EMAIL env var)
#   -k API_KEY     NCBI API key for higher rate limits (or set NCBI_API_KEY env var)
#   -h             Show this help message
#
# Set your API key once (recommended): copy config/ncbi_credentials.example.sh to
# config/ncbi_credentials.sh and paste your key. This script sources it on every
# run. Key precedence: -k flag > NCBI_API_KEY env var > credentials file.
#
# On HPC/SLURM: wrap the command in an sbatch job. Example:
#   sbatch --time=24:00:00 --mem=8G --wrap="download-assemblies -t Vibrionaceae -o assemblies"
#
# Output:
#   OUTDIR/*.fna          — one FASTA file per assembly (multi-replicon genomes concatenated)
#   OUTDIR/metadata.csv   — assembly metadata for AssayVal's metadata config key
#
set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
TAXA=()
OUTDIR="assemblies"
LEVELS="complete,chromosome,scaffold,contig"
ASSEMBLY_SOURCE="refseq"
LIMIT=""       # empty = keep all (no subsampling)
SEED=0

# ── Credentials ("set once") ────────────────────────────────────────────────────
# Load an optional credentials file so an NCBI API key is applied on every run.
# Override its location with the ASSAYVAL_CREDENTIALS environment variable.
# Key precedence: -k flag > NCBI_API_KEY env var > credentials file.
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
_REPO_ROOT="$(cd "${_SCRIPT_DIR}/../.." && pwd)"
CREDENTIALS_FILE="${ASSAYVAL_CREDENTIALS:-${_REPO_ROOT}/config/ncbi_credentials.sh}"
_ENV_API_KEY="${NCBI_API_KEY:-}"
_ENV_EMAIL="${NCBI_EMAIL:-}"
if [[ -f "${CREDENTIALS_FILE}" ]]; then
    # shellcheck source=/dev/null
    source "${CREDENTIALS_FILE}"
fi
# A pre-existing environment variable takes precedence over the file's value.
EMAIL="${_ENV_EMAIL:-${NCBI_EMAIL:-}}"
API_KEY="${_ENV_API_KEY:-${NCBI_API_KEY:-}}"

# ── Argument parsing ───────────────────────────────────────────────────────────
usage() {
    # Print the header comment block from "# Usage:" to the first non-comment
    # line, stripping the leading "# ". awk is used (not GNU-sed block syntax)
    # so this works with the BSD sed/awk shipped on macOS.
    awk '/^# Usage:/{p=1}
         p && /^[^#]/{exit}
         p{sub(/^# ?/, ""); print}' "$0"
    exit 1
}

while getopts "t:o:l:s:n:r:e:k:h" opt; do
    case $opt in
        t) TAXA+=("$OPTARG") ;;
        o) OUTDIR="$OPTARG" ;;
        l) LEVELS="$OPTARG" ;;
        s) ASSEMBLY_SOURCE="$OPTARG" ;;
        n) LIMIT="$OPTARG" ;;
        r) SEED="$OPTARG" ;;
        e) EMAIL="$OPTARG" ;;
        k) API_KEY="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

if [[ ${#TAXA[@]} -eq 0 ]]; then
    echo "ERROR: at least one -t TAXON is required." >&2
    usage
fi

case "${ASSEMBLY_SOURCE}" in
    refseq|genbank|all) ;;
    *) echo "ERROR: -s must be one of: refseq, genbank, all (got '${ASSEMBLY_SOURCE}')." >&2; exit 1 ;;
esac

if [[ -n "${LIMIT}" ]]; then
    if ! [[ "${LIMIT}" =~ ^[0-9]+$ ]] || [[ "$((10#${LIMIT}))" -lt 1 ]]; then
        echo "ERROR: -n must be a positive integer (got '${LIMIT}')." >&2
        exit 1
    fi
fi
if ! [[ "${SEED}" =~ ^-?[0-9]+$ ]]; then
    echo "ERROR: -r (seed) must be an integer (got '${SEED}')." >&2
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PARSER="${SCRIPT_DIR}/parse_metadata.py"
if [[ ! -f "$PARSER" ]]; then
    echo "ERROR: parse_metadata.py not found at ${PARSER}" >&2
    exit 1
fi
SUBSAMPLE="${SCRIPT_DIR}/subsample.py"
if [[ -n "${LIMIT}" && ! -f "${SUBSAMPLE}" ]]; then
    echo "ERROR: subsample.py not found at ${SUBSAMPLE}" >&2
    exit 1
fi

# ── Setup ──────────────────────────────────────────────────────────────────────
WORK_TMP="${OUTDIR}/.tmp"
ACCESSION_LIST="${WORK_TMP}/accessions.txt"
SUMMARY_ALL="${WORK_TMP}/summary_all.jsonl"
JSONL_RAW="${WORK_TMP}/summary_raw.jsonl"
METADATA_CSV="${OUTDIR}/metadata.csv"
FAILED_LIST="${WORK_TMP}/failed_accessions.txt"
LOG="${OUTDIR}/download.log"

mkdir -p "${OUTDIR}" "${WORK_TMP}"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG}"; }

log "===== AssayVal assembly download ====="
log "Taxa         : ${#TAXA[@]} (${TAXA[*]})"
log "Output dir   : ${OUTDIR}"
log "Assembly levels: ${LEVELS}"
log "Assembly source: ${ASSEMBLY_SOURCE}"
if [[ -n "${LIMIT}" ]]; then
    log "Subsampling  : up to ${LIMIT} per taxon (seed=${SEED})"
fi

# ── Step 1: Fetch accession list ───────────────────────────────────────────────
log "Fetching assembly list from NCBI ..."

# Apply the API key (if any) to every datasets call. Empty-array expansion is
# guarded so it is safe under `set -u` on bash 3.2 (macOS default).
APIKEY_ARGS=()
if [[ -n "${API_KEY}" ]]; then
    APIKEY_ARGS=(--api-key "${API_KEY}")
    log "Using NCBI API key (higher rate limit)."
fi

# Query each taxon. When -n is set, randomly cap each taxon's summaries via
# subsample.py before accumulating. Overlapping taxa are de-duplicated by
# accession in the next step (so a genome shared between taxa is kept once).
: > "${SUMMARY_ALL}"
for TAX in "${TAXA[@]}"; do
    log "  querying: ${TAX}"
    TAX_JSONL="${WORK_TMP}/taxon_summary.jsonl"
    datasets summary genome taxon "${TAX}" \
        --assembly-source "${ASSEMBLY_SOURCE}" \
        --assembly-level "${LEVELS}" \
        --as-json-lines \
        ${APIKEY_ARGS[@]+"${APIKEY_ARGS[@]}"} > "${TAX_JSONL}"
    AVAIL=$(grep -c . "${TAX_JSONL}" 2>/dev/null || true)
    if [[ -n "${LIMIT}" ]]; then
        TAX_SAMPLED="${WORK_TMP}/taxon_sampled.jsonl"
        python3 "${SUBSAMPLE}" --jsonl "${TAX_JSONL}" --n "${LIMIT}" \
            --seed "${SEED}" --out "${TAX_SAMPLED}"
        KEPT=$(grep -c . "${TAX_SAMPLED}" 2>/dev/null || true)
        log "  ${TAX}: sampled ${KEPT} of ${AVAIL} available (seed=${SEED})"
        cat "${TAX_SAMPLED}" >> "${SUMMARY_ALL}"
    else
        log "  ${TAX}: ${AVAIL} available (no subsampling)"
        cat "${TAX_JSONL}" >> "${SUMMARY_ALL}"
    fi
done

# De-duplicate by accession → canonical JSONL (for metadata) + unique accession
# list (for downloading). Keeps the first occurrence of each accession.
python3 - "${SUMMARY_ALL}" "${JSONL_RAW}" "${ACCESSION_LIST}" << 'PYEOF'
import sys, json
in_path, jsonl_out, acc_out = sys.argv[1], sys.argv[2], sys.argv[3]
seen = set()
kept_lines = []
accessions = []
with open(in_path) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        try:
            acc = json.loads(line).get("accession", "")
        except json.JSONDecodeError:
            continue
        if not acc or acc in seen:
            continue
        seen.add(acc)
        accessions.append(acc)
        kept_lines.append(line)
with open(jsonl_out, "w") as out:
    out.write("\n".join(kept_lines))
    if kept_lines:
        out.write("\n")
with open(acc_out, "w") as out:
    out.writelines(a + "\n" for a in accessions)
print(f"Found {len(accessions)} unique accessions.", flush=True)
PYEOF

TOTAL=$(wc -l < "${ACCESSION_LIST}")
log "Total assemblies found: ${TOTAL}"

# ── Step 2: Download genome FASTA files (resume-aware) ────────────────────────
log "Downloading genome FASTA files ..."

DOWNLOADED=0; SKIPPED=0; FAILED=0
> "${FAILED_LIST}"

while IFS= read -r ACC; do
    TARGET="${OUTDIR}/${ACC}.fna"

    if [[ -s "${TARGET}" ]]; then
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    ACC_TMP="${WORK_TMP}/${ACC}"
    mkdir -p "${ACC_TMP}"

    if datasets download genome accession "${ACC}" \
            --include genome \
            --filename "${ACC_TMP}/download.zip" \
            --no-progressbar \
            ${APIKEY_ARGS[@]+"${APIKEY_ARGS[@]}"} \
            2>> "${LOG}"; then

        unzip -q -o "${ACC_TMP}/download.zip" -d "${ACC_TMP}/unzipped" 2>> "${LOG}"
        # Portable array fill (macOS ships bash 3.2, which lacks `mapfile`)
        FNA_FILES=()
        while IFS= read -r _fna; do FNA_FILES+=("${_fna}"); done \
            < <(find "${ACC_TMP}/unzipped" -name "*.fna" | sort)

        if [[ ${#FNA_FILES[@]} -eq 0 ]]; then
            log "WARNING: no .fna files found for ${ACC}"
            echo "${ACC}" >> "${FAILED_LIST}"
            FAILED=$((FAILED + 1))
        else
            cat "${FNA_FILES[@]}" > "${TARGET}"
            DOWNLOADED=$((DOWNLOADED + 1))
            log "OK  ${ACC}  ($(grep -c '^>' "${TARGET}") contig(s))"
        fi
    else
        log "WARNING: download failed for ${ACC}"
        echo "${ACC}" >> "${FAILED_LIST}"
        FAILED=$((FAILED + 1))
    fi

    rm -rf "${ACC_TMP}"
    sleep 0.3   # polite pause for NCBI

done < "${ACCESSION_LIST}"

log "Downloads complete — Downloaded: ${DOWNLOADED}  Skipped: ${SKIPPED}  Failed: ${FAILED}"
[[ -s "${FAILED_LIST}" ]] && log "Failed accessions: ${FAILED_LIST} (re-run to retry)"

# ── Step 3: Build metadata CSV ─────────────────────────────────────────────────
log "Building metadata CSV ..."

python3 "${PARSER}" --jsonl "${JSONL_RAW}" --out "${METADATA_CSV}" --merge

META_ROWS=$(( $(wc -l < "${METADATA_CSV}") - 1 ))
if [[ "${META_ROWS}" -lt 1 ]]; then
    log "ERROR: metadata CSV has no data rows."
    exit 1
fi
log "Metadata written: ${META_ROWS} rows → ${METADATA_CSV}"

# ── Done ───────────────────────────────────────────────────────────────────────
log "===== Done ====="
log "FNA files on disk : $(find "${OUTDIR}" -maxdepth 1 -name '*.fna' | wc -l)"
log "Metadata CSV      : ${METADATA_CSV}"
log ""
log "Next, from your analysis directory:"
log "  assay-val init            # if you have not yet: drops config.yaml + assay_table.csv"
log "  (edit assay_table.csv; check assembly_dir and metadata in config.yaml)"
log "  assay-val --run-name <name>"
