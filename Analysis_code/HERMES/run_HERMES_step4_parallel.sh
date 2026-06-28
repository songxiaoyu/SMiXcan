#!/usr/bin/env bash

# Run HERMES PWAS step 4 by chromosome in parallel.
#
# Defaults:
#   - Runs chromosomes 1-22.
#   - Starts up to 4 chromosomes at a time.
#   - Skips chromosome result files that already exist.
#   - Writes the combined hermes_result_pwas.csv only after all workers finish.
#
# Useful overrides:
#   HERMES_WORKSPACE_DIR=/path/to/workspace
#   HERMES_CHR_LIST=1,2,3
#   HERMES_PARALLEL_JOBS=4
#   HERMES_SKIP_EXISTING=true
#   HERMES_ASSOC_REG_MODE=fixed
#   HERMES_ASSOC_REG_SCALES=0.005

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PAPER_DIR="${PAPER_SMIXCAN_DIR:-/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan}"
WORKSPACE_DIR="${HERMES_WORKSPACE_DIR:-${PAPER_DIR}/Results/hermes_pwas/hermes_workspace}"
RESULT_DIR="${WORKSPACE_DIR}/hermes_result"
LOG_DIR="${WORKSPACE_DIR}/logs/step4_parallel"
RSCRIPT_BIN="${RSCRIPT_BIN:-/Library/Frameworks/R.framework/Resources/bin/Rscript}"
JOBS="${HERMES_PARALLEL_JOBS:-4}"
CHR_LIST_RAW="${HERMES_CHR_LIST:-$(printf '%s,' {1..22} | sed 's/,$//')}"
CHR_LIST="${CHR_LIST_RAW//,/ }"

mkdir -p "${RESULT_DIR}" "${LOG_DIR}"

if ! [[ "${JOBS}" =~ ^[0-9]+$ ]] || [ "${JOBS}" -lt 1 ]; then
  echo "HERMES_PARALLEL_JOBS must be a positive integer." >&2
  exit 1
fi

run_chr() {
  local chr="$1"
  local log_file="${LOG_DIR}/chr${chr}.log"
  echo "Starting chr${chr}; log: ${log_file}"
  (
    cd "${REPO_DIR}"
    HERMES_CHR_LIST="${chr}" \
    HERMES_WORKSPACE_DIR="${WORKSPACE_DIR}" \
    HERMES_SKIP_EXISTING="${HERMES_SKIP_EXISTING:-true}" \
    HERMES_WRITE_COMBINED=false \
    "${RSCRIPT_BIN}" Analysis_code/HERMES/4_HERMES_run_analysis_pwas.R
  ) >"${log_file}" 2>&1
}

for chr in ${CHR_LIST}; do
  while [ "$(jobs -pr | wc -l | tr -d ' ')" -ge "${JOBS}" ]; do
    sleep 5
  done
  run_chr "${chr}" &
done

status=0
for pid in $(jobs -pr); do
  if ! wait "${pid}"; then
    status=1
  fi
done

if [ "${status}" -ne 0 ]; then
  echo "At least one chromosome failed. Check logs in ${LOG_DIR}." >&2
  exit "${status}"
fi

echo "Combining chromosome results..."
(
  cd "${REPO_DIR}"
  HERMES_CHR_LIST="${CHR_LIST_RAW}" \
  HERMES_WORKSPACE_DIR="${WORKSPACE_DIR}" \
  HERMES_SKIP_EXISTING=true \
  HERMES_WRITE_COMBINED=true \
  HERMES_MAX_GENES=0 \
  "${RSCRIPT_BIN}" -e '
    result_dir <- file.path(Sys.getenv("HERMES_WORKSPACE_DIR"), "hermes_result")
    chr_list <- as.integer(strsplit(Sys.getenv("HERMES_CHR_LIST"), ",")[[1]])
    files <- file.path(result_dir, sprintf("hermes_chr%d_result_pwas.csv", chr_list))
    missing <- files[!file.exists(files)]
    if (length(missing)) {
      stop("Missing chromosome result files:\n", paste(missing, collapse = "\n"), call. = FALSE)
    }
    combined <- do.call(rbind, lapply(files, read.csv))
    out <- file.path(result_dir, "hermes_result_pwas.csv")
    write.csv(combined, out, row.names = FALSE)
    cat("Wrote:", out, "\n")
    cat("Rows:", nrow(combined), "\n")
  '
)

echo "Done."
