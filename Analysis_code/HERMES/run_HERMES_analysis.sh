#!/usr/bin/env bash

# Run the HERMES PWAS association + Primo + plotting workflow.
#
# Assumes step 1-3 inputs already exist:
#   - hermes_input/*.rds
#   - hermes_filtered_id/*.bim/*.raw
#
# Defaults use estimated regularization.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
PAPER_DIR="${PAPER_SMIXCAN_DIR:-/Users/zhusinan/Library/CloudStorage/Dropbox/Paper_SMiXcan}"
WORKSPACE_DIR="${HERMES_WORKSPACE_DIR:-${PAPER_DIR}/Results/hermes_pwas/hermes_workspace}"
FIGURE_DIR="${HERMES_FIGURE_DIR:-${PAPER_DIR}/Figure/hermes_pwas}"
RSCRIPT_BIN="${RSCRIPT_BIN:-/Library/Frameworks/R.framework/Resources/bin/Rscript}"
RESULT_TAG="${HERMES_RESULT_TAG:-estimate}"
PLOT_PREFIX="${HERMES_PLOT_PREFIX:-HERMES_PWAS_estimate}"

export HERMES_WORKSPACE_DIR="${WORKSPACE_DIR}"
export HERMES_FIGURE_DIR="${FIGURE_DIR}"
export HERMES_ASSOC_REG_MODE="${HERMES_ASSOC_REG_MODE:-estimate}"
export HERMES_SKIP_EXISTING="${HERMES_SKIP_EXISTING:-true}"
export HERMES_PARALLEL_JOBS="${HERMES_PARALLEL_JOBS:-4}"

cd "${REPO_DIR}"

echo "Step 4: parallel association"
bash Analysis_code/HERMES/run_HERMES_step4_parallel.sh

RESULT_FILE="${WORKSPACE_DIR}/hermes_result/hermes_result_pwas.csv"
ANNOTATED_FILE="${WORKSPACE_DIR}/hermes_result/hermes_result_pwas_${RESULT_TAG}_annotated.csv"

echo "Step 5: Primo"
HERMES_RESULT_FILE="${RESULT_FILE}" \
HERMES_P_JOIN_COL="p_join" \
HERMES_RESULT_TAG="${RESULT_TAG}" \
"${RSCRIPT_BIN}" Analysis_code/HERMES/5_run_Primo_hermes_pwas.R

echo "Step 6: plot"
HERMES_RESULT_FILE="${RESULT_FILE}" \
HERMES_ANNOTATED_FILE="${ANNOTATED_FILE}" \
HERMES_PLOT_PREFIX="${PLOT_PREFIX}" \
"${RSCRIPT_BIN}" Analysis_code/HERMES/6_plot_hermes_pwas.R

echo "Done."
echo "Result: ${RESULT_FILE}"
echo "Annotated: ${ANNOTATED_FILE}"
echo "Figure dir: ${FIGURE_DIR}"
