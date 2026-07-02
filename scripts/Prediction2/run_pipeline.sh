#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Handle clean flag
if [[ "${1:-}" == "--clean" ]]; then
    rm -rf data/predictathon/*/processed
    rm -rf data/processed/global_union
    rm -rf trained_models
    rm -rf results/cv0_predictions
    rm -rf results/cv00_predictions
    rm -rf submission
    echo "Clean complete."
    exit 0
fi

TRIALS=(
    AWY1_DVPWA_2024
    TCAP_2025_MANKS
    25_Big6_SVREC_SVREC
    OHRWW_2025_SPO
    CornellMaster_2025_McGowan
    24Crk_AY2-3
    2025_AYT_Aurora
    YT_Urb_25
    STP1_2025_MCG
)

# Preprocess genotypes for all trials
for trial in "${TRIALS[@]}"; do
    PROC_DIR="data/predictathon/$trial/processed"
    GRM_PATH="$PROC_DIR/GRM.npy"

    if [ ! -f "$GRM_PATH" ]; then
        python -m src.genotypes.preprocess_genotypes "$trial"
    fi
done

# Build global grm once
GLOBAL_GRM="data/processed/global_union/GRM_global_union.npy"

if [ ! -f "$GLOBAL_GRM" ]; then
    python src/model/build_global_grm_union.py
fi

# Train models and generate predictions
for trial in "${TRIALS[@]}"; do
    python -m src.model.train_model "$trial"

    python "$SCRIPT_DIR/src/model/cv0_predict_global.py" \
        --config "$SCRIPT_DIR/config.yaml" \
        --trial "$trial" \
        --out "$SCRIPT_DIR/results/cv0_predictions/${trial}.csv"

    python "$SCRIPT_DIR/src/model/cv00_predict_global.py" \
        --config "$SCRIPT_DIR/config.yaml" \
        --trial "$trial" \
        --out "$SCRIPT_DIR/results/cv00_predictions/${trial}.csv"
done

# Print completion message
echo "Pipeline complete."