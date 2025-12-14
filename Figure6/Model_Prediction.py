#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
External dataset prediction with a stacked ensemble:
- Uses precomputed DNA-BERT embeddings from training (mapped by peak_id)
- Applies the saved preprocessor to additional features
- Predicts with 5-fold MLP models + rebuilt ExtraTrees + rebuilt CatBoost
- Combines predictions with saved meta-learner weights
"""

import os
import time
import joblib
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from tqdm import tqdm
from sklearn.ensemble import ExtraTreesRegressor
from catboost import CatBoostRegressor

# =============================================================================
# 0) Logging helper
# =============================================================================
_START_TIME = time.time()
_LAST_STEP_TIME = _START_TIME

def log_step(message: str, major: bool = False) -> None:
    """Print a timestamped log message with elapsed time since last step."""
    global _LAST_STEP_TIME
    now = time.time()
    elapsed = now - _LAST_STEP_TIME

    if major:
        print("\n" + "=" * 80)

    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message} (elapsed: {elapsed:.2f}s)")

    if major:
        print("=" * 80)

    _LAST_STEP_TIME = now


# =============================================================================
# 1) Configuration (VERIFY ALL PATHS)
# =============================================================================
log_step("Initializing configuration", major=True)

# Core input paths
EXTERNAL_DATA_PATH = "/home/yzhao19/m6a/Model/external_3485samples_mlp_new.csv"
ORIGINAL_TRAIN_CSV_PATH = "/home/yzhao19/m6a/Model/training_data_n_5_rppa_mrna_all_12372_seq_TCGAbam.csv"
TRAIN_EMBEDDINGS_NPY_PATH = "/home/yzhao19/m6a/Model/dna_bert_mlp_results_v13/train_embeddings.npy"

# Model artifacts (must match the training run)
MODEL_ARTIFACTS_DIR = "/home/yzhao19/m6a/Model/dna_bert_mlp_results_v16/test_2"

# Output
FINAL_PRED_OUTPUT_PATH = "/home/yzhao19/m6a/Model/Prediction/New/external_3485_predictions_v16.csv"
os.makedirs(os.path.dirname(FINAL_PRED_OUTPUT_PATH), exist_ok=True)

# Model/data settings (must match training)
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Device: {DEVICE}")

HIDDEN_DIM_MLP = 256
OUTPUT_DIM_MLP = 1
DROPOUT_MLP = 0.3
N_MLP_FOLDS = 5

TARGET_COL = "m6A"
SEQUENCE_COLUMN_NAME = "sequence"
IDENTIFIER_COLS = ["peak_id", "gene", "sample_id"]

required_paths = {
    "External dataset CSV": EXTERNAL_DATA_PATH,
    "Original training CSV": ORIGINAL_TRAIN_CSV_PATH,
    "Training embeddings NPY": TRAIN_EMBEDDINGS_NPY_PATH,
    "Model artifacts dir": MODEL_ARTIFACTS_DIR,
}
for name, path in required_paths.items():
    if not os.path.exists(path):
        raise FileNotFoundError(f"{name} not found at: {path}")


# =============================================================================
# 2) Define the MLP (must match training exactly)
# =============================================================================
class MLP(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int, output_dim: int, dropout: float):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),

            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout * 0.8),

            nn.Linear(hidden_dim // 2, output_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.network(x).squeeze(-1)


# =============================================================================
# 3) Load preprocessor and build peak_id -> embedding index map
# =============================================================================
log_step("Loading preprocessor and building peak_id -> embedding index map", major=True)

preprocessor_path = os.path.join(MODEL_ARTIFACTS_DIR, "preprocessor.joblib")
preprocessor = joblib.load(preprocessor_path)
log_step("Loaded preprocessor.joblib")

train_df_for_map = pd.read_csv(ORIGINAL_TRAIN_CSV_PATH, usecols=["peak_id"])
X_train_embed_full = np.load(TRAIN_EMBEDDINGS_NPY_PATH)
log_step(f"Loaded training CSV rows: {len(train_df_for_map):,}")
log_step(f"Loaded training embeddings shape: {X_train_embed_full.shape}")

peak_to_idx_map = (
    train_df_for_map.reset_index()
    .drop_duplicates(subset="peak_id", keep="first")
    .set_index("peak_id")["index"]
    .to_dict()
)
log_step(f"Map size (unique peaks): {len(peak_to_idx_map):,}")
del train_df_for_map


# =============================================================================
# 4) Load and prepare external data
# =============================================================================
log_step("Loading and preparing external dataset", major=True)
external_df = pd.read_csv(EXTERNAL_DATA_PATH)
log_step(f"Loaded external dataset rows: {len(external_df):,}")

missing_peaks = set(external_df["peak_id"].unique()) - set(peak_to_idx_map.keys())
if missing_peaks:
    raise ValueError(
        f"Found {len(missing_peaks)} peak_id(s) in external data without training embeddings. "
        f"Examples: {list(sorted(missing_peaks))[:20]}"
    )

log_step("Matching embeddings for each external row...")
row_indices = [peak_to_idx_map[pid] for pid in tqdm(external_df["peak_id"], desc="Embedding lookup")]
X_external_embed = X_train_embed_full[row_indices]

log_step("Selecting additional feature columns based on training CSV header...")
all_train_cols = pd.read_csv(ORIGINAL_TRAIN_CSV_PATH, nrows=0).columns.tolist()
additional_feature_cols = [
    c for c in all_train_cols
    if c not in IDENTIFIER_COLS + [TARGET_COL, SEQUENCE_COLUMN_NAME]
]

log_step("Transforming additional features with preprocessor...")
X_external_additional_raw = external_df[additional_feature_cols]
X_external_additional = preprocessor.transform(X_external_additional_raw)
log_step(f"Additional features transformed shape: {X_external_additional.shape}")

X_external_comb = np.concatenate([X_external_embed, X_external_additional], axis=1)
log_step(f"Final external feature matrix shape: {X_external_comb.shape}")


# =============================================================================
# 5) Predict with trained models
# =============================================================================
log_step("Running MLP fold models for prediction", major=True)

mlp_preds_folds = []
input_dim_mlp = X_external_comb.shape[1]

for i in tqdm(range(N_MLP_FOLDS), desc="MLP 5-fold prediction"):
    mlp_model = MLP(input_dim_mlp, HIDDEN_DIM_MLP, OUTPUT_DIM_MLP, DROPOUT_MLP).to(DEVICE)
    model_path = os.path.join(MODEL_ARTIFACTS_DIR, f"best_model_A_fold{i}.pth")
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"MLP model file not found: {model_path}")

    mlp_model.load_state_dict(torch.load(model_path, map_location=DEVICE))
    mlp_model.eval()

    with torch.no_grad():
        preds = (
            mlp_model(torch.tensor(X_external_comb, dtype=torch.float32, device=DEVICE))
            .cpu()
            .numpy()
        )
        mlp_preds_folds.append(preds)

final_mlp_preds = np.mean(mlp_preds_folds, axis=0)
log_step("MLP prediction complete")

# ---- Rebuild ExtraTrees and CatBoost using full training data ----
log_step("Rebuilding ExtraTrees and CatBoost from the original training data", major=True)

train_df_orig = pd.read_csv(ORIGINAL_TRAIN_CSV_PATH)
X_train_additional_orig = preprocessor.transform(train_df_orig[additional_feature_cols])
X_train_comb_orig = np.concatenate([X_train_embed_full, X_train_additional_orig], axis=1)
y_train_orig = train_df_orig[TARGET_COL].values
log_step(f"Prepared training matrix for rebuilding: {X_train_comb_orig.shape}")

log_step("Fitting ExtraTreesRegressor (n_estimators=1000)...")
et_full = ExtraTreesRegressor(n_estimators=1000, random_state=42, n_jobs=-1, verbose=0)
et_full.fit(X_train_comb_orig, y_train_orig)
log_step("ExtraTrees fitted; predicting on external data...")
final_et_preds = et_full.predict(X_external_comb)
log_step("ExtraTrees prediction complete")

log_step("Fitting CatBoostRegressor (iterations=1500)...")
cb_params = {
    "iterations": 1500,
    "learning_rate": 0.01,
    "depth": 7,
    "loss_function": "RMSE",
    "random_seed": 42,
}
cb_full = CatBoostRegressor(**cb_params)
cb_full.fit(X_train_comb_orig, y_train_orig, verbose=200)
log_step("CatBoost fitted; predicting on external data...")
final_cb_preds = cb_full.predict(X_external_comb)
log_step("CatBoost prediction complete")

del train_df_orig, X_train_additional_orig, X_train_comb_orig, y_train_orig


# =============================================================================
# 6) Meta-learner weighted fusion
# =============================================================================
log_step("Loading meta-learner weights and fusing predictions", major=True)

weights_path = os.path.join(MODEL_ARTIFACTS_DIR, "meta_model_weights.txt")
if not os.path.exists(weights_path):
    raise FileNotFoundError(f"Meta-learner weights file not found: {weights_path}")

weights = {}
with open(weights_path, "r") as f:
    for line in f:
        key, value = line.strip().split(": ")
        weights[key] = float(value)

log_step(f"Loaded weights: {weights}")

final_predictions = (
    final_mlp_preds * weights["MLP_weight"]
    + final_et_preds * weights["ET_weight"]
    + final_cb_preds * weights["CB_weight"]
    + weights["intercept"]
)
log_step("Final fused predictions computed")


# =============================================================================
# 7) Save results
# =============================================================================
log_step("Saving final prediction output", major=True)

output_df = external_df.copy()
output_df["pred_MLP"] = final_mlp_preds
output_df["pred_ET"] = final_et_preds
output_df["pred_CB"] = final_cb_preds
output_df["pred_Meta_final"] = final_predictions
output_df.to_csv(FINAL_PRED_OUTPUT_PATH, index=False)

total_runtime = time.time() - _START_TIME
log_step(f"Done. Output written to: {FINAL_PRED_OUTPUT_PATH}", major=True)
print(f"Total runtime: {total_runtime / 60:.2f} minutes")
