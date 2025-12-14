#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Train a stacked ensemble for m6A prediction:
- Inputs: precomputed DNA-BERT embeddings + additional tabular features
- Base models: 5-fold MLP, ExtraTrees, CatBoost
- Meta learner: Ridge on OOF predictions
- Saves: preprocessor, fold MLP models, OOF predictions, meta model, weights, test predictions,
         full ET/CB models, feature importances
"""

# =============================================================================
# 0) Headless matplotlib setup
# =============================================================================
import matplotlib
matplotlib.use("Agg")  # Enables saving figures on servers without a display

# =============================================================================
# 1) Imports
# =============================================================================
import os
import gc
import numpy as np
import pandas as pd
from tqdm import tqdm

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from torch.optim.lr_scheduler import ReduceLROnPlateau

from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score

from catboost import CatBoostRegressor
from scipy.stats import pearsonr

import joblib

# =============================================================================
# 2) Configuration (update paths as needed)
# =============================================================================
BERT_MODEL_NAME = "zhihan1996/DNA_bert_6"  # Not used here (embeddings are precomputed)
MAX_SEQ_LENGTH = 201                      # Not used here (embeddings are precomputed)

BATCH_SIZE_MLP = 128
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# MLP & stacking hyperparameters
EMBEDDING_DIM = 768
HIDDEN_DIM_MLP = 256
OUTPUT_DIM_MLP = 1
DROPOUT_MLP = 0.3
LEARNING_RATE_MLP = 3e-4
EPOCHS_MLP = 150
EARLY_STOPPING_PATIENCE = 5
MIN_DELTA = 1e-4

N_FOLDS = 5

# Data paths
TRAIN_CSV_PATH = "/home/yzhao19/m6a/Model/training_data_n_5_rppa_mrna_all_12372_seq_TCGAbam.csv"
TEST_CSV_PATH  = "/home/yzhao19/m6a/Model/test_data_n_5_rppa_mrna_all_12372_seq_TCGAbam.csv"
TRAIN_EMBED_PATH = "/home/yzhao19/m6a/Model/dna_bert_mlp_results_v13/train_embeddings.npy"
TEST_EMBED_PATH  = "/home/yzhao19/m6a/Model/dna_bert_mlp_results_v13/test_embeddings.npy"

SEQUENCE_COLUMN_NAME = "sequence"
IDENTIFIER_COLS = ["peak_id", "gene", "sample_id"]
TARGET_COL = "m6A"

# Output directory
OUTPUT_DIR = "/home/yzhao19/m6a/Model/dna_bert_mlp_results_v16/test_2"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Columns to remove before modeling
COLS_TO_REMOVE = ["has_cluster50", "CancerType", "DEL", "drach_count"]


# =============================================================================
# 3) Load data and remove specified columns
# =============================================================================
print("Loading data...")
train_df = pd.read_csv(TRAIN_CSV_PATH)
test_df = pd.read_csv(TEST_CSV_PATH)
print(f"Original training data shape: {train_df.shape}")
print(f"Original test     data shape: {test_df.shape}")

print(f"\nRemoving columns (if present): {COLS_TO_REMOVE}")
train_df = train_df.drop(columns=COLS_TO_REMOVE, errors="ignore")
test_df  = test_df.drop(columns=COLS_TO_REMOVE, errors="ignore")

y_train = train_df[TARGET_COL].values
y_test  = test_df[TARGET_COL].values
print(f"\nTraining data shape after removal: {train_df.shape}")
print(f"Test     data shape after removal: {test_df.shape}\n")


# =============================================================================
# 4) Prepare additional features and preprocessing
# =============================================================================
all_cols = train_df.columns.tolist()
additional_feature_cols = [
    col for col in all_cols
    if col not in IDENTIFIER_COLS + [TARGET_COL, SEQUENCE_COLUMN_NAME]
]

X_train_additional_raw = train_df[additional_feature_cols].copy()
X_test_additional_raw  = test_df[additional_feature_cols].copy()

# Numerical / categorical split
all_numerical_cols = X_train_additional_raw.select_dtypes(include=np.number).columns.tolist()
categorical_cols   = X_train_additional_raw.select_dtypes(exclude=np.number).columns.tolist()

# Drop all-NaN columns
null_counts = X_train_additional_raw.isnull().sum()
n_rows = X_train_additional_raw.shape[0]
all_nan_cols = null_counts[null_counts == n_rows].index.tolist()
if all_nan_cols:
    print("Warning: dropping columns that are 100% NaN:", all_nan_cols)
    X_train_additional_raw = X_train_additional_raw.drop(columns=all_nan_cols)
    X_test_additional_raw  = X_test_additional_raw.drop(columns=all_nan_cols)
    all_numerical_cols = X_train_additional_raw.select_dtypes(include=np.number).columns.tolist()
    categorical_cols   = X_train_additional_raw.select_dtypes(exclude=np.number).columns.tolist()

# Identify zero-variance numerical columns
variances = X_train_additional_raw[all_numerical_cols].var(skipna=True)
numerical_cols_to_scale = variances[variances > 0].index.tolist()
numerical_cols_constant = variances[variances == 0].index.tolist()

num_scale_pipeline = Pipeline(steps=[
    ("imputer", SimpleImputer(strategy="median")),
    ("scaler", StandardScaler()),
])
num_const_pipeline = Pipeline(steps=[
    ("imputer", SimpleImputer(strategy="median")),
])
cat_pipeline = Pipeline(steps=[
    ("imputer", SimpleImputer(strategy="most_frequent")),
    ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
])

transformers_list = []
if numerical_cols_to_scale:
    transformers_list.append(("num_scale", num_scale_pipeline, numerical_cols_to_scale))
if numerical_cols_constant:
    transformers_list.append(("num_const", num_const_pipeline, numerical_cols_constant))
if categorical_cols:
    transformers_list.append(("cat", cat_pipeline, categorical_cols))

preprocessor = ColumnTransformer(transformers=transformers_list, remainder="passthrough")

print("Fitting ColumnTransformer ...")
preprocessor.fit(X_train_additional_raw)

preprocessor_path = os.path.join(OUTPUT_DIR, "preprocessor.joblib")
joblib.dump(preprocessor, preprocessor_path)
print(f"Preprocessor saved to: {preprocessor_path}")

X_train_additional = preprocessor.transform(X_train_additional_raw)
X_test_additional  = preprocessor.transform(X_test_additional_raw)
print(f"Shapes after preprocess: train={X_train_additional.shape}, test={X_test_additional.shape}")

try:
    processed_feature_names = preprocessor.get_feature_names_out()
except Exception:
    processed_feature_names = None


# =============================================================================
# 5) Load precomputed embeddings
# =============================================================================
print("\nLoading precomputed embeddings...")
X_train_embed = np.load(TRAIN_EMBED_PATH)
X_test_embed  = np.load(TEST_EMBED_PATH)
print(f"Train embeddings shape: {X_train_embed.shape}")
print(f"Test  embeddings shape: {X_test_embed.shape}")


# =============================================================================
# 6) Define MLP and Pearson correlation loss
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


def pearson_corr_loss(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    x_mean, y_mean = torch.mean(x), torch.mean(y)
    xm, ym = x - x_mean, y - y_mean
    cov = torch.mean(xm * ym)
    x_std = torch.sqrt(torch.mean(xm ** 2) + 1e-8)
    y_std = torch.sqrt(torch.mean(ym ** 2) + 1e-8)
    return -cov / (x_std * y_std + 1e-8)


class CombinedDataset(Dataset):
    def __init__(self, embed: np.ndarray, add_feats: np.ndarray, labels: np.ndarray):
        self.embed = torch.tensor(embed, dtype=torch.float32)
        self.add = torch.tensor(add_feats, dtype=torch.float32)
        self.labels = torch.tensor(labels, dtype=torch.float32)

    def __len__(self) -> int:
        return len(self.labels)

    def __getitem__(self, idx: int):
        return self.embed[idx], self.add[idx], self.labels[idx]


# =============================================================================
# 7) Prepare containers for stacking
# =============================================================================
N_train = X_train_embed.shape[0]

oof_preds_A = np.zeros(N_train, dtype=np.float32)  # MLP
oof_preds_B = np.zeros(N_train, dtype=np.float32)  # ExtraTrees
oof_preds_C = np.zeros(N_train, dtype=np.float32)  # CatBoost

test_preds_A_folds = []
test_preds_B_folds = []
test_preds_C_folds = []


# =============================================================================
# 8) 5-fold stacking training
# =============================================================================
kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=42)

for fold_idx, (train_idx, val_idx) in enumerate(kf.split(np.arange(N_train))):
    print(f"\n{'=' * 24} Fold {fold_idx + 1}/{N_FOLDS} {'=' * 24}\n")

    # Split data
    X_tr_embed, X_vl_embed = X_train_embed[train_idx], X_train_embed[val_idx]
    X_tr_add, X_vl_add = X_train_additional[train_idx], X_train_additional[val_idx]
    y_tr, y_vl = y_train[train_idx], y_train[val_idx]

    X_tr_comb = np.concatenate((X_tr_embed, X_tr_add), axis=1)
    X_vl_comb = np.concatenate((X_vl_embed, X_vl_add), axis=1)
    X_test_comb = np.concatenate((X_test_embed, X_test_additional), axis=1)

    # ---------------------------
    # 8.1 Base model A: MLP
    # ---------------------------
    print("Training base model A: MLP")

    mlp_model = MLP(
        input_dim=X_tr_comb.shape[1],
        hidden_dim=HIDDEN_DIM_MLP,
        output_dim=OUTPUT_DIM_MLP,
        dropout=DROPOUT_MLP,
    ).to(DEVICE)

    optimizer = optim.AdamW(mlp_model.parameters(), lr=LEARNING_RATE_MLP, weight_decay=1e-5)
    scheduler = ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=2, min_lr=1e-6, verbose=False)
    criterion = pearson_corr_loss

    train_loader = DataLoader(
        CombinedDataset(X_tr_embed, X_tr_add, y_tr),
        batch_size=BATCH_SIZE_MLP,
        shuffle=True,
        num_workers=2,
        pin_memory=True,
    )
    val_loader = DataLoader(
        CombinedDataset(X_vl_embed, X_vl_add, y_vl),
        batch_size=BATCH_SIZE_MLP,
        shuffle=False,
        num_workers=2,
        pin_memory=True,
    )

    best_val_loss = float("inf")
    epochs_no_improve = 0
    fold_model_path = os.path.join(OUTPUT_DIR, f"best_model_A_fold{fold_idx}.pth")

    for epoch in range(EPOCHS_MLP):
        # Train
        mlp_model.train()
        train_loss = 0.0
        for emb, add, yb in train_loader:
            inputs = torch.cat([emb.to(DEVICE), add.to(DEVICE)], dim=1)
            targets = yb.to(DEVICE)

            optimizer.zero_grad()
            preds = mlp_model(inputs)
            loss = criterion(preds, targets)
            loss.backward()
            optimizer.step()

            train_loss += loss.item() * len(yb)

        epoch_train_loss = train_loss / len(train_loader.dataset)

        # Validate
        mlp_model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for emb, add, yb in val_loader:
                inputs = torch.cat([emb.to(DEVICE), add.to(DEVICE)], dim=1)
                targets = yb.to(DEVICE)
                preds = mlp_model(inputs)
                val_loss += criterion(preds, targets).item() * len(yb)

        epoch_val_loss = val_loss / len(val_loader.dataset)

        print(
            f"Epoch {epoch + 1:3d}/{EPOCHS_MLP} | "
            f"TrainLoss={epoch_train_loss:.5f} | ValLoss={epoch_val_loss:.5f}"
        )
        scheduler.step(epoch_val_loss)

        # Early stopping
        if epoch_val_loss < best_val_loss - MIN_DELTA:
            best_val_loss = epoch_val_loss
            epochs_no_improve = 0
            torch.save(mlp_model.state_dict(), fold_model_path)
        else:
            epochs_no_improve += 1
            if epochs_no_improve >= EARLY_STOPPING_PATIENCE:
                print(f"Early stopping triggered at epoch {epoch + 1}")
                break

    print(f"MLP fold complete. Best model saved to: {fold_model_path}")

    # Load best model
    mlp_model.load_state_dict(torch.load(fold_model_path, map_location=DEVICE))
    mlp_model.eval()

    # OOF prediction for validation indices
    with torch.no_grad():
        vl_inputs = torch.cat(
            [
                torch.tensor(X_vl_embed, dtype=torch.float32, device=DEVICE),
                torch.tensor(X_vl_add, dtype=torch.float32, device=DEVICE),
            ],
            dim=1,
        )
        oof_preds_A[val_idx] = mlp_model(vl_inputs).cpu().numpy()

        # Test prediction for this fold
        test_dataset = CombinedDataset(X_test_embed, X_test_additional, np.zeros(X_test_embed.shape[0], dtype=np.float32))
        test_loader = DataLoader(test_dataset, batch_size=BATCH_SIZE_MLP, shuffle=False)
        preds_test_A = []
        for emb, add, _ in test_loader:
            inputs = torch.cat([emb.to(DEVICE), add.to(DEVICE)], dim=1)
            preds_test_A.extend(mlp_model(inputs).cpu().numpy())
        test_preds_A_folds.append(np.array(preds_test_A))

    # ---------------------------
    # 8.2 Base model B: ExtraTrees
    # ---------------------------
    print("\nTraining base model B: ExtraTreesRegressor")
    et = ExtraTreesRegressor(n_estimators=1000, random_state=42, n_jobs=-1, verbose=0)
    et.fit(X_tr_comb, y_tr)
    oof_preds_B[val_idx] = et.predict(X_vl_comb)
    test_preds_B_folds.append(et.predict(X_test_comb))
    print("ExtraTrees fold complete.")

    # ---------------------------
    # 8.3 Base model C: CatBoost
    # ---------------------------
    print("\nTraining base model C: CatBoostRegressor")
    cb_params = {
        "iterations": 1500,
        "learning_rate": 0.01,
        "depth": 7,
        "loss_function": "RMSE",
        "random_seed": 42,
    }
    cb = CatBoostRegressor(**cb_params)
    cb.fit(X_tr_comb, y_tr, eval_set=(X_vl_comb, y_vl), early_stopping_rounds=50, verbose=0)
    oof_preds_C[val_idx] = cb.predict(X_vl_comb)
    test_preds_C_folds.append(cb.predict(X_test_comb))
    print("CatBoost fold complete.")

    # Cleanup
    del mlp_model, et, cb, X_tr_comb, X_vl_comb
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


# =============================================================================
# 9) Save OOF predictions and meta-features
# =============================================================================
print("\nCreating and saving OOF predictions and meta-features...")

oof_df = pd.DataFrame({
    "oof_pred_MLP": oof_preds_A,
    "oof_pred_ET":  oof_preds_B,
    "oof_pred_CB":  oof_preds_C,
    "true_m6A":     y_train,
})
oof_path = os.path.join(OUTPUT_DIR, "oof_predictions.csv")
oof_df.to_csv(oof_path, index=False)
print(f"OOF predictions saved to: {oof_path}")

X_meta_train = oof_df[["oof_pred_MLP", "oof_pred_ET", "oof_pred_CB"]].values

avg_test_pred_A = np.mean(np.vstack(test_preds_A_folds), axis=0)
avg_test_pred_B = np.mean(np.vstack(test_preds_B_folds), axis=0)
avg_test_pred_C = np.mean(np.vstack(test_preds_C_folds), axis=0)
X_meta_test = np.vstack([avg_test_pred_A, avg_test_pred_B, avg_test_pred_C]).T


# =============================================================================
# 10) Train and save meta-learner (Ridge)
# =============================================================================
print("\nTraining and saving meta-learner: Ridge")
meta_model = Ridge(alpha=1e-3, random_state=42)
meta_model.fit(X_meta_train, y_train)

meta_model_path = os.path.join(OUTPUT_DIR, "meta_model_ridge.joblib")
joblib.dump(meta_model, meta_model_path)
print(f"Meta-learner saved to: {meta_model_path}")

weights_path = os.path.join(OUTPUT_DIR, "meta_model_weights.txt")
with open(weights_path, "w") as f:
    f.write(f"MLP_weight: {meta_model.coef_[0]:.6f}\n")
    f.write(f"ET_weight:  {meta_model.coef_[1]:.6f}\n")
    f.write(f"CB_weight:  {meta_model.coef_[2]:.6f}\n")
    f.write(f"intercept:  {meta_model.intercept_:.6f}\n")
print(f"Meta-learner weights saved to: {weights_path}")


# =============================================================================
# 11) Evaluate and save final test predictions
# =============================================================================
print("\nFinal evaluation...")
final_test_pred = meta_model.predict(X_meta_test)

r_oof_A = pearsonr(oof_preds_A, y_train)[0]
r_oof_B = pearsonr(oof_preds_B, y_train)[0]
r_oof_C = pearsonr(oof_preds_C, y_train)[0]

r_test_A = pearsonr(avg_test_pred_A, y_test)[0]
r_test_B = pearsonr(avg_test_pred_B, y_test)[0]
r_test_C = pearsonr(avg_test_pred_C, y_test)[0]

r_test_meta = pearsonr(final_test_pred, y_test)[0]
mse_meta = mean_squared_error(y_test, final_test_pred)
r2_meta = r2_score(y_test, final_test_pred)

print(f"OOF Pearson r:  MLP={r_oof_A:.4f}, ET={r_oof_B:.4f}, CB={r_oof_C:.4f}")
print(f"Test Pearson r: MLP={r_test_A:.4f}, ET={r_test_B:.4f}, CB={r_test_C:.4f}")
print(f"Test stacking:  Pearson r={r_test_meta:.4f}, RMSE={np.sqrt(mse_meta):.5f}, R2={r2_meta:.5f}")

test_df["pred_MLP"] = avg_test_pred_A
test_df["pred_ET"] = avg_test_pred_B
test_df["pred_CB"] = avg_test_pred_C
test_df["pred_Meta"] = final_test_pred

final_pred_path = os.path.join(OUTPUT_DIR, "test_peak_sample_predictions.csv")
test_df[["peak_id", "sample_id", TARGET_COL, "pred_MLP", "pred_ET", "pred_CB", "pred_Meta"]].to_csv(
    final_pred_path, index=False
)
print(f"Final test predictions saved to: {final_pred_path}")


# =============================================================================
# 12) Save full models for inference
# =============================================================================
print("\nTraining and saving full models for inference...")
X_all_comb = np.concatenate((X_train_embed, X_train_additional), axis=1)

print("Training full ExtraTrees model...")
et_full = ExtraTreesRegressor(n_estimators=1000, random_state=42, n_jobs=-1).fit(X_all_comb, y_train)
et_full_path = os.path.join(OUTPUT_DIR, "et_full.joblib")
joblib.dump(et_full, et_full_path)
print(f"Full ExtraTrees model saved to: {et_full_path}")

print("Training full CatBoost model...")
cb_full = CatBoostRegressor(**cb_params).fit(X_all_comb, y_train, verbose=0)
cb_full_path = os.path.join(OUTPUT_DIR, "cb_full.cbm")
cb_full.save_model(cb_full_path)
print(f"Full CatBoost model saved to: {cb_full_path}")

print("MLP fold models were already saved during cross-validation.")


# =============================================================================
# 13) Feature importances (optional)
# =============================================================================
print("\nComputing feature importances on full training data...")

if processed_feature_names is not None:
    combined_feature_names = [f"embed_{i}" for i in range(EMBEDDING_DIM)] + list(processed_feature_names)
else:
    combined_feature_names = [f"embed_{i}" for i in range(EMBEDDING_DIM)] + [
        f"feat_{i}" for i in range(X_train_additional.shape[1])
    ]

importance_et_df = (
    pd.DataFrame({"feature": combined_feature_names, "importance": et_full.feature_importances_})
    .sort_values("importance", ascending=False)
)
importance_et_path = os.path.join(OUTPUT_DIR, "feature_importances_ET.csv")
importance_et_df.to_csv(importance_et_path, index=False)
print(f"Saved ExtraTrees feature importances to: {importance_et_path}")

importance_cb_df = (
    pd.DataFrame({"feature": combined_feature_names, "importance": cb_full.get_feature_importance()})
    .sort_values("importance", ascending=False)
)
importance_cb_path = os.path.join(OUTPUT_DIR, "feature_importances_CB.csv")
importance_cb_df.to_csv(importance_cb_path, index=False)
print(f"Saved CatBoost feature importances to: {importance_cb_path}")

print(f"\nAll done. Models and outputs saved under: {OUTPUT_DIR}\n")
