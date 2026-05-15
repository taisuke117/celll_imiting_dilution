# 🔬 Limiting Dilution Calculator

**細胞の限界希釈（Limiting Dilution）プロトコールを自動計算するPythonツール**

手持ちの細胞濃度とボリュームを入力するだけで、段階希釈のステップ・各チューブへの分注量・ポアソン統計を自動計算します。

---

## ✨ 機能

- **段階希釈プロトコールの自動計算** — 総希釈倍率を実用的なステップに分解（例: 1:10 × 5回 + 最終調整）
- **ボリューム計算** — 各チューブの分注量・培地量・合計量を表示、最低分注量も考慮
- **ポアソン統計** — P(0), P(1), P(≥2) の自動計算と棒グラフ表示
- **ストック不足の警告** — 入力ボリュームが足りない場合にアラート
- **Colab widget UI** — ipywidgets による対話型フォーム
- **CLI モード** — ローカルPythonからの実行にも対応

---

## 🚀 クイックスタート

### Google Colab（推奨）

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_USERNAME/limiting-dilution/blob/main/limiting_dilution.ipynb)

1. 上のバッジをクリック → Google Colab で開く
2. 「ランタイム → すべてのセルを実行」（`Ctrl+F9`）
3. フォームに値を入力して **「計算する」** をクリック

### ローカル実行

```bash
pip install scipy matplotlib numpy
python limiting_dilution.py
```

### モジュールとして使用

```python
from limiting_dilution import calculate_protocol, print_protocol, plot_result

result = calculate_protocol(
    input_conc   = 1e6,   # cells/mL
    input_vol_mL = 1.0,   # mL
    lambda_val   = 0.6,   # cells/well (Poisson λ)
    well_vol_uL  = 100,   # µL/well
    n_wells      = 96,
    max_step_df  = 10,    # max dilution per step (1:10)
)

print_protocol(result)
plot_result(result)
```

---

## 📊 出力例

```
====================================================================
  🔬  LIMITING DILUTION PROTOCOL
====================================================================
  Input : 1.00e+06 cells/mL  (1.0 mL available)
  Target: 6.00e+00 cells/mL  (λ = 0.6 / 100 µL well)
  Total dilution factor: 1/1.67e+05  →  6 step(s)

  Step  From (c/mL)    DF        Transfer    + Medium    = Total   To (c/mL)
  ──────────────────────────────────────────────────────────────────
  1     1.000e+06     1:10       30.0 µL     270.0 µL    300.0 µL  1.000e+05
  2     1.000e+05     1:10       30.0 µL     270.0 µL    300.0 µL  1.000e+04
  3     1.000e+04     1:10       30.0 µL     270.0 µL    300.0 µL  1.000e+03
  4     1.000e+03     1:10      126.5 µL    1138.9 µL   1265.5 µL  1.000e+02
  5     1.000e+02     1:10      973.4 µL    8761.0 µL   9734.4 µL  1.000e+01
  6     1.000e+01     1:2      7488.0 µL    4992.0 µL  12480.0 µL  6.000e+00

  PLATING: Take 100 µL from Tube 6 → each of 96 wells

  POISSON STATISTICS
  P(0 cells/well) = 54.9%   ← empty wells
  P(1 cell/well)  = 32.9%   ← single-cell cloning ✓
  P(≥2 cells/well)= 12.2%   ← multi-cell contamination
  Expected single-cell wells: ~32 / 96
```

---

## ⚙️ パラメータ

| パラメータ | デフォルト | 説明 |
|-----------|-----------|------|
| `input_conc` | — | 手持ち細胞濃度 (cells/mL) |
| `input_vol_mL` | — | 手持ちボリューム (mL) |
| `lambda_val` | `0.6` | 目標 cells/well (Poisson λ) |
| `well_vol_uL` | `100` | ウェル容量 (µL) |
| `n_wells` | `96` | プレートのウェル数 |
| `max_step_df` | `10` | 1ステップあたりの最大希釈倍率 |
| `min_transfer_uL` | `20` | 最小分注量 (µL) |
| `safety_factor` | `1.3` | 余裕係数（30%余分に確保）|

---

## 📐 理論背景

Limiting dilution（限界希釈）では、各ウェルに入る細胞数が **ポアソン分布** に従うことを利用します。

$$P(k; \lambda) = \frac{\lambda^k e^{-\lambda}}{k!}$$

| λ | P(0) | **P(1)** | P(≥2) | 用途 |
|---|------|----------|-------|------|
| 0.3 | 74.1% | **22.2%** | 3.7% | 厳格なシングルセル確認 |
| 0.5 | 60.7% | **30.3%** | 9.0% | 保守的 |
| **0.6** | **54.9%** | **32.9%** | **12.2%** | **一般推奨** |
| 1.0 | 36.8% | **36.8%** | 26.4% | 効率優先 |

> **参考:** Lefkovits & Waldmann (1984) *Limiting Dilution Analysis of Cells of the Immune System.* Cambridge University Press.

---

## 🗂 ファイル構成

```
limiting-dilution/
├── limiting_dilution.py      # メインモジュール（計算 + 可視化 + CLI）
├── limiting_dilution.ipynb   # Google Colab ノートブック
└── README.md
```

---

## 📋 Requirements

- Python ≥ 3.8
- `numpy`, `scipy`, `matplotlib`
- `ipywidgets` (Colab / Jupyter のみ)

---

## 📄 License

MIT License
