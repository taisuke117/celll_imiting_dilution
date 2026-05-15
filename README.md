# 🔬 Limiting Dilution Calculator

**何プレートをどのλで播くか → 必要な細胞数・懸濁液量を計算するPythonツール**

限界希釈（Limiting Dilution）実験で「何細胞を何mLに調製すればよいか」を  
プレート数・λ値から逆算します。

---

## 🚀 クイックスタート

### Google Colab（推奨）

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_USERNAME/limiting-dilution/blob/main/limiting_dilution.ipynb)

1. バッジをクリック → Colab で開く  
2. `Ctrl+F9`（すべてのセルを実行）  
3. フォームに条件を入力 → **「計算する」**

### ローカル実行

```bash
pip install scipy matplotlib numpy
python limiting_dilution.py
```

### モジュールとして使用

```python
from limiting_dilution import calculate, render_html, plot_result

# 3プレート、λ=0.6、96-well、100 µL/well
r = calculate(lambda_val=0.6, well_vol_uL=100, n_wells=96, n_plates=3)

print(f"必要細胞数: {r['total_cells']:.2e} cells")
print(f"懸濁液量:   {r['total_vol_mL']:.2f} mL")
print(f"目標濃度:   {r['target_conc']:.1f} cells/mL")

plot_result(r)   # Poisson グラフ + λ 比較表
```

---

## 📊 計算内容

| 計算項目 | 式 |
|---------|---|
| 目標濃度 (cells/mL) | λ ÷ ウェル容量 (mL) |
| 必要懸濁液量 (mL) | プレート数 × ウェル数 × ウェル容量 × 安全係数 |
| 必要細胞数 | 目標濃度 × 必要懸濁液量 |

### 出力例（3 plates × 96 wells × 100 µL、λ = 0.6）

```
必要細胞数   :  2.07 × 10² cells
懸濁液量     :  34.56 mL
目標濃度     :  6.0 cells/mL
P(1) 単細胞率:  32.9%
予想シングル :  ~95 / 288 wells
```

> **手順イメージ:** カウントした細胞を `6.0 cells/mL` になるよう培地で希釈 →  
> 各ウェルに `100 µL` ずつ分注 → 完了

---

## ⚙️ パラメータ

| パラメータ | デフォルト | 説明 |
|-----------|-----------|------|
| `lambda_val` | `0.6` | Poisson λ = cells/well |
| `well_vol_uL` | `100` | ウェル容量 (µL) |
| `n_wells` | `96` | ウェル数/プレート |
| `n_plates` | `1` | プレート数 |
| `safety_factor` | `1.2` | 余裕係数（20% 余分） |

---

## 📐 理論背景

$$P(k; \lambda) = \frac{\lambda^k e^{-\lambda}}{k!}$$

| λ | P(0) | **P(1)** | P(≥2) |
|---|------|----------|-------|
| 0.3 | 74.1% | **22.2%** | 3.7% |
| 0.5 | 60.7% | **30.3%** | 9.0% |
| **0.6** | **54.9%** | **32.9%** | **12.2%** ← 一般推奨 |
| 1.0 | 36.8% | **36.8%** | 26.4% |

> Lefkovits & Waldmann (1984) *Limiting Dilution Analysis of Cells of the Immune System.* Cambridge University Press.

---

## 🗂 ファイル構成

```
limiting-dilution/
├── limiting_dilution.py      # メインモジュール（計算 + 可視化 + CLI）
├── limiting_dilution.ipynb   # Google Colab ノートブック
└── README.md
```

## 📋 Requirements

- Python ≥ 3.10
- `numpy`, `scipy`, `matplotlib`
- `ipywidgets` (Colab / Jupyter のみ)

## 📄 License

MIT
