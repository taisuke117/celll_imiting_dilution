"""
limiting_dilution.py
====================
Limiting Dilution — Cell Preparation Calculator

「何プレートをλいくつで播くのに、何細胞・何mL用意するか」を計算するツール。

Usage (Google Colab / Jupyter):
    %run limiting_dilution.py   # または
    from limiting_dilution import launch_ui; launch_ui()

Usage (standalone):
    python limiting_dilution.py

Reference:
    Lefkovits & Waldmann (1984) Limiting Dilution Analysis of Cells
    of the Immune System. Cambridge University Press.
"""

import math
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import poisson

warnings.filterwarnings("ignore")

try:
    from IPython.display import display, HTML, clear_output
    import ipywidgets as widgets
    _COLAB = True
except ImportError:
    _COLAB = False


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. Core calculation
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def calculate(
    lambda_val: float = 0.6,
    well_vol_uL: float = 100.0,
    n_wells: int = 96,
    n_plates: int = 1,
    safety_factor: float = 1.2,
) -> dict:
    """
    Calculate required cell number and suspension volume for limiting dilution.

    Parameters
    ----------
    lambda_val    : Poisson λ = target cells per well
    well_vol_uL   : volume per well (µL)
    n_wells       : wells per plate (96 or 384)
    n_plates      : number of plates to seed
    safety_factor : extra buffer (default 1.2 = 20% extra)

    Returns
    -------
    dict with all derived quantities
    """
    well_vol_mL   = well_vol_uL / 1000.0
    target_conc   = lambda_val / well_vol_mL          # cells/mL — add this to each well

    total_wells   = n_wells * n_plates
    total_vol_mL  = total_wells * well_vol_mL * safety_factor
    total_cells   = target_conc * total_vol_mL

    po = poisson(lambda_val)
    p0     = po.pmf(0)
    p1     = po.pmf(1)
    p2plus = 1.0 - po.cdf(1)

    expected_empty  = total_wells * p0
    expected_single = total_wells * p1
    expected_multi  = total_wells * p2plus

    return {
        "lambda_val":       lambda_val,
        "well_vol_uL":      well_vol_uL,
        "well_vol_mL":      well_vol_mL,
        "n_wells":          n_wells,
        "n_plates":         n_plates,
        "safety_factor":    safety_factor,
        "target_conc":      target_conc,        # cells/mL
        "total_wells":      total_wells,
        "total_vol_mL":     total_vol_mL,       # mL of suspension to prepare
        "total_cells":      total_cells,        # cells to prepare
        "p0":               p0,
        "p1":               p1,
        "p2plus":           p2plus,
        "expected_empty":   expected_empty,
        "expected_single":  expected_single,
        "expected_multi":   expected_multi,
    }


def table_by_lambda(
    well_vol_uL: float = 100.0,
    n_wells: int = 96,
    n_plates: int = 1,
    safety_factor: float = 1.2,
    lambdas: tuple = (0.3, 0.5, 0.6, 1.0),
) -> list[dict]:
    """Return calculate() results for multiple λ values."""
    return [calculate(lam, well_vol_uL, n_wells, n_plates, safety_factor)
            for lam in lambdas]


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. Visualization
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

_LAMBDA_COLORS = {
    0.3: "#2980b9",
    0.5: "#27ae60",
    0.6: "#e67e22",
    1.0: "#e74c3c",
}
_DEFAULT_COLOR = "#8e44ad"

RECOMMENDED_LAMBDA = 0.6


def _bar_color(lam):
    return _LAMBDA_COLORS.get(lam, _DEFAULT_COLOR)


def plot_result(result: dict, fig=None, axes=None):
    """
    Two-panel figure:
      Left  — Poisson PMF bar chart
      Right — λ comparison table (cells & volume for λ = 0.3, 0.5, 0.6, 1.0)
    """
    if fig is None:
        fig = plt.figure(figsize=(13, 5.5), facecolor="white")

    gs = gridspec.GridSpec(1, 2, figure=fig,
                           width_ratios=[1, 1.4], wspace=0.10)
    ax_poi = fig.add_subplot(gs[0])
    ax_tbl = fig.add_subplot(gs[1])

    _plot_poisson(ax_poi, result)
    _plot_comparison_table(ax_tbl, result)

    fig.suptitle(
        f"Limiting Dilution  —  {result['n_plates']} plate{'s' if result['n_plates'] > 1 else ''} "
        f"× {result['n_wells']} wells, {result['well_vol_uL']:.0f} µL/well",
        fontsize=12, fontweight="bold", y=1.01,
    )
    plt.tight_layout(pad=2.2)
    plt.show()


def _plot_poisson(ax, result):
    lam    = result["lambda_val"]
    k_vals = list(range(7))
    probs  = [poisson.pmf(k, lam) * 100 for k in k_vals]
    probs[6] = (1 - poisson.cdf(5, lam)) * 100
    xlabels = ["0", "1", "2", "3", "4", "5", "≥6"]

    base_col = _bar_color(lam)
    colors = [
        "#95a5a6",    # 0  — empty
        "#27ae60",    # 1  — target ✓
        "#e67e22",    # 2
        "#e74c3c",    # 3
        "#c0392b",    # 4
        "#922b21",    # 5
        "#7b241c",    # ≥6
    ]
    bars = ax.bar(range(7), probs, color=colors,
                  edgecolor="white", linewidth=1.4, width=0.62, zorder=2)

    for bar, p in zip(bars, probs):
        if p > 0.5:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5,
                    f"{p:.1f}%", ha="center", va="bottom",
                    fontsize=9, fontweight="bold")

    ax.set_xticks(range(7))
    ax.set_xticklabels(xlabels, fontsize=11)
    ax.set_xlabel("Cells per well (k)", fontsize=11)
    ax.set_ylabel("Probability (%)", fontsize=11)
    ax.set_title(f"Poisson  (λ = {lam})", fontsize=12, fontweight="bold", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_facecolor("white")
    ax.set_ylim(0, max(probs) * 1.30)
    ax.yaxis.grid(True, linestyle="--", alpha=0.35, zorder=0)

    tw = result["total_wells"]
    info = (
        f"P(0) = {result['p0']:.1%}  →  ~{result['expected_empty']:.0f} empty wells\n"
        f"P(1) = {result['p1']:.1%}  →  ~{result['expected_single']:.0f} single-cell ✓\n"
        f"P(≥2) = {result['p2plus']:.1%}  →  ~{result['expected_multi']:.0f} multi-cell"
    )
    ax.text(0.97, 0.97, info, transform=ax.transAxes,
            ha="right", va="top", fontsize=8.8, fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#eaf4fb",
                      edgecolor="#a9cce3", alpha=0.92))


def _plot_comparison_table(ax, result):
    """Right panel: λ comparison with highlighted selected row."""
    ax.axis("off")
    ax.set_facecolor("white")
    ax.set_title("Required cells & volume by λ", fontsize=12, fontweight="bold", pad=8)

    lambdas      = [0.3, 0.5, 0.6, 1.0]
    rows_data    = table_by_lambda(
        well_vol_uL   = result["well_vol_uL"],
        n_wells       = result["n_wells"],
        n_plates      = result["n_plates"],
        safety_factor = result["safety_factor"],
        lambdas       = lambdas,
    )

    col_labels = ["λ", "Target conc.\n(cells/mL)", "Volume\n(mL)", "Cells needed", "P(1)\nsingle-cell",
                  "Expected\nsingle-cell wells"]
    col_widths = [0.08, 0.20, 0.14, 0.22, 0.14, 0.18]

    # Header
    y_header = 0.93
    x_start  = 0.02
    xs = [x_start + sum(col_widths[:i]) for i in range(len(col_widths))]

    ax.add_patch(plt.Rectangle((x_start - 0.01, y_header - 0.045),
                                sum(col_widths) + 0.02, 0.055,
                                transform=ax.transAxes, zorder=2,
                                facecolor="#1a4f7a", edgecolor="none"))
    for x, label in zip(xs, col_labels):
        ax.text(x + col_widths[col_labels.index(label)] / 2,
                y_header - 0.015,
                label, transform=ax.transAxes,
                ha="center", va="center", fontsize=8.5,
                color="white", fontweight="bold")

    # Data rows
    row_h  = 0.115
    y0     = y_header - 0.045

    for i, (lam, row) in enumerate(zip(lambdas, rows_data)):
        y_top  = y0 - i * row_h
        y_mid  = y_top - row_h / 2
        is_sel = abs(lam - result["lambda_val"]) < 1e-9

        bg = "#fff8ee" if is_sel else ("#f5faff" if i % 2 == 0 else "white")
        edge_col = "#e67e22" if is_sel else "none"
        lw = 2.0 if is_sel else 0

        ax.add_patch(plt.Rectangle((x_start - 0.01, y_top - row_h),
                                    sum(col_widths) + 0.02, row_h,
                                    transform=ax.transAxes, zorder=1,
                                    facecolor=bg, edgecolor=edge_col, linewidth=lw))

        cells = row["total_cells"]
        cells_str = (f"{cells/1e6:.2f}×10⁶" if cells >= 1e6
                     else f"{cells/1e3:.1f}×10³" if cells >= 1e3
                     else f"{cells:.0f}")

        conc = row["target_conc"]
        conc_str = (f"{conc:.0f}" if conc >= 1
                    else f"{conc:.2f}")

        values = [
            f"{'★ ' if is_sel else ''}{lam}",
            conc_str,
            f"{row['total_vol_mL']:.2f}",
            cells_str,
            f"{row['p1']:.1%}",
            f"~{row['expected_single']:.0f}",
        ]

        fw = "bold" if is_sel else "normal"
        fs = 9.2 if is_sel else 8.8

        for j, (x, val) in enumerate(zip(xs, values)):
            ax.text(x + col_widths[j] / 2, y_mid, val,
                    transform=ax.transAxes,
                    ha="center", va="center",
                    fontsize=fs, fontweight=fw,
                    color="#1a1a1a")

    # Footer note
    sf_pct = int((result["safety_factor"] - 1) * 100)
    note = (f"★ = selected lambda   /   +{sf_pct}% safety buffer   "
            f"/   {result['well_vol_uL']:.0f} uL per well")
    ax.text(0.5, 0.01, note, transform=ax.transAxes,
            ha="center", va="bottom", fontsize=8, color="#666",
            style="italic")


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. HTML table (Colab)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

_CSS = """
<style>
.ld { font-family: 'Helvetica Neue', Arial, sans-serif; max-width: 820px; }
.ld-card { background:#f8fafd; border:1.5px solid #c8dff5; border-radius:10px;
           padding:16px 22px; margin:10px 0; }
.ld-h1   { font-size:1.2em; font-weight:700; color:#1a4f7a; margin:0 0 4px 0; }
.ld-sub  { color:#5d6d7e; font-size:0.88em; margin:0 0 14px 0; }
.ld-big  { font-size:2em; font-weight:800; color:#1a4f7a; letter-spacing:-0.5px; }
.ld-unit { font-size:0.85em; color:#5d6d7e; margin-left:4px; }
.ld-kpi  { display:inline-block; background:#eaf4fb; border:1px solid #aed6f1;
           border-radius:8px; padding:10px 18px; margin:4px 6px; text-align:center; }
.ld-tbl  { border-collapse:collapse; width:100%; font-size:0.9em; margin-top:12px; }
.ld-tbl th { background:#1a4f7a; color:#fff; padding:8px 12px; text-align:center; }
.ld-tbl td { padding:7px 12px; text-align:center; border-bottom:1px solid #e0eaf5; }
.ld-tbl tr:nth-child(even) td { background:#f0f8ff; }
.ld-tbl tr.sel td { background:#fff3e0; font-weight:700;
                    border-top:2px solid #e67e22; border-bottom:2px solid #e67e22; }
.ld-tbl td:first-child { font-weight:700; }
.ld-badge { display:inline-block; background:#1a4f7a; color:#fff; border-radius:6px;
            padding:3px 11px; margin:2px 3px; font-size:0.88em; }
.ld-note  { color:#888; font-size:0.82em; margin-top:10px; }
</style>
"""


def render_html(result: dict) -> str:
    """Styled HTML summary card + λ comparison table."""
    lam  = result["lambda_val"]
    rows = table_by_lambda(
        well_vol_uL   = result["well_vol_uL"],
        n_wells       = result["n_wells"],
        n_plates      = result["n_plates"],
        safety_factor = result["safety_factor"],
    )

    def _cells_str(c):
        if c >= 1e6:   return f"{c/1e6:.3f} × 10<sup>6</sup>"
        if c >= 1e3:   return f"{c/1e3:.2f} × 10<sup>3</sup>"
        return f"{c:.0f}"

    # KPI cards for selected λ
    r = result
    kpi_cells = _cells_str(r["total_cells"])
    kpi_vol   = f"{r['total_vol_mL']:.2f}"
    kpi_conc  = f"{r['target_conc']:.1f}" if r["target_conc"] >= 1 else f"{r['target_conc']:.3f}"

    kpi_html = f"""
<div class="ld-kpi">
  <div style="font-size:0.8em;color:#888;margin-bottom:2px">必要細胞数</div>
  <span class="ld-big">{kpi_cells}</span>
  <span class="ld-unit">cells</span>
</div>
<div class="ld-kpi">
  <div style="font-size:0.8em;color:#888;margin-bottom:2px">懸濁液量</div>
  <span class="ld-big">{kpi_vol}</span>
  <span class="ld-unit">mL</span>
</div>
<div class="ld-kpi">
  <div style="font-size:0.8em;color:#888;margin-bottom:2px">懸濁液濃度</div>
  <span class="ld-big">{kpi_conc}</span>
  <span class="ld-unit">cells/mL</span>
</div>
"""

    # Comparison table rows
    tbl_rows = ""
    for row in rows:
        is_sel = abs(row["lambda_val"] - lam) < 1e-9
        cls    = ' class="sel"' if is_sel else ""
        star   = "★ " if is_sel else ""
        tbl_rows += (
            f"<tr{cls}>"
            f"<td>{star}{row['lambda_val']}</td>"
            f"<td>{row['target_conc']:.1f}" + (" (= λ / vol)" if is_sel else "") + "</td>"
            f"<td>{row['total_vol_mL']:.2f}</td>"
            f"<td>{_cells_str(row['total_cells'])}</td>"
            f"<td>{row['p1']:.1%}</td>"
            f"<td>~{row['expected_single']:.0f} / {result['total_wells']}</td>"
            f"</tr>"
        )

    sf_pct = int((result["safety_factor"] - 1) * 100)

    return f"""
{_CSS}
<div class="ld">
<div class="ld-card">
  <p class="ld-h1">🔬 Limiting Dilution — 準備量の計算</p>
  <p class="ld-sub">
    {result['n_plates']} plate{'s' if result['n_plates'] > 1 else ''} ×
    {result['n_wells']} wells,
    {result['well_vol_uL']:.0f} µL/well,
    λ = {lam}　（安全係数 +{sf_pct}%）
  </p>

  <div style="margin-bottom:14px">{kpi_html}</div>

  <table class="ld-tbl">
    <thead><tr>
      <th>λ</th>
      <th>目標濃度 (cells/mL)</th>
      <th>必要量 (mL)</th>
      <th>必要細胞数</th>
      <th>P(1) 単細胞率</th>
      <th>予想シングルwells</th>
    </tr></thead>
    <tbody>{tbl_rows}</tbody>
  </table>

  <p class="ld-note">
    ★ 選択中の λ ／
    必要量 = {result['n_plates']} plate × {result['n_wells']} wells ×
    {result['well_vol_uL']:.0f} µL × {result['safety_factor']} (安全係数) ／
    必要濃度 = λ ÷ {result['well_vol_uL']:.0f} µL = {kpi_conc} cells/mL
  </p>
</div>
</div>
"""


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. Widget UI (Colab / Jupyter)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def launch_ui():
    if not _COLAB:
        _cli(); return

    S = {"description_width": "180px"}
    L = widgets.Layout(width="400px")

    w_plates = widgets.BoundedIntText(
        value=1, min=1, max=50,
        description="プレート数:", style=S, layout=L)

    w_lambda = widgets.Dropdown(
        options=[("0.3  — very conservative", 0.3),
                 ("0.5  — conservative",       0.5),
                 ("0.6  — recommended ✓",      0.6),
                 ("1.0  — higher seeding",     1.0)],
        value=0.6,
        description="λ (cells/well):", style=S, layout=L)

    w_well_vol = widgets.Dropdown(
        options=[("100 µL  (96-well standard)", 100.0),
                 ("200 µL",                     200.0),
                 ("50 µL",                       50.0)],
        value=100.0,
        description="ウェル容量 (µL):", style=S, layout=L)

    w_n_wells = widgets.Dropdown(
        options=[96, 384, 48, 24], value=96,
        description="ウェル数:", style=S, layout=L)

    w_safety = widgets.FloatSlider(
        value=1.2, min=1.0, max=2.0, step=0.05,
        description="安全係数:",
        readout_format=".2f",
        style=S, layout=L)

    btn = widgets.Button(
        description="　計算する　",
        button_style="primary",
        icon="calculator",
        layout=widgets.Layout(width="190px", height="42px",
                              margin="14px 0 4px 180px"))

    out = widgets.Output()
    sep = widgets.HTML("<hr style='margin:8px 0;border-color:#c8dff5;'>")

    ui = widgets.VBox([
        widgets.HTML(
            "<h2 style='color:#1a4f7a;margin:4px 0 2px'>🔬 Limiting Dilution Calculator</h2>"
            "<p style='color:#5d6d7e;font-size:0.9em;margin:0 0 10px'>"
            "何プレートをどのλで播くか → 必要な細胞数・懸濁液量を計算</p>"),
        widgets.HTML("<b style='color:#2c3e50'>▸ 実験条件</b>"),
        w_plates,
        sep,
        widgets.HTML("<b style='color:#2c3e50'>▸ プレート設定</b>"),
        w_lambda, w_well_vol, w_n_wells,
        sep,
        widgets.HTML("<b style='color:#2c3e50'>▸ オプション</b>"),
        w_safety,
        btn,
        out,
    ])

    def on_click(_):
        with out:
            clear_output(wait=True)
            try:
                result = calculate(
                    lambda_val   = w_lambda.value,
                    well_vol_uL  = w_well_vol.value,
                    n_wells      = int(w_n_wells.value),
                    n_plates     = w_plates.value,
                    safety_factor= w_safety.value,
                )
                display(HTML(render_html(result)))
                plot_result(result)
            except Exception as e:
                display(HTML(
                    f"<div style='background:#fef9e7;border-left:4px solid #e74c3c;"
                    f"padding:10px;border-radius:4px'>❌ {e}</div>"))

    btn.on_click(on_click)
    display(ui)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. CLI mode
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _cli():
    print("=== Limiting Dilution Calculator ===\n")
    try:
        n_plates = int(input("プレート数 [default 1]: ") or 1)
        lam      = float(input("λ (cells/well) [default 0.6]: ") or 0.6)
        well_vol = float(input("ウェル容量 µL [default 100]: ") or 100)
        n_wells  = int(input("ウェル数 [default 96]: ") or 96)
        sf       = float(input("安全係数 [default 1.2]: ") or 1.2)
    except ValueError:
        print("入力値が無効です。"); return

    result = calculate(lam, well_vol, n_wells, n_plates, sf)
    r = result

    print(f"\n{'='*55}")
    print(f"  {n_plates} plate{'s' if n_plates > 1 else ''} × {n_wells} wells, "
          f"{well_vol:.0f} µL/well, λ = {lam}")
    print(f"{'='*55}")
    print(f"  目標濃度   : {r['target_conc']:.2f} cells/mL")
    print(f"  必要量     : {r['total_vol_mL']:.2f} mL")
    print(f"  必要細胞数 : {r['total_cells']:.2e} cells")
    print(f"\n  P(0) = {r['p0']:.1%}   P(1) = {r['p1']:.1%}   P(≥2) = {r['p2plus']:.1%}")
    print(f"  予想シングルセルwell: ~{r['expected_single']:.0f} / {r['total_wells']}")
    print(f"{'='*55}\n")

    print("λ 比較:")
    rows = table_by_lambda(well_vol, n_wells, n_plates, sf)
    print(f"  {'λ':>5}  {'濃度 (c/mL)':>12}  {'量 (mL)':>9}  {'細胞数':>14}  P(1)")
    print("  " + "-" * 52)
    for row in rows:
        cells = row["total_cells"]
        mark  = " ★" if abs(row["lambda_val"] - lam) < 1e-9 else ""
        print(f"  {row['lambda_val']:>5}  {row['target_conc']:>12.2f}  "
              f"{row['total_vol_mL']:>9.2f}  {cells:>14.2e}  {row['p1']:.1%}{mark}")

    plot_result(result)


if __name__ == "__main__":
    _cli()
