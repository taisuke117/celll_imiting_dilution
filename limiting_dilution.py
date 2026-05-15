"""
limiting_dilution.py
====================
Limiting Dilution (限界希釈) Protocol Calculator
for single-cell cloning experiments.

Usage (Google Colab / Jupyter):
    from limiting_dilution import launch_ui
    launch_ui()

Usage (standalone):
    python limiting_dilution.py

References:
    Lefkovits & Waldmann (1984) Limiting Dilution Analysis of Cells of
    the Immune System. Cambridge University Press.
"""

import math
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from scipy.stats import poisson

warnings.filterwarnings("ignore")

# Try to import Colab-specific modules
try:
    from IPython.display import display, HTML, clear_output
    import ipywidgets as widgets
    _COLAB_AVAILABLE = True
except ImportError:
    _COLAB_AVAILABLE = False


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 1. Core calculation
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def calculate_protocol(
    input_conc: float,
    input_vol_mL: float,
    lambda_val: float = 0.6,
    well_vol_uL: float = 100.0,
    n_wells: int = 96,
    max_step_df: int = 10,
    min_transfer_uL: float = 20.0,
    safety_factor: float = 1.3,
) -> dict:
    """
    Calculate a limiting dilution protocol.

    Parameters
    ----------
    input_conc      : float  – Starting cell concentration (cells/mL)
    input_vol_mL    : float  – Available volume of starting suspension (mL)
    lambda_val      : float  – Target Poisson lambda = cells per well (default 0.6)
    well_vol_uL     : float  – Volume per well in µL (default 100)
    n_wells         : int    – Number of wells to seed (default 96)
    max_step_df     : int    – Maximum dilution factor per step (default 10 = 1:10)
    min_transfer_uL : float  – Minimum pipetting volume in µL (default 20)
    safety_factor   : float  – Extra volume buffer factor (default 1.3 = 30% extra)

    Returns
    -------
    dict with protocol steps, Poisson statistics, and volume flags.
    """
    if input_conc <= 0:
        raise ValueError("細胞濃度には正の値を入力してください。")
    if input_vol_mL <= 0:
        raise ValueError("ボリュームには正の値を入力してください。")

    well_vol_mL = well_vol_uL / 1000.0
    target_conc = lambda_val / well_vol_mL  # cells/mL

    if input_conc < target_conc:
        raise ValueError(
            f"入力濃度 ({input_conc:.2e} cells/mL) が目標濃度 "
            f"({target_conc:.2e} cells/mL) より低いです。\n"
            "細胞を濃縮するか、λ値を下げてください。"
        )

    total_df = input_conc / target_conc

    # ── Factorize total dilution factor into serial steps ──────────
    step_dfs = []
    remaining = total_df
    while remaining > max_step_df * 1.001:
        step_dfs.append(float(max_step_df))
        remaining /= max_step_df
    if remaining > 1.001:
        step_dfs.append(remaining)
    if not step_dfs:
        step_dfs = [total_df]

    # ── Calculate tube volumes (backwards from final tube) ─────────
    plate_vol_uL  = n_wells * well_vol_uL * safety_factor  # volume needed for plating
    min_tube_vol  = 300.0                                   # µL, minimum per tube

    tube_vols = [max(plate_vol_uL, min_tube_vol)]

    for df in reversed(step_dfs[1:]):
        transfer_needed = tube_vols[0] / df
        # Enforce minimum transfer volume by scaling up downstream volumes
        if transfer_needed < min_transfer_uL:
            scale = min_transfer_uL / transfer_needed
            tube_vols = [v * scale for v in tube_vols]
            transfer_needed = min_transfer_uL
        this_vol = max(transfer_needed * safety_factor, min_tube_vol)
        tube_vols.insert(0, this_vol)

    # ── Build step list ────────────────────────────────────────────
    steps = []
    conc = input_conc
    for i, (df, tvol) in enumerate(zip(step_dfs, tube_vols)):
        transfer = tvol / df
        media    = tvol - transfer
        new_conc = conc / df
        steps.append({
            "step":        i + 1,
            "from_conc":   conc,
            "df":          df,
            "transfer_uL": round(transfer, 1),
            "media_uL":    round(media, 1),
            "total_uL":    round(tvol, 1),
            "to_conc":     new_conc,
        })
        conc = new_conc

    # ── Poisson statistics ─────────────────────────────────────────
    po     = poisson(lambda_val)
    p0     = po.pmf(0)
    p1     = po.pmf(1)
    p2plus = 1.0 - po.cdf(1)

    vol_needed_uL = steps[0]["transfer_uL"]
    volume_ok     = input_vol_mL * 1000 >= vol_needed_uL

    return {
        "input_conc":      input_conc,
        "input_vol_mL":    input_vol_mL,
        "target_conc":     target_conc,
        "lambda_val":      lambda_val,
        "well_vol_uL":     well_vol_uL,
        "n_wells":         n_wells,
        "total_df":        total_df,
        "n_steps":         len(steps),
        "steps":           steps,
        "final_conc":      conc,
        "plate_vol_uL":    plate_vol_uL,
        "vol_needed_uL":   vol_needed_uL,
        "volume_ok":       volume_ok,
        "p0":              p0,
        "p1":              p1,
        "p2plus":          p2plus,
        "expected_single": n_wells * p1,
    }


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 2. Visualization
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _fmt(c: float) -> str:
    """Format concentration as compact scientific notation string."""
    if c == 0:
        return "0"
    exp = int(math.floor(math.log10(abs(c))))
    m   = c / 10 ** exp
    if abs(m - round(m)) < 0.05:
        return f"{m:.0f}e{exp}"
    return f"{m:.2f}e{exp}"


def _box(ax, x, y, w, h, color, line1="", line2="", label=""):
    """Draw a rounded box (tube) on a matplotlib axes."""
    ax.add_patch(FancyBboxPatch(
        (x - w, y - h), 2 * w, 2 * h,
        boxstyle="round,pad=0.01",
        facecolor=color, edgecolor="white", linewidth=2, alpha=0.92,
        transform=ax.transData, zorder=2,
    ))
    ax.text(x, y + h * 0.18, line1, ha="center", va="center",
            fontsize=7.5, color="white", fontweight="bold", zorder=3, clip_on=True)
    ax.text(x, y - h * 0.28, line2, ha="center", va="center",
            fontsize=6.5, color="white", alpha=0.88, zorder=3, clip_on=True)
    if label:
        ax.text(x, y - h - 0.07, label, ha="center", va="top",
                fontsize=7.5, fontweight="bold", color="#2c3e50", zorder=3)


def _plate(ax, x, y, w, h, n_wells):
    """Draw a micro-plate icon."""
    ax.add_patch(FancyBboxPatch(
        (x - w, y - h), 2 * w, 2 * h,
        boxstyle="round,pad=0.01",
        facecolor="#ecf0f1", edgecolor="#95a5a6", linewidth=2, zorder=2,
    ))
    rows, cols = (8, 12) if n_wells == 96 else (16, 24) if n_wells == 384 else (4, 6)
    r_s, c_s = min(rows, 8), min(cols, 10)
    dxs = np.linspace(x - w * 0.80, x + w * 0.80, c_s)
    dys = np.linspace(y - h * 0.68, y + h * 0.68, r_s)
    ms  = 2.5 if n_wells <= 96 else 1.2
    for dx in dxs:
        for dy in dys:
            ax.plot(dx, dy, "o", color="#95a5a6", markersize=ms, zorder=3)
    ax.text(x, y - h - 0.07, f"{n_wells}-well plate",
            ha="center", va="top", fontsize=8.5, fontweight="bold", color="#2c3e50")


def _plot_poisson(ax, result):
    """Bar chart of Poisson PMF."""
    lam    = result["lambda_val"]
    k_vals = list(range(6))
    probs  = [poisson.pmf(k, lam) * 100 for k in k_vals]
    probs[5] = (1 - poisson.cdf(4, lam)) * 100
    labels = ["0", "1", "2", "3", "4", "≥5"]
    colors = ["#95a5a6", "#27ae60", "#e67e22", "#e74c3c", "#c0392b", "#7b241c"]

    bars = ax.bar(range(6), probs, color=colors, edgecolor="white",
                  linewidth=1.5, width=0.62, zorder=2)
    ax.set_xticks(range(6))
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_xlabel("Cells per well (k)", fontsize=11)
    ax.set_ylabel("Probability (%)", fontsize=11)
    ax.set_title(f"Poisson Distribution  (λ = {lam})",
                 fontsize=12, fontweight="bold", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_facecolor("white")
    ax.set_ylim(0, max(probs) * 1.28)
    ax.yaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)

    for bar, p in zip(bars, probs):
        if p > 0.5:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.5,
                    f"{p:.1f}%", ha="center", va="bottom",
                    fontsize=9, fontweight="bold")

    info = (f"P(0) = {result['p0']:.1%}   empty\n"
            f"P(1) = {result['p1']:.1%}   single ✓\n"
            f"P(≥2) = {result['p2plus']:.1%}  multi")
    ax.text(0.97, 0.97, info, transform=ax.transAxes,
            ha="right", va="top", fontsize=9, fontfamily="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#eaf4fb",
                      edgecolor="#a9cce3", alpha=0.9))


def _plot_scheme(ax, result):
    """Dilution flow diagram: Stock → Tube1 → … → TubeN → Plate."""
    ax.axis("off")
    ax.set_facecolor("white")
    ax.set_title("Dilution Scheme", fontsize=12, fontweight="bold", pad=8)

    steps = result["steps"]
    n_obj = len(steps) + 2  # Stock + tubes + plate

    xs  = np.linspace(0.06, 0.94, n_obj)
    y   = 0.60
    bw  = min(0.072, 0.78 / (2 * n_obj))
    bh  = 0.20

    COLORS = (["#2e86de"]
              + ["#8e44ad"] * max(len(steps) - 1, 0)
              + ["#e74c3c"])

    # Stock
    _box(ax, xs[0], y, bw, bh, "#2e86de",
         line1=_fmt(result["input_conc"]),
         line2="cells/mL", label="Stock")

    for i, step in enumerate(steps):
        x0  = xs[i] + bw
        x1  = xs[i + 1] - bw
        mid = (x0 + x1) / 2

        ax.annotate("", xy=(x1, y), xytext=(x0, y),
                    arrowprops=dict(arrowstyle="-|>", color="#444",
                                   lw=1.8, mutation_scale=14), zorder=4)

        n_s = len(steps)
        fsa = max(6.5, 8.5 - 0.4 * max(0, n_s - 3))
        ax.text(mid, y + bh + 0.09, f"{step['transfer_uL']:.1f} µL",
                ha="center", va="bottom", fontsize=fsa,
                color="#e67e22", fontweight="bold")
        ax.text(mid, y + bh + 0.01, f"(1:{step['df']:.0f})",
                ha="center", va="bottom", fontsize=max(6.0, fsa-1), color="#666")
        ax.text(mid, y - bh - 0.03, f"+ {step['media_uL']:.1f} µL",
                ha="center", va="top", fontsize=max(6.0, fsa-0.5), color="#27ae60")

        col = COLORS[min(i + 1, len(COLORS) - 1)]
        is_last = (i == len(steps) - 1)
        _box(ax, xs[i + 1], y, bw, bh, col,
             line1=_fmt(step["to_conc"]),
             line2="cells/mL",
             label=f"Tube {i+1}" + (" (final)" if is_last else ""))

    # Arrow to plate
    x0  = xs[-2] + bw
    x1  = xs[-1] - bw
    mid = (x0 + x1) / 2
    ax.annotate("", xy=(x1, y), xytext=(x0, y),
                arrowprops=dict(arrowstyle="-|>", color="#444",
                                lw=1.8, mutation_scale=14), zorder=4)
    ax.text(mid, y + bh + 0.09, f"{result['well_vol_uL']:.0f} µL/well",
            ha="center", va="bottom", fontsize=8.5,
            color="#e67e22", fontweight="bold")
    _plate(ax, xs[-1], y, bw, bh, result["n_wells"])

    # Footer summary
    foot = (f"Total dilution: 1/{result['total_df']:.2e}    "
            f"Target: {result['target_conc']:.2e} cells/mL    "
            f"Expected single-cell wells: ~{result['expected_single']:.0f} / {result['n_wells']}")
    ax.text(0.5, 0.06, foot, transform=ax.transAxes,
            ha="center", va="center", fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#eaf4fb",
                      edgecolor="#aed6f1", alpha=0.9))


def plot_result(result, figsize=None):
    """Main visualization: Poisson chart + dilution scheme side by side."""
    n = result["n_steps"]
    if figsize is None:
        figsize = (max(13, 5 + n * 2), 6)
    ratio = max(1.5, 0.65 * (n + 2))

    fig = plt.figure(figsize=figsize, facecolor="white")
    gs  = fig.add_gridspec(1, 2, width_ratios=[1, ratio], wspace=0.06)

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    _plot_poisson(ax1, result)
    _plot_scheme(ax2, result)

    plt.tight_layout(pad=2.5)
    plt.show()


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 3. Text / CLI output
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def print_protocol(result):
    """Print a human-readable protocol to stdout."""
    W = 68
    print("=" * W)
    print("  🔬  LIMITING DILUTION PROTOCOL")
    print("=" * W)
    print(f"  Input : {result['input_conc']:.2e} cells/mL  "
          f"({result['input_vol_mL']} mL available)")
    print(f"  Target: {result['target_conc']:.2e} cells/mL  "
          f"(λ = {result['lambda_val']} / {result['well_vol_uL']} µL well)")
    print(f"  Total dilution factor: 1/{result['total_df']:.2e}  "
          f"→  {result['n_steps']} step(s)")

    if not result["volume_ok"]:
        need = result["vol_needed_uL"] / 1000
        print(f"\n  ⚠️  WARNING: input volume ({result['input_vol_mL']} mL) "
              f"may be insufficient (need ≥ {need:.3f} mL)")

    print()
    print("─" * W)
    print(f"  {'Step':<5} {'From (c/mL)':<14} {'DF':<8} "
          f"{'Transfer':>10} {'+ Medium':>10} {'= Total':>10} {'To (c/mL)':<14}")
    print("  " + "─" * 64)

    for s in result["steps"]:
        print(f"  {s['step']:<5} {s['from_conc']:<14.3e} 1:{s['df']:<6.0f} "
              f"{s['transfer_uL']:>8.1f} µL {s['media_uL']:>8.1f} µL "
              f"{s['total_uL']:>8.1f} µL {s['to_conc']:<14.3e}")

    n = result["n_steps"]
    print("─" * W)
    print(f"\n  PLATING: Take {result['well_vol_uL']:.0f} µL from Tube {n} "
          f"→ each of {result['n_wells']} wells")
    print(f"  Volume needed: {result['plate_vol_uL']:.0f} µL "
          f"({result['plate_vol_uL']/1000:.2f} mL)\n")

    print("─" * W)
    print("  POISSON STATISTICS")
    print("─" * W)
    print(f"  P(0 cells/well) = {result['p0']:.1%}   ← empty wells")
    print(f"  P(1 cell/well)  = {result['p1']:.1%}   ← single-cell cloning ✓")
    print(f"  P(≥2 cells/well)= {result['p2plus']:.1%}   ← multi-cell contamination")
    print(f"\n  Expected single-cell wells: ~{result['expected_single']:.0f} / {result['n_wells']}")
    print("=" * W)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 4. HTML output (Colab)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

_CSS = """
<style>
.ld-wrap  { font-family: 'Helvetica Neue', Arial, sans-serif; max-width: 900px; }
.ld-card  { background: #f8fafd; border: 1px solid #c8dff5; border-radius: 10px;
            padding: 16px 20px; margin: 8px 0; }
.ld-h1    { font-size: 1.25em; font-weight: 700; color: #1a4f7a; margin: 0 0 4px 0; }
.ld-sub   { color: #5d6d7e; font-size: 0.88em; margin: 0 0 12px 0; }
.ld-tbl   { border-collapse: collapse; width: 100%; font-size: 0.9em; margin-top: 8px; }
.ld-tbl th { background: #2471a3; color: #fff; padding: 7px 10px; text-align: center; }
.ld-tbl td { padding: 6px 10px; text-align: center; border-bottom: 1px solid #e0eaf5; }
.ld-tbl tr:nth-child(even) td { background: #eaf4fb; }
.ld-tbl tr.final td { background: #fdebd0; font-weight: 600; }
.ld-warn  { background: #fef9e7; border-left: 4px solid #f39c12; padding: 8px 12px;
            border-radius: 4px; margin: 8px 0; font-size: 0.9em; }
.ld-ok    { background: #eafaf1; border-left: 4px solid #27ae60; padding: 8px 12px;
            border-radius: 4px; margin: 8px 0; font-size: 0.9em; }
.ld-badge { display: inline-block; background: #1a4f7a; color: #fff; border-radius: 6px;
            padding: 3px 11px; margin: 2px 3px; font-size: 0.88em; }
</style>
"""

def render_html(result) -> str:
    """Return styled HTML string of the protocol."""
    steps = result["steps"]
    rows  = ""
    for s in steps:
        cls = ' class="final"' if s["step"] == len(steps) else ""
        rows += (f"<tr{cls}>"
                 f"<td>{s['step']}</td>"
                 f"<td>{s['from_conc']:.3e}</td>"
                 f"<td>1 : {s['df']:.0f}</td>"
                 f"<td><b>{s['transfer_uL']:.1f}</b></td>"
                 f"<td>{s['media_uL']:.1f}</td>"
                 f"<td>{s['total_uL']:.1f}</td>"
                 f"<td>{s['to_conc']:.3e}</td>"
                 f"</tr>")

    warn = ""
    if not result["volume_ok"]:
        need = result["vol_needed_uL"] / 1000
        warn = (f'<div class="ld-warn">⚠️ 入力ボリューム ({result["input_vol_mL"]:.2f} mL) が不足する可能性があります。'
                f'最低 <b>{need:.3f} mL</b> 必要です。</div>')

    plating = (f'<div class="ld-ok">'
               f'✅ <b>プレーティング:</b> Tube {len(steps)} から '
               f'<b>{result["well_vol_uL"]:.0f} µL</b> を各ウェルに分注 '
               f'({result["n_wells"]} wells)  ／  '
               f'必要量: {result["plate_vol_uL"]:.0f} µL ({result["plate_vol_uL"]/1000:.2f} mL)'
               f'</div>')

    badges = "".join([
        f'<span class="ld-badge">λ = {result["lambda_val"]}</span>',
        f'<span class="ld-badge">P(0) = {result["p0"]:.1%}</span>',
        f'<span class="ld-badge">P(1) = {result["p1"]:.1%} ✓</span>',
        f'<span class="ld-badge">P(≥2) = {result["p2plus"]:.1%}</span>',
        f'<span class="ld-badge">予想シングルwell: ~{result["expected_single"]:.0f} / {result["n_wells"]}</span>',
    ])

    return f"""
{_CSS}
<div class="ld-wrap">
<div class="ld-card">
  <p class="ld-h1">🔬 Limiting Dilution Protocol</p>
  <p class="ld-sub">
    Input: {result["input_conc"]:.2e} cells/mL ({result["input_vol_mL"]} mL) &nbsp;|&nbsp;
    Target: {result["target_conc"]:.2e} cells/mL &nbsp;|&nbsp;
    総希釈倍率: 1/{result["total_df"]:.2e} &nbsp;({result["n_steps"]} step{'s' if result["n_steps"] > 1 else ''})
  </p>

  {warn}

  <table class="ld-tbl">
    <thead>
      <tr>
        <th>Step</th><th>From (cells/mL)</th><th>希釈倍率</th>
        <th>分注量 (µL)</th><th>培地量 (µL)</th><th>合計 (µL)</th><th>To (cells/mL)</th>
      </tr>
    </thead>
    <tbody>{rows}</tbody>
  </table>

  {plating}
  <div style="margin-top:10px">{badges}</div>
</div>
</div>
"""


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 5. Colab Widget UI
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def launch_ui():
    """Launch interactive widget UI (requires Google Colab / Jupyter)."""
    if not _COLAB_AVAILABLE:
        print("ipywidgets / IPython not available. Using CLI mode.\n")
        _cli_mode()
        return

    style  = {"description_width": "200px"}
    layout = widgets.Layout(width="420px")

    # ── Input widgets ──────────────────────────────────────────────
    w_conc = widgets.FloatText(
        value=1e6, description="細胞濃度 (cells/mL):",
        style=style, layout=layout)
    w_vol = widgets.FloatText(
        value=1.0, description="手持ちボリューム (mL):",
        style=style, layout=layout)

    w_lambda = widgets.Dropdown(
        options=[("0.3  — very conservative", 0.3),
                 ("0.5  — conservative",       0.5),
                 ("0.6  — recommended ✓",      0.6),
                 ("1.0  — higher seeding",     1.0)],
        value=0.6, description="λ (cells/well):",
        style=style, layout=layout)

    w_well_vol = widgets.Dropdown(
        options=[("100 µL  (96-well standard)", 100),
                 ("200 µL",                     200),
                 ("50 µL",                       50)],
        value=100, description="ウェル容量 (µL):",
        style=style, layout=layout)

    w_n_wells = widgets.Dropdown(
        options=[96, 384, 48, 24], value=96,
        description="ウェル数:", style=style, layout=layout)

    w_max_df = widgets.Dropdown(
        options=[("1:10 ずつ — 推奨",  10),
                 ("1:5  ずつ",          5),
                 ("1:20 ずつ",         20),
                 ("1:50 ずつ",         50),
                 ("1:100 ずつ",       100)],
        value=10, description="1ステップ最大希釈倍率:",
        style=style, layout=layout)

    w_min_tr = widgets.FloatText(
        value=20.0, description="最小分注量 (µL):",
        style=style, layout=layout)

    btn = widgets.Button(
        description=" 計算する ",
        button_style="primary",
        icon="calculator",
        layout=widgets.Layout(width="200px", height="42px", margin="14px 0 4px 200px"))

    out = widgets.Output()
    sep = widgets.HTML("<hr style='margin:10px 0; border-color:#c8dff5;'>")

    ui = widgets.VBox([
        widgets.HTML(
            "<h2 style='color:#1a4f7a;margin:4px 0 2px 0'>🔬 Limiting Dilution Calculator</h2>"
            "<p style='color:#5d6d7e;font-size:0.9em;margin:0 0 12px 0'>"
            "限界希釈プロトコール自動計算ツール</p>"),
        widgets.HTML("<b style='color:#2c3e50'>▸ 細胞サンプル</b>"),
        w_conc, w_vol,
        sep,
        widgets.HTML("<b style='color:#2c3e50'>▸ プレート設定</b>"),
        w_lambda, w_well_vol, w_n_wells,
        sep,
        widgets.HTML("<b style='color:#2c3e50'>▸ 希釈設定</b>"),
        w_max_df, w_min_tr,
        btn,
        out,
    ])

    def on_click(_):
        with out:
            clear_output(wait=True)
            try:
                result = calculate_protocol(
                    input_conc=w_conc.value,
                    input_vol_mL=w_vol.value,
                    lambda_val=w_lambda.value,
                    well_vol_uL=w_well_vol.value,
                    n_wells=int(w_n_wells.value),
                    max_step_df=int(w_max_df.value),
                    min_transfer_uL=w_min_tr.value,
                )
                display(HTML(render_html(result)))
                plot_result(result)
            except ValueError as e:
                display(HTML(
                    f"<div style='background:#fef9e7;border-left:4px solid #e74c3c;"
                    f"padding:10px;border-radius:4px'>❌ <b>エラー:</b> {e}</div>"))
            except Exception as e:
                display(HTML(
                    f"<div style='background:#fef9e7;border-left:4px solid #e74c3c;"
                    f"padding:10px;border-radius:4px'>❌ 予期しないエラー: {e}</div>"))

    btn.on_click(on_click)
    display(ui)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 6. CLI mode (standalone script)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def _cli_mode():
    print("=== Limiting Dilution Calculator (CLI mode) ===\n")
    try:
        conc   = float(input("細胞濃度 (cells/mL) [e.g. 1e6]: "))
        vol    = float(input("手持ちボリューム (mL) [e.g. 1.0]: "))
        lam    = float(input("λ (cells/well) [default 0.6]: ") or 0.6)
        wvol   = float(input("ウェル容量 µL [default 100]: ")  or 100)
        nw     = int(input("ウェル数 [default 96]: ")          or 96)
        maxdf  = int(input("1ステップ最大希釈倍率 [default 10]: ") or 10)
    except ValueError:
        print("入力値が無効です。終了します。")
        return

    result = calculate_protocol(conc, vol, lam, wvol, nw, maxdf)
    print_protocol(result)
    plot_result(result)


if __name__ == "__main__":
    _cli_mode()
