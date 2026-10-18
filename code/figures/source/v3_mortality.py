# Release adaptation: paths only unless noted in docs/RELEASE_CHANGES.md.
import os as _release_os
def _release_path(env, relative=''):
    root = _release_os.environ.get(env)
    if not root: raise RuntimeError('Set environment variable '+env+'; see docs/REPRODUCIBILITY.md')
    return _release_os.path.join(root, relative)
"""Render the frozen long-term mortality figure; no model fitting."""

import re
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

import base_plotting as b


def make_mortality():
    """Create Candidate_02_LongTerm_Mortality PNG (300 dpi) and PDF."""
    groups = ["G1", "G2", "G3"]
    counts = {
        g: b.one("main", "Group counts", threshold_mgl=30, group=g)
        for g in b.GROUPS
    }
    effects = {
        contrast: b.one(
            "main", "Cox evidence", scenario="primary", model="B",
            design="survey_Taylor", contrast=contrast
        )
        for contrast in ["G1_G0", "G2_G0", "G3_G0", "G2_G1"]
    }
    assert sum(r["ccB"] for r in counts.values()) == 5155
    assert sum(r["death_ccB"] for r in counts.values()) == 1115
    assert all(r["n"] == 5155 and r["events"] == 1115 for r in effects.values())
    periods = ["0_60", "after60"]
    period_rows = {
        (g, p): b.one("main", "Period Cox", group=g, period=p)
        for g in groups for p in periods
    }
    # Same frozen displayed P, read from a public aggregate instead of private handoff text.
    import json
    interaction_p = json.loads((b.ROOT / 'Inputs' / 'period_interaction.json').read_text())['display_input_p']
    assert abs(interaction_p - 0.17983) < 1e-10

    fig = plt.figure(figsize=(12, 8.4))
    b.header(
        fig, "Directional classification and long-term mortality",
        "NHANES 1999–2018 · Survey-weighted Cox Model B · n = 5,155; deaths = 1,115"
    )
    b.panel(fig, .055, .859, "A", "Whole follow-up: adjusted associations")

    # Explicit counts refer to the complete-case population used by Model B.
    ax = fig.add_axes([.555, .535, .245, .275])
    ax.set_ylim(-1.02, 4.68)
    b.logx(ax, [.5, 1, 2, 3], (.5, 3.2))
    ax.set_yticks([])
    ax.set_xlabel("Hazard ratio (log scale)", labelpad=6)
    trans = ax.get_yaxis_transform()
    label_x = (.055 - .555) / .245
    n_x = (.398 - .555) / .245
    deaths_x = (.472 - .555) / .245
    value_x = (.835 - .555) / .245
    for x, text in [(label_x, "Classification"), (n_x, "n"),
                    (deaths_x, "Deaths"), (value_x, "HR (95% CI)")]:
        ax.text(x, 4.56, text, transform=trans, fontsize=10,
                weight="bold", va="center", clip_on=False)
    for g, y in zip(b.GROUPS, [4, 3, 2, 1]):
        r = counts[g]
        ax.text(label_x, y, f"{g}  {b.NAMES[g]}", transform=trans,
                color=b.C[g], fontsize=11, va="center", clip_on=False,
                weight="bold" if g == "G2" else "normal")
        ax.text(n_x, y, f"{r['ccB']:,}", transform=trans,
                fontsize=10.5, va="center", clip_on=False)
        ax.text(deaths_x, y, f"{r['death_ccB']:,}", transform=trans,
                fontsize=10.5, va="center", clip_on=False)
        if g == "G0":
            ax.plot(1, y, "D", color=b.C[g], ms=5.8)
            label = "1.00 (reference)"
        else:
            effect = effects[f"{g}_G0"]
            b.interval(ax, y, effect, color=b.C[g])
            label = b.fmt(effect)
        ax.text(value_x, y, label, transform=trans, fontsize=11,
                va="center", clip_on=False,
                weight="bold" if g == "G2" else "normal")
    ax.axhline(.36, color=b.GRID, lw=.9)
    direct_y = -.42
    ax.text(label_x, direct_y, "G2 vs G1 · direct contrast", transform=trans,
            fontsize=11, va="center", clip_on=False)
    b.interval(ax, direct_y, effects["G2_G1"], color=b.INK,
               marker="D", openmark=True)
    ax.text(value_x, direct_y, b.fmt(effects["G2_G1"]), transform=trans,
            fontsize=11, va="center", clip_on=False)

    b.panel(fig, .055, .443, "B", "Time-specific associations: frozen sensitivity model")
    fig.text(.084, .403, "Each cell compares the indicated group with G0 in that period.",
             fontsize=10, color=b.MUTED, va="top")
    heat_ax = fig.add_axes([.327, .161, .495, .194])
    heat_ax.set_xlim(0, 2)
    heat_ax.set_ylim(3, 0)
    cmap = LinearSegmentedColormap.from_list(
        "mortality_hr", ["#F1F7F9", "#BAD7E2", "#76A7BB", "#2A647E"]
    )
    norm = Normalize(vmin=.75, vmax=2.5)
    values = np.array([[period_rows[(g, p)]["HR"] for p in periods] for g in groups])
    heat_ax.imshow(values, extent=(0, 2, 3, 0), aspect="auto",
                   cmap=cmap, norm=norm, interpolation="nearest")
    heat_ax.set_xticks([.5, 1.5], ["0–60 months", ">60 months"])
    heat_ax.xaxis.tick_top()
    heat_ax.tick_params(axis="x", length=0, pad=10, labelsize=10.5)
    heat_ax.set_yticks(np.arange(3) + .5,
                      [f"{g}  {b.NAMES[g]}" for g in groups])
    heat_ax.tick_params(axis="y", length=0, pad=13, labelsize=10.5)
    for tick, g in zip(heat_ax.get_yticklabels(), groups):
        tick.set_color(b.C[g])
        tick.set_weight("bold" if g == "G2" else "normal")
    for spine in heat_ax.spines.values():
        spine.set_visible(False)
    heat_ax.axvline(1, color="white", lw=3)
    for boundary in [1, 2]:
        heat_ax.axhline(boundary, color="white", lw=3)
    for i, g in enumerate(groups):
        for j, p in enumerate(periods):
            r = period_rows[(g, p)]
            text_color = "white" if r["HR"] >= 1.9 else b.INK
            heat_ax.text(j + .5, i + .34, f"{r['HR']:.2f}",
                         ha="center", va="center", fontsize=13,
                         weight="bold", color=text_color)
            heat_ax.text(j + .5, i + .71, f"({r['low']:.2f}–{r['high']:.2f})",
                         ha="center", va="center", fontsize=10, color=text_color)
    cax = fig.add_axes([.865, .161, .015, .194])
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                           ticks=[.75, 1, 1.5, 2, 2.5])
    colorbar.ax.set_yticklabels(["0.75", "1.0", "1.5", "2.0", "2.5"])
    colorbar.ax.tick_params(labelsize=9, width=.6, length=3)
    colorbar.set_label("Hazard ratio", fontsize=10, labelpad=8)
    colorbar.outline.set_visible(False)
    fig.text(.327, .128,
             f"Joint group × period interaction: P = {interaction_p:.3f}",
             fontsize=10.5, color=b.MUTED, ha="left", va="center")
    b.note(
        fig,
        "G2 vs G0 is the original primary comparison. Panel B reports HRs and 95% CIs from a separate frozen sensitivity model.\n"
        "G0: Alb−/UACR−; G1: Alb+/UACR−; G2: Alb−/UACR+; G3: Alb+/UACR+. Thresholds: 30 mg/L and 30 mg/g.",
        y=.037
    )
    stem = b.MAIN / "Candidate_02_LongTerm_Mortality"
    temporary_png = stem.with_suffix(".tmp.png")
    temporary_pdf = stem.with_suffix(".tmp.pdf")
    fig.savefig(temporary_png, dpi=300, facecolor="white")
    with Image.open(temporary_png) as check:
        check.verify()
    fig.savefig(temporary_pdf, facecolor="white")
    temporary_png.replace(stem.with_suffix(".png"))
    temporary_pdf.replace(stem.with_suffix(".pdf"))
    plt.close(fig)
    return stem.with_suffix(".png"), stem.with_suffix(".pdf")


if __name__ == "__main__":
    make_mortality()
