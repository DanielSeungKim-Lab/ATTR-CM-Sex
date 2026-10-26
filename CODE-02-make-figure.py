"""Regenerate the letter's Figure (two panels) from the CODE-01 outputs.

Usage: python CODE-02-make-figure.py
Requires pandas, numpy, matplotlib. Writes figure1_sex_pipeline.png at 300 dpi.
Run CODE-01 first. Color marks disease; in panel B marker shape marks the kind of benchmark.
"""
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt

S = pd.read_csv("RESULT-01-pipeline-decomposition.csv")
RES = pd.read_csv("RESULT-02-benchmark-ratios.csv")
RES["short_label"] = RES.short_label.replace({"Diagnosed claims, all": "Diagnosed claims cohort"})
assert len(RES) == 5, "expected the 5 benchmarks named in the text"
DIS = {"HCM": "#1b6ca8", "ATTR-CM": "#c0392b"}
MK = {"population": "s", "disease": "o"}

fig = plt.figure(figsize=(7.4, 6.4))
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.72], hspace=0.34)

axA = fig.add_subplot(gs[0]); x = np.arange(3)
VOFF = {"HCM": (-13, -15, 10), "ATTR-CM": (10, 12, 10)}
DSTEP = {"HCM": 19, "ATTR-CM": -21}
for _, r in S.iterrows():
    c = DIS[r.disease]; y = [r.population_pct, r.referral_pct, r.trial_pct]
    axA.plot(x, y, "-o", color=c, lw=1.8, ms=7, zorder=3, clip_on=False)
    for xi, yi, off in zip(x, y, VOFF[r.disease]):
        axA.annotate(f"{yi:.1f}%", (xi, yi), xytext=(0, off), textcoords="offset points",
                     ha="center", fontsize=6.8, color=c, fontweight="bold", zorder=5)
    for xm, lab, val, lo, hi, off in [
            (0.5, "referral-to-population", r.referral_to_population, r.ref_pop_lo, r.ref_pop_hi,
             DSTEP[r.disease]),
            (1.5, "trial-to-referral", r.trial_to_referral, r.trial_ref_lo, r.trial_ref_hi, 15)]:
        ym = (y[0] + y[1]) / 2 if xm == 0.5 else (y[1] + y[2]) / 2
        axA.annotate(f"{lab}\n{val:.2f} ({lo:.2f} to {hi:.2f})", (xm, ym), xytext=(0, off),
                     textcoords="offset points",
                     ha="center", va="center", fontsize=6.6, color=c, zorder=5,
                     bbox=dict(fc="white", ec="none", pad=0.9, alpha=0.9))
    axA.annotate(r.disease, (2, y[2]), xytext=(11, 0), textcoords="offset points",
                 va="center", fontsize=7.6, color=c, fontweight="bold")
axA.set_xticks(x)
axA.set_xticklabels(["General population,\ndisease-specific age band", "Specialty-center\nreferral registry",
                     "Landmark\ntrials"], fontsize=6.9)
axA.set_xlim(-0.28, 2.6); axA.set_ylim(0, 68)
axA.set_ylabel("Women (% of the population)")

axB = fig.add_subplot(gs[1]); rows, y0 = [], 0
for dis in ["HCM", "ATTR-CM"]:
    for _, r in RES[RES.disease == dis].sort_values("ppr").iterrows():
        rows.append((y0, r, dis)); y0 += 1
    y0 += 1.3
ymax = y0
axB.axvline(1.0, color="0.35", lw=1.0, zorder=1)
for yy, r, dis in rows:
    c = DIS[dis]
    axB.plot([r.ppr_lo, r.ppr_hi], [ymax - yy] * 2, color=c, lw=1.5, solid_capstyle="round", zorder=2)
    axB.plot(r.ppr, ymax - yy, MK[r.kind], color=c, ms=5.5, zorder=3,
             mfc=c, mew=1.3)
    # labels go right of the interval unless they would sit on the reference line at 1
    left = 0.8 < r.ppr_hi < 1.0
    axB.annotate(f"{r.ppr:.2f}", (r.ppr_lo if left else r.ppr_hi, ymax - yy), xytext=(-4 if left else 4, 0),
                 textcoords="offset points", ha="right" if left else "left", va="center", fontsize=6.4, color=c)
axB.set_yticks([ymax - yy for yy, _, _ in rows])
axB.set_yticklabels([f"{r.short_label}  ({r.bench_pct:.1f}%)" for _, r, _ in rows], fontsize=6.4)
for t, (_, _, dis) in zip(axB.get_yticklabels(), rows):
    t.set_color(DIS[dis])
axB.set_xscale("log"); axB.set_xlim(0.105, 3.4)
axB.set_xticks([0.125, 0.25, 0.5, 1, 2])
axB.get_xaxis().set_major_formatter(mpl.ticker.FixedFormatter(["0.125", "0.25", "0.5", "1", "2"]))
axB.get_xaxis().set_minor_formatter(mpl.ticker.NullFormatter())
axB.set_xlabel("Trial-to-benchmark ratio for female representation")
hs = [plt.Line2D([], [], ls="none", marker=m, color="0.35", ms=5.5, label=l) for m, l in
      [("s", "General population"), ("o", "Referral registry or diagnosed patients")]]
axB.legend(handles=hs, loc="upper left", bbox_to_anchor=(0.01, 0.99), frameon=False,
           fontsize=6.3, handletextpad=0.5, labelspacing=0.55)
for ax in (axA, axB):
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
for ax, L in [(axA, "A"), (axB, "B")]:
    ax.text(-0.20, 1.04, L, transform=ax.transAxes, fontweight="bold", fontsize=16, va="bottom")
fig.savefig("figure1_sex_pipeline.png", dpi=300, bbox_inches="tight")
print("wrote figure1_sex_pipeline.png")
