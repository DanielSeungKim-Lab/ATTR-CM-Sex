"""Reproduce every number reported in the research letter on female representation in landmark
hypertrophic cardiomyopathy (HCM) and transthyretin amyloid cardiomyopathy (ATTR-CM) trials.

Usage:  python CODE-01-reproduce-results.py      (requires numpy, pandas)
Reads DATA-01 (trial sex counts), DATA-02 (benchmarks) and DATA-04 (ATTR-CM trial design). Writes
RESULT-01 and RESULT-02, which CODE-02 uses to draw the Figure.

Primary estimator, both diseases: the participant-weighted proportion of women (women summed over
trials divided by participants summed over trials). Its 95% interval comes from a nonparametric
cluster bootstrap that resamples trials with replacement and ENUMERATES every ordered resample
(7^7 = 823,543 for HCM; 5^5 = 3,125 for ATTR-CM), so the trial distribution has no Monte Carlo
error. The interval describes uncertainty over trials that could have been run, not over the
observed counts, and is approximate with only 5 or 7 clusters.

Ratio intervals: 2.5th and 97.5th percentiles of 200,000 paired draws (seed 20261021).
  trial          resampled from the enumerated bootstrap distribution
  census         female and male totals drawn normal, SE = published margin of error / 1.645
                 (variances summed across age bands)
  binomial       counted benchmark, female share drawn Beta(f + 0.5, n - f + 0.5) (Jeffreys)
Within each disease the same referral draw enters both component ratios, so
referral-to-population x trial-to-referral = trial-to-population in every draw.
"""
import itertools
import numpy as np
import pandas as pd

SEED, NDRAW = 20261021, 200_000
trials = pd.read_csv("DATA-01-trial-sex-counts.csv")
bench = pd.read_csv("DATA-02-benchmarks.csv")
COHORT = {d: trials[trials.disease == d].reset_index(drop=True) for d in ["HCM", "ATTR-CM"]}
REFERRAL = {"HCM": "Specialty-center referral registry (SHaRe)",
            "ATTR-CM": "Specialty-center registry, ATTR stratum"}
POPULATION = {"HCM": "US adults 20 y and over", "ATTR-CM": "US adults 75 y and over"}
# ACS 2019-2023 B01001 national margins of error, combined across age bands; SE = MOE / 1.645
CENSUS_SE = {"US adults 75 y and over": (15498, 12233), "US adults 20 y and over": (38651, 36718)}


def pw(df):
    return 100 * df.female.sum() / df.sex_denom.sum()


def enumerated_bootstrap(df):
    f, n, k = df.female.to_numpy(float), df.sex_denom.to_numpy(float), len(df)
    idx = np.array(list(itertools.product(range(k), repeat=k)), dtype=np.int8)
    return 100 * f[idx].sum(1) / n[idx].sum(1)


def pct_ci(x):
    return np.percentile(x, 2.5), np.percentile(x, 97.5)


# ---- trial cohorts: participant-weighted proportions (Results) -----------------------------------
boot, table = {}, {}
for d, df in COHORT.items():
    boot[d] = enumerated_bootstrap(df)
    lo, hi = pct_ci(boot[d])
    per = 100 * df.female / df.sex_denom
    table[d] = dict(k=len(df), women=int(df.female.sum()), total=int(df.sex_denom.sum()),
                    pw=pw(df), pw_lo=lo, pw_hi=hi, trial_min=per.min(), trial_max=per.max())
    t = table[d]
    print(f"{d}: {t['women']} of {t['total']} enrollment records were women, {t['pw']:.1f}% "
          f"(95% CI, {t['pw_lo']:.1f}%-{t['pw_hi']:.1f}%), k={t['k']}; trial range {t['trial_min']:.1f}%-{t['trial_max']:.1f}%")

# ---- ratios: paired draws --------------------------------------------------------------------------
rng = np.random.default_rng(SEED)


def bench_draws(row):
    if row.benchmark in CENSUS_SE:
        fse, mse = CENSUS_SE[row.benchmark]
        fd = rng.normal(row.bench_female, fse, NDRAW)
        md = rng.normal(row.bench_denom - row.bench_female, mse, NDRAW)
        return 100 * fd / (fd + md), "census"
    return 100 * rng.beta(row.bench_female + 0.5, row.bench_denom - row.bench_female + 0.5, NDRAW), "binomial"


trial_draws = {d: rng.choice(boot[d], NDRAW, replace=True) for d in COHORT}
draws, ratio_rows = {}, []
for _, b in bench.iterrows():
    den, model = bench_draws(b)
    draws[(b.disease, b.benchmark)] = den
    obs = table[b.disease]["pw"]
    lo, hi = pct_ci(trial_draws[b.disease] / den)
    ratio_rows.append((b.disease, b.benchmark, b.short_label, b.kind, model, b.bench_pct, obs / b.bench_pct, lo, hi))
R = pd.DataFrame(ratio_rows, columns=["disease", "benchmark", "short_label", "kind", "den_model", "bench_pct",
                                      "ppr", "ppr_lo", "ppr_hi"])
strata = R[R.kind == "disease_stratum"]
R = R[R.kind != "disease_stratum"]
R.to_csv("RESULT-02-benchmark-ratios.csv", index=False)
print("\nTrial-to-benchmark ratios (Figure B):")
print(R[["disease", "short_label", "bench_pct", "ppr", "ppr_lo", "ppr_hi"]].round(3).to_string(index=False))

# ---- decomposition: population -> referral registry -> trials (Figure A) --------------------------
dec = []
for d in COHORT:
    pop_row = bench[bench.benchmark == POPULATION[d]].iloc[0]
    ref_row = bench[bench.benchmark == REFERRAL[d]].iloc[0]
    pop, ref, tri = draws[(d, POPULATION[d])], draws[(d, REFERRAL[d])], trial_draws[d]
    p, r, t = pop_row.bench_pct, ref_row.bench_pct, table[d]["pw"]
    rp, tr, tp = ref / pop, tri / ref, tri / pop
    assert np.allclose(rp * tr, tp)
    dec.append((d, p, r, ref_row.benchmark.split("(")[-1].rstrip(")") if d == "HCM" else "CARS",
                int(ref_row.bench_female), int(ref_row.bench_denom), t,
                r / p, *pct_ci(rp), t / r, *pct_ci(tr), t / p, *pct_ci(tp)))
    print(f"\n{d}: population {p:.1f}% -> referral {r:.1f}% -> trials {t:.1f}%")
    print(f"   referral-to-population {r/p:.2f} ({np.percentile(rp,2.5):.2f}-{np.percentile(rp,97.5):.2f}); "
          f"trial-to-referral {t/r:.2f} ({np.percentile(tr,2.5):.2f}-{np.percentile(tr,97.5):.2f}); "
          f"trial-to-population {t/p:.2f} ({np.percentile(tp,2.5):.2f}-{np.percentile(tp,97.5):.2f})")
pd.DataFrame(dec, columns=["disease", "population_pct", "referral_pct", "referral_source", "referral_female", "referral_n",
                           "trial_pct", "referral_to_population", "ref_pop_lo", "ref_pop_hi", "trial_to_referral",
                           "trial_ref_lo", "trial_ref_hi", "trial_to_population", "trial_pop_lo", "trial_pop_hi"]
             ).to_csv("RESULT-01-pipeline-decomposition.csv", index=False)

# ---- claims cohort by genotype (Results) -----------------------------------------------------------
print("\nClaims cohort by genotype:")
print(strata[["short_label", "bench_pct", "ppr", "ppr_lo", "ppr_hi"]].round(3).to_string(index=False))

# ---- spread across the ATTR-CM disease benchmarks named in the letter ------------------------------
named = R[(R.disease == "ATTR-CM") & (R.kind == "disease")]
assert len(named) == 2, "expected the 2 ATTR-CM disease benchmarks"
print(f"\nATTR-CM disease benchmarks differ by a factor of {named.ppr.max() / named.ppr.min():.1f}")

# ---- variant (hereditary) disease: trials vs referral registry (Limitations) ---------------------
G = pd.read_csv("DATA-05-genotype-mix.csv")
tv, tg = G[G.kind == "trial"].variant.sum(), G[G.kind == "trial"].genotyped.sum()
reg = G[G.kind == "benchmark"].iloc[0]
print(f"Variant ATTR-CM: trials {tv} of {tg} ({100 * tv / tg:.0f}%), "
      f"referral registry {reg.variant} of {reg.genotyped} ({100 * reg.variant / reg.genotyped:.0f}%)")


# ---- ATTR-CM trial design and efficacy in women (Table) ----------------------------------------------
# Table cells are abstracted from registry records and published methods (DOC-01); this block recomputes
# the values derived from them.
D = pd.read_csv("DATA-04-attr-trial-design.csv")
comp, ong = D[D.status == "completed"], D[D.status == "ongoing"]
print("\nATTR-CM trial design")
for r in comp.itertuples():
    print(f"   {r.trial}: {int(r.women)}/{int(r.n):,} women ({100 * r.women / r.n:.1f}%)")
print(f"   ongoing pivotal trials: {len(ong)}, planned participants {int(ong.n.sum()):,}")
print(f"   unindexed wall-thickness criterion: completed {int(comp.wall_rule.ne('Not stated').sum())} of {len(comp)}, "
      f"ongoing {int(ong.wall_rule.ne('Not stated').sum())} of {len(ong)}")
print(f"   upper age limit: completed {int(comp.upper_age.notna().sum())} of {len(comp)}, ongoing {int(ong.upper_age.notna().sum())} of {len(ong)}")
print(f"   enrollment target for women: {int(D.target_women.eq('yes').sum())} of {len(D)}")
pub = comp.dropna(subset=["women_est"])
print(f"   primary end-point estimates in women published: {len(pub)} of {len(comp)}")
for r in pub.itertuples():
    w = np.log(r.women_hi / r.women_lo) / np.log(r.men_hi / r.men_lo)
    print(f"   {r.trial}: women {r.women_est:.2f} ({r.women_lo:.2f}-{r.women_hi:.2f}) vs men {r.men_est:.2f} "
          f"({r.men_lo:.2f}-{r.men_hi:.2f}); interval width ratio on the log scale {w:.1f}")
