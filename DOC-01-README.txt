SEX REPRESENTATION IN LANDMARK HYPERTROPHIC AND TRANSTHYRETIN AMYLOID CARDIOMYOPATHY TRIALS
DATA AND CODE

RUN
  python CODE-01-reproduce-results.py   prints every computed value in the letter; writes RESULT-01
                                        and RESULT-02
  python CODE-02-make-figure.py         draws the Figure from RESULT-01 and RESULT-02
  Python 3.11 with numpy, pandas and matplotlib. Output is deterministic
  (enumerated bootstrap; ratio draws use seed 20261021).

FILES
  DATA-01-trial-sex-counts.csv   one row per included trial (7 HCM, 5 ATTR-CM). Columns: disease,
      trial, nct, female, male, sex_denom (participants with sex reported), sex_source
      (registry baseline module or primary-report table).
  DATA-02-benchmarks.csv   the 5 benchmarks plus the claims cohort's wild-type and variant strata.
      Columns: disease, benchmark, kind (population, disease, or disease_stratum),
      bench_female, bench_denom, bench_pct, source (citation and table), short_label (Figure label).
      Census margins of error are in CODE-01 (CENSUS_SE).
  DATA-03-trial-screening.csv   every ClinicalTrials.gov trial meeting the design criteria in the
      HCM and ATTR searches, with its disposition. Columns: disease, stage, disposition, nct,
      acronym, title, phases, enroll, status, us_sites, has_results.
  DATA-04-attr-trial-design.csv   the Table: completed and ongoing pivotal ATTR-CM trials.
      Columns: trial, nct, status, women, n (enrolled, or planned for ongoing trials), upper_age
      (blank if none), wall_rule, target_women, other_target, efficacy_measure,
      women_est/women_lo/women_hi and men_est/men_lo/men_hi (primary end point and 95% CI by sex;
      blank if not reported), efficacy_source, criteria_source.
  DATA-05-genotype-mix.csv   variant (hereditary) ATTR-CM counts for the 5 ATTR-CM trials and the
      referral registry. Columns: population, kind (trial or benchmark), variant, genotyped
      (participants with genotype reported), source (citation and table).
  LICENSE   MIT License, for the code.
  LICENSE-DATA   CC BY 4.0, for the DATA files.

SOURCES
  ClinicalTrials.gov API v2. Completed-trial eligibility also uses
  published methods (doi:10.1056/NEJMoa1805689, doi:10.1056/NEJMoa2300757,
  doi:10.1056/NEJMoa2305434, doi:10.1056/NEJMoa2409134, doi:10.1056/NEJMoa2608510,
  doi:10.1038/s41591-025-03851-z) and the protocols and statistical analysis plans
  published with the first three. Genotype counts are in DATA-05.
  Efficacy in women: doi:10.1056/NEJMoa2608510 (Figure 2) and doi:10.1093/ejhf/xuag208.
  Benchmark sources are in DATA-02. 