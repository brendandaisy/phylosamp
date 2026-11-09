# phylosamp (development version)

* Added `varfreq_time_logistic()` to calculate the time at which a variant reaches a given prevalence under logistic growth (the inverse of `varfreq_freq_logistic()`)
* Added `vartrack_samplesize_prev_cont()` and enabled `sampling_freq = "cont"` in `vartrack_samplesize_prev()` to calculate the sample size needed for variant prevalence estimation under periodic sampling
* Added `vartrack_samplesize_growth()` to calculate the sample size needed for estimating variant growth advantage under periodic sampling
* Added `vartrack_tmin_growth()` and `vartrack_tmin_prev()` to calculate the time needed to estimate variant growth advantage or prevalence given a per-timestep sample size
* Added `vartrack_samplesize_dominance()` and `vartrack_tmin_dominance()` to calculate the sample size or sampling time needed to estimate when a variant will become dominant (reach 50% prevalence) under periodic sampling
* Added `vartrack_samplesize_growth_power()` and `vartrack_power_growth()` for sample size and power calculations when testing for a variant growth advantage under periodic sampling

# phylosamp 1.0.1

* Added a `NEWS.md` file to track changes to the package.
* Added `relR_*` functions to the package to allow for sample size and power calculations when estimating differential transmission between groups