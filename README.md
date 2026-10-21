
# EpiStrainDynamics

[![codecov](https://codecov.io/gh/ropensci/EpiStrainDynamics/graph/badge.svg)](https://app.codecov.io/gh/ropensci/EpiStrainDynamics)
[![](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-CMD-check](https://github.com/ropensci/EpiStrainDynamics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ropensci/EpiStrainDynamics/actions/workflows/R-CMD-check.yaml)
[![Status at rOpenSci Software Peer
Review](https://badges.ropensci.org/763_status.svg)](https://github.com/ropensci/software-review/issues/763)

`EpiStrainDynamics` is a statistical modelling framework capable of
inferring trends of multiple pathogens. Estimating the temporal trends
in infectious disease activity is crucial for monitoring disease spread
and the impact of interventions. Surveillance indicators routinely
collected to monitor these trends are often a composite of multiple
pathogens. For example, ‘influenza-like illness’ — routinely monitored
as a proxy for influenza infections — is a symptom definition that could
be caused by a wide range of pathogens, including multiple subtypes of
influenza, SARS-CoV-2, and RSV. Inferred trends from such composite time
series may not reflect the trends of any one of the component pathogens,
each of which can exhibit distinct dynamics. Although many surveillance
systems routinely test a subset of individuals contributing to a
surveillance indicator — providing information on the relative
contribution of the component pathogens — trends may be obscured by
time-varying testing rates or substantial noise in the observation
process.

`EpiStrainDynamics` builds on existing modelling frameworks built to
handle two pathogens ([Eales et
al. 2022](https://doi.org/10.1038/s41467-022-32096-4); [Eales et
al. 2022](https://doi.org/10.1016/j.epidem.2022.100604)), and extends
them be to able to:

- infer the trends of any number of pathogens
- fit to time series data of counts (eg, daily number of cases)
- incorporate influenza testing data in which the subtype for influenza
  A samples may be undetermined
- account for day-of-the-week effects in daily data
- include options for fitting penalized splines or random walks
- support additional (optional) correlation structures in the parameters
  describing the smoothness of the penalized splines (or random walks),
  and
- account for additional (optional) sources of noise in the observation
  process.

## Installation

You can install the development version of EpiStrainDynamics from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("ropensci/EpiStrainDynamics")
```

As there is not yet a CRAN release, this always builds
`EpiStrainDynamics` from source, so users will need to configure their
C++ toolchain first. This is because `EpiStrainDynamics` implements the
underlying models in Stan (a statistical modelling language built on
C++).

Alternatively, thanks to
[R-universe](https://ropensci.r-universe.dev/EpiStrainDynamics), you can
install a pre-built binary without needing to configure a C++ toolchain:

``` r
install.packages("EpiStrainDynamics", repos = c('https://ropensci.r-universe.dev', 'https://cloud.r-project.org'))
```

Each operating system has a different set up procedure. Windows users
need to install an appropriate version of
[RTools](https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows).
Mac users can [follow these
steps](https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Mac),
and Linux users can use [this
guide](https://github.com/stan-dev/rstan/wiki/Configuring-C-Toolchain-for-Linux).

## Using `EpiStrainDynamics`

Detailed instructions can be found on the
[vignette](https://docs.ropensci.org/EpiStrainDynamics/articles/Using-EpiStrainDynamics.html).
Here we provide a short overview.

A full analysis proceeds through four phases, each producing an object
consumed by the next: preparing the pathogen data, configuring the
model, fitting it, and deriving epidemiological quantities. The sections
below walk through each phase using the `sarscov2` data bundled with the
package.

### Phase 1: Preparing pathogen data

There are three main types of pathogen structure available to model:
`single()`, `multiple()`, and `subtyped()`. These functions require the
name of the dataset itself and the column names for different data
elements.

The `single()` pathogen structure is the simplest and models a single
pathogen timeseries. The user must specify the dataset name, the name of
the column with total case data, and the name of the column of time
data. The `multiple()` pathogen structure allows modelling of different
component pathogens. In addition to specifying the dataset, total case
data, and time, the additional pathogens are specified as a vector of
column names. The `subtyped()` pathogen structure enables additional
complexity for a scenario where the user has a total timeseries for a
parent pathogen alongside data on some of its underlying subtypes — most
commonly an influenza modelling scenario with testing data for different
influenza A subtypes. The user specifies columns containing the total
case count for the parent pathogen, including cases that were
subsequently subtyped, as well as the subtyped cases, and any additional
pathogens to be modelled. See the vignette for further detail.

Each of these functions returns a pathogen structure object, which is
then passed into `construct_model()` in Phase 2.

### Phase 2: Configuring the model

Modelling specifications are provided using the `construct_model()`
function. The correct stan model is then applied based on the
specifications provided. `construct_model()` takes the pathogen
structure object from Phase 1, plus these arguments describing the
model: `method`, `smoothing_params`, `dispersion_params`,
`pathogen_noise`, and `dow_effect`.

#### Method

EpiStrainDynamics has pre-compiled stan models that fit either with
bayesian penalised splines or random walks. These are specified using
the `method` argument of `construct_model()` as functions, either with
`random_walk()` or `p_spline()`. The penalised spline model has two
further options to specify: `spline_degree` is the polynomial degree of
the individual spline segments used to construct the overall curve (must
be a positive whole number) and `days_per_knot`, which is the number of
days for each knot (must also be a positive whole number).

So we may specify:

    method = random_walk(),
    # OR
    method = p_spline(spline_degree = 3, days_per_knot = 2)

#### Smoothing parameters

The argument `smoothing_params` allows users to modify the correlation
structures in the parameters describing smoothness and to set related
priors. These are specified with the function `smoothing_structure()`,
which requires the user to specify a `smoothing_type` that is either
`shared` (all pathogens have the same smoothing parameter),
`independent` (each pathogen has a completely independent smoothing
parameter), or `correlated` (smoothing structure is correlated among
pathogens). For smoothign types `shared` and `independent` the priors
for the mean and standard deviation on tau can also be specified. The
number of values specified for each parameter depend on how many
pathogens are provided. As below:

    smoothing_structure(
      smoothing_type = 'shared',
      tau_mean = 0, tau_sd = 1
    )
    # for a model with 4 parameters:
    smoothing_structure(
      smoothing_type = 'independent',
      tau_mean = c(0, 0, 0, 0),
      tau_sd = c(1, 1, 1, 1)
    )

#### Dispersion parameters

The argument `dispersion_params` optionally allows users to set a prior
for the overdispersion parameter of the negative binomial likelihood for
the case timeseries. If specified, each parameter only ever needs a
single value. It is specified using `dispersion_structure()` as below:

    dispersion_structure(
      phi_mean = 0, phi_sd = 1
    )

#### Pathogen noise

A logical (`TRUE` or `FALSE`) value indicating whether to include noise
between individual pathogens in addition to the observation noise.

#### Day of week effect

Day of week effect is specified as a logical (`TRUE` or `FALSE`) to the
`dow_effect` argument. In plotting, the day of week effect can be
selectively removed.

#### Example

Altogether, an example constructed model for a random walk model with
multiple pathogens might look as follows, illustrated using data that
has been provided with the package - `sarscov2`:

    sarscov2_multi <- multiple(                 # multiple pathogen structure
       data = sarscov2,
       case_timeseries = 'cases',               # timeseries of case data
       time = 'date',                           # date or time variable labels
       component_pathogen_timeseries = c(       # component pathogens
         'alpha', 'delta', 'omicron', 'other'
       )
    )

    mod <- construct_model(
      pathogen_structure = sarscov2_multi,

      method = random_walk(),                   # random_walk method

      smoothing_params = smoothing_structure(   # independent smoothing structure
        'independent',                          # with four values for each prior
        tau_mean = c(0, 0.1, 0.3, 0),           # parameter - one for each pathogen
        tau_sd = rep(1, times = 4)
      ),
      dispersion_params = dispersion_structure(phi_mean = 0, phi_sd = 1),
      pathogen_noise = FALSE,
      dow_effect = TRUE
    )

### Phase 3: Fitting

The model estimates the expected value of the time series (eg, a
smoothed trend in the daily number of cases accounting for noise) for
each individual pathogen. In this step additional fitting parameters
related to stan models can be modified such as `n_chain`, `n_iter`,
`n_warmup`, `thin`, `adapt_delta`, `multi_cores`, `verbose`, and `seed`,
which are described in further detail in the documentation. Model
parameterisation decisions specified when configuring the model in Phase
2 mean the correct stan model will be applied at this stage by simply
calling:

    fit <- fit_model(mod)

Convergence can then be checked with `diagnose_model(fit)`; see the
vignette for guidance on interpreting and troubleshooting these
diagnostics.

### Phase 4: Deriving epidemiological quantities

Calculate epidemic growth rate with `growth_rate(fit)`, effective
reproduction number over time with `Rt(fit, gi_dist = X)` (requiring
specification of a generation interval distribution X), incidence with
or without a day of week effect with
`incidence(fit, dow_effect = TRUE)`, and proportions of different
combinations of cases attributable to different pathogens/subtypes using
`proportion()`.

    prop <- proportion(fit)
    plot(prop)

For a more detailed discussion, check out the
[vignette](https://docs.ropensci.org/EpiStrainDynamics/articles/Using-EpiStrainDynamics.html).

## Citation

When using this package, please cite both the package and the research
article underlying the statistical model developments:

Windecker S, Eales O (2025). EpiStrainDynamics: Infer temporal trends of
multiple pathogens. R package version 0.0.1,
<https://docs.ropensci.org/EpiStrainDynamics/>.

Oliver Eales, Saras M Windecker, James M McCaw, Freya M Shearer,
Inferring temporal trends of multiple pathogens, variants, subtypes or
serotypes from routine surveillance data, American Journal of
Epidemiology, 2025;, kwaf119, <https://doi.org/10.1093/aje/kwaf119>

For code corresponding to the AJE paper, see branch
[`paper_analysis`](https://github.com/ropensci/EpiStrainDynamics/tree/paper_analysis).

## Contribution

`EpiStrainDynamics` is in a stable state of development, with some
degree of active subsequent development as envisioned by the primary
authors. Authors are committed to maintaining the package. If you spot a
bug, have a feature request, or want to contribute an improvement to the
package (branch `main`) or the code associated with the paper analyses
(branch `paper_analysis`), please open an
[issue](https://github.com/ropensci/EpiStrainDynamics/issues) or pull
request. The [contributing
guide](https://github.com/ropensci/EpiStrainDynamics/blob/main/.github/CONTRIBUTING.md)
describes the recommended workflow and scope for contributions. Test
cases and improvements are especially welcome.

## Code of Conduct

Please note that the `EpiStrainDynamics` project is released with a
[Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Support

This project was supported by the Australia-Aotearoa Consortium for
Epidemic Forecasting and Analytics.

<a href="https://acefa-hubs.github.io/"><img src="man/figures/ACEFA.png" align = "center" height="150" alt="EpiStrainDynamics website" /></a>

<!-- `devtools::build_readme()` is handy for this. -->
