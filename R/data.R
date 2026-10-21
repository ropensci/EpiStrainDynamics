#' United Kingdom Health Security Agency SARS-CoV-2 case data
#'
#' Daily SARS-CoV-2 case numbers by specimen date for the United Kingdom from
#' 2020 - 2022 retrieved from the UK Health Security Agency's data dashboard.
#' downloaded data describing the daily number of variants detected by
#' collection date. The data classified all sequences based on ‘major lineage
#' calls’. We considered groupings of the variants based on their major
#' lineage calls or other, consisting of all lineages with a designation not
#' consistent with any of the major lineage calls.
#'
#' @format ## `sarscov2`
#' A data frame with 830 rows and 6 columns:
#' \describe{
#'   \item{date}{Date, between 23 September 2020 and 31 December 2022}
#'   \item{cases}{Numeric, daily number of cases}
#'   \item{alpha}{Integer, daily number of cases of B.1.1.7 (Alpha variant)}
#'   \item{delta}{Integer, daily number of cases of B.1.617.2 (Delta variant)}
#'   \item{omicron}{Integer, daily number of cases of BA.1, BA.2, BA.2.75,
#'      BA.4, BA.5, BQ.1 (Omicron variants)}
#'   \item{other}{Integer, daily number of cases of B.1.177, XBB (a recombinant
#'      of omicron sub-variants), and all other lineages with a designation
#'      not consistent with any of the major lineage calls}
#' }
#' @source <https://ukhsa-dashboard.data.gov.uk/covid-19-archive-data-download;
#'   https://datadryad.org/dataset/doi:10.5061/dryad.hx3ffbgm2>
#' @srrstats {G5.1} Dataset `sarscov2` is exported and described in detail here.
#'   It is used in both example code and in testing.
#' @srrstats {G1.4} uses `Roxygen2` documentation
"sarscov2"

#' World Health Organisation Global Influenza Programme for Australia
#'
#' Influenza-like illness data were retrieved from the World Health
#' Organization’s Global Influenza Programme for Australia, and subset to the
#' weeks starting January 2, 2012 to February 24, 2020, before the COVID-19
#' pandemic disrupted influenza transmission and surveillance. The data
#' described: (1) the weekly number of cases of influenza-like illness; and
#' (2) the weekly number of laboratory specimens tested, and of those positive
#' for influenza by type and subtype. The specimen counts come from laboratory
#' surveillance separate from the influenza-like illness counts, and are not
#' a subset of them. We grouped the specimens into: influenza A (all
#' influenza A positives, whether or not subtyped), of which some were
#' subtyped as influenza A H3N2 or influenza A H1N1 (influenza A H1N1,
#' influenza A H1N1pdm09); influenza B (influenza B Yamagata, influenza B
#' Victoria, influenza B lineage not determined); and other (specimens
#' tested that were not positive for influenza).
#'
#' @format ## `influenza`
#' A data frame with 426 rows and 7 columns:
#' \describe{
#'   \item{ili}{Integer, weekly number of cases of influenza-like illness}
#'   \item{week}{Date, start of the week, from 2 January 2012 to 24 February 2020}
#'   \item{inf_A}{Integer, weekly number of specimens positive for influenza
#'     A, including those subsequently subtyped as H3N2 or H1N1}
#'   \item{inf_B}{Integer, weekly number of specimens positive for influenza B}
#'   \item{inf_H3N2}{Integer, weekly number of specimens positive for
#'     influenza A subtype H3N2}
#'   \item{inf_H1N1}{Integer, weekly number of specimens positive for
#'     influenza A subtype H1N1}
#'   \item{other}{Integer, weekly number of specimens tested that were not
#'     positive for influenza, calculated as the number of specimens tested
#'     minus the number positive for influenza A or B}
#' }
#' @source <https://www.who.int/teams/global-influenza-programme/surveillance-and-monitoring/influenza-surveillance-outputs>
#' @srrstats {G5.1} Dataset `influenza` is exported and described in detail
#'   here. It is used in both example code and in testing.
#' @srrstats {G1.4} uses `Roxygen2` documentation
"influenza"
