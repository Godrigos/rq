#' @title
#' Calculates the Relative Quantity of a qPCR set
#'
#' @description
#' It takes a tibble of five columns: Sample, Target Name, Normalizer Name
#' and DCq (Delta Cq), to calculate the relative quantity and return another
#' tibble containing such data. Take a look at qpcr data available with this
#' package for and example.
#'
#' @param data A tibble
#' @param normalizer The name of the normalizer samples.
#' @param reference The name of the reference target.
#' @param treatment The name of the treatment target.
#' @param delimiter The delimiter character separating sample name and replica. Defaults to "-".
#' @return A tibble
#' @importFrom glue glue
#' @importFrom dplyr summarise group_by n
#' @importFrom tidyr separate_wider_delim pivot_longer
#' @importFrom cli cli_alert cli_alert_danger bg_red
#' @importFrom stats sd
#' @importFrom tidyselect all_of
#' @export
#' @examples
#' rq(qpcr, "pUC18", "WT", "ITS2")
rq <- function(data, normalizer, reference, treatment,
                    delimiter = '-') {
  Sample <- Cq <- DCq <- SD <- NULL
  if (is.null(normalizer)) {
    stop(bg_red('You must provide the normalizer samples name.'),
         call. = F)
  }
  if (is.null(reference)) {
    stop(bg_red('You must provide the reference target name.'),
         call. = F)
  }
  if (is.null(treatment)) {
    stop(bg_red('You must provide the treatment target name.'),
         call. = F)
  }
  # tidy the data to a format suitable for analyses
  tryCatch(
    {
      cat("\n")
      cli_alert("Tidying data...")
      qpcr <- pivot_longer(qpcr, cols = all_of(c(treatment, normalizer)),
                           names_to = 'Target', values_to = 'Cq')
    },
    error = function(cond) {
      cat("\n")
      cli_alert_danger(bg_red("Error tidying data."))
      stop(conditionMessage(cond), call. = F)
    }
  )
  # Sample means ignoring normalizer data
  cli_alert("Calculating table means...")
  qpcr_mean <- qpcr[qpcr$Target != normalizer,] |>
    group_by(Sample) |>
    summarise(
      `Mean Cq` = mean(`Cq`, na.rm = T),
      `Mean DCq` = mean(`DCq`, na.rm = T),
      .groups = "keep"
    )
  ## Sample Names Adjustment
  tryCatch(
    {
      cli_alert("Adjusting sample naming...")
      qpcr_mean <- qpcr_mean |>
        separate_wider_delim(
          Sample, delim = delimiter,
          names = c(
            "Sample", "Replica"
          )
        )
    },
    error = function(cond) {
      cat("\n")
      cli_alert_danger(bg_red("Error adjusting sample names."))
      stop(conditionMessage(cond), call. = F)
    }
  )
  ## Defining Base Value
  one <- mean(qpcr_mean$`Mean DCq`[
    qpcr_mean$Sample == reference
  ])
  if (is.na(one)) {
    cli_alert_danger(bg_red("Error with reference."))
    stop(reference, call. = F)
  }
  ## Relative Quantity
  cli_alert("Calculating Relative Quantity...")
  qpcr_mean$RQ <- 2^-(qpcr_mean$`Mean DCq` - (one))

  return(qpcr_mean)
}
