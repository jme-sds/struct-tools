#' Right
#'
#' Return the last n characters from a string.
#'
#' @param x the string.
#' @param n the number of characters to keep.
#' @return String of last n characters of string x.
#' @export
right <- function(x,n) {
  substr(x,nchar(x)-n+1,nchar(x))
}


#' Left
#'
#' Return the first n characters from a string.
#'
#' @param x the string.
#' @param n the number of characters to keep.
#' @return String of first n characters of string x.
#' @export
left <- function(x,n) {
  substr(x,1,n)
}

right <- Vectorize(right)
left <- Vectorize(left)


#' Get Data Files
#'
#' Interactively select directory containing data files, returns a dataframe of sample names and dataframes containing raw data.
#'
#' @param No inputs.
#' @return A dataframe of sample names and dataframes containing raw data for each sample.
#' @export
get_data_files <- function(caption = "Select a directory") {
  if (exists('utils::choose.dir')) {
    dir <- choose.dir(caption = caption)
  } else {
    dir <- tk_choose.dir(caption = caption)
  }
  filepaths <- as.list(list.files(dir,recursive = T,pattern = '.txt',full.names = T))
  names <- lapply(filepaths,get_sample_name)
  data.frame(Sample = do.call(rbind,names),filepath = do.call(rbind,filepaths)) %>%
    rowwise() %>%
    mutate(data = list(read.table(filepath,skip = 7, header = T))) %>%
    select(-filepath)
}

#' Sort Plots
#'
#' Sort a list of ggplots by plot title.
#'
#' @param plots: a list of ggplots.
#' @return A list of ggplots sorted by plot title.
#' @export
sort_plots <- function(plots) {
  titles <- list()
  for (i in 1:length(plots)) {
    titles[i] <- plots[[i]]$labels$title
  }

  names(plots) <- titles
  sorted <- plots[order(names(plots))]
  return(sorted)
}

#' Get Sample Name
#'
#' Extract the sample name from an MTS data file.
#'
#' @param path: the path to the MTS data file from which to extract the sample name.
#' @return The sample name.
#' @export
get_sample_name <- function(path) {
  name <- readr::read_lines(path,n_max = 7)
  unlist(strsplit(name[3],': '))[2]
}

#' Get Failure Points for All Data
#'
#' Retrieve the maximum load as well as the deflection and time at maximum load for each specimen in a dataframe containing dataframes of specimen data.
#'
#' @param result_df: The data frame containing data frames of test data.
#' @return Data frame with failure points for each specimen.
#' @export
get_failure <- function(result_df) {
  result_df %>%
    unnest(data) %>%
    group_by(Sample) %>%
    slice_max(lbf,n=1) %>%
    select(lbf,in.,sec) %>%
    rename(failload = lbf,faildefl = in., failtime = sec)
}

#' Load vs. Deflection Scatter Plot.
#'
#' Produce a ggplot scatter plot of Load vs. Deflection for all data in a list of load deflection data sets.
#'
#' @param result_df: The data frame containing data frames of test data.
#' @param invert: Default is False. When true, a plot with Load on the x-axis and Deflection on y-axis is returned.
#' @return A list of ggplot objects.
#' @export
load_defl_plot <- function(results_df,invert = F) {

  names(results_df$data) <- results_df$Sample

  if (invert == F) {
    plots <- lapply(names(results_df$data),
                    function (k) ggplot(results_df$data[[k]], aes(x=in., y= lbf))+geom_point()+labs(title=k,x='Deflection [in]',y='Load [lbf]'))
  } else {
    plots <- lapply(names(results_df$data),
                    function (k) ggplot(results_df$data[[k]], aes(x=in., y= lbf))+geom_point()+labs(title=k,x='Deflection [in]',y='Load [lbf]')+coord_flip())
  }
  return(plots)
}

#' Span Data
#'
#' Finds the data point in a data set which is some percentage different from a given reference point.
#' e.g. if the reference point is 100 lbf, and the span is set to -5, the first point which is greater than or equal to 95% of the reference (95 lbf) is returned.
#'
#' @param data the data set, a data frame
#' @param ref the reference point, a numeric
#' @param span the amount you want to span in percent.
span_data <- function(data,ref,span = -1) {
  pct <- (100+span)/100

  data %>%
    filter(lbf>=pct*ref) %>%
    select(sec,lbf,in.) %>%
    first()
}
