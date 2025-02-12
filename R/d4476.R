
#' Crack Detection
#'
#' Find possible cracking prior to structural failure in results from mechanical testing.
#'
#' @param data the data frame containing the load, deflection, and time data from test.
#' @param threshold default = 2, the drop in load to classify data point as a crack.
#' @param after default = 100, the load after which to begin looking for cracks.
#' @return A data frame of places where possible cracks occurred.
#' @export
crack.detect <- function(data, threshold = 2, after = 100) {

  data <- data %>%
      mutate(P2 = abs(lbf),
             Y2 = abs(in.),
             dlbf = c(NA,diff(P2)),
             din = c(NA,diff(Y2))
      ) %>%
      rename(t2=sec)

  cracks <- data %>%
    filter( (Y2 <= data[data$P2==max(data$P2),]$Y2) & (dlbf< -threshold) & (P2 >= after)) %>%
    select(-dlbf,-din,-in.,-lbf) %>%
    arrange(t2)

  if (nrow(cracks)==0) {
    cracks <- data[data$P2==max(data$P2),] %>%
      select(t2,P2,Y2)
  } else {
    cracks <- cracks %>%
      select(t2,P2,Y2) %>%
      first()
  }

  return(cracks)
}

#' D4476 Analysis
#'
#' Analyze data from an MTS data file to produce data frame of results required by D4476.
#'
#' @param data.file the input file path to the MTS data file.
#' @param meaurements the data frame containing the specimen name, Span Length, original diameter, and thickness measurements for each specimen.
#' @return A data frame of results for the test.
#' @export
d4476.analysis <- function(file_list,measurements,thresh = 2) {

  resultsdf <- merge(file_list,measurements,by='Sample') %>%
    mutate(R = D/2,
           gamma = Thickness/R,
           A = sqrt(gamma*(2-gamma)),
           B = 1-gamma,
           G = asin(A),
           H = 2*A*B,
           I = (R^4)*((1/4)*(G-A*B+2*(A^3)*B)-(4/9)*((A^6)/(G-A*B))),
           C = R*(1-((4*A^3)/(6*G-3*H)))
    ) %>%
    unnest(data) %>%
    mutate(lbf = abs(lbf),
           in. = abs(in.),
           inst.mod = ( (lbf*L^3)/(48*I*in.) )/1000) %>%
    nest(data = c(lbf,sec,in.,inst.mod)) %>%
    rowwise() %>%
    mutate(cracks = list(crack.detect(data,threshold = thresh))) %>%
    unnest(cracks) %>%
    rowwise() %>%
    mutate(pre = span_data(data = data,ref = P2)) %>%
    unnest(pre) %>%
    rename(P1 = lbf, t1 = sec, Y1 = in.)

  merge(resultsdf %>%
          unnest(data) %>%
          group_by(Sample) %>%
          mutate(E_b = mean(ifelse(sec >= t1 & sec <= t2,inst.mod,NA),na.rm=T)) %>%
          nest(data = c(lbf,sec,in.,inst.mod))  %>%
          na.omit(),
        get_failure(resultsdf)
  )%>%
    mutate(max.stress = (failload*L*C/(4*I))/1000,
           max.strain = 12*C*faildefl/(L^2)
    )
}

#' D4476 Plots
#'
#' Produce a list of plots for D4476 Testing with highlighted Modulus and Failure points.
#'
#' @param data.file the input file path to the MTS data file.
#' @param results the data frame containing the results from D4476 Analysis.
#' @return A ggplot of the test.
#' @export
d4476.plots <- function(results_df,invert = F) {
  plotlist <- load_defl_plot(results_df,invert)
  lapply(1:length(plotlist), function (i) add_d4476_plot_details(results_df,plotlist,i,invert))
}

#' Add D4476 Plot Details
#'
#' Adds failure point as well as a line between start and stop of modulus slope to a ggplot of load vs. deflection.
#'
#' @param results_df the data frame containing the results from D4476 Analysis.
#' @param plotlist the list of plots to add to.
#' @param invert default is False, when True a plot of deflection vs. load is returned.
#' @return A ggplot of the test with additional points for modulus and failure.
#' @export
add_d4476_plot_details <- function(results_df,plotlist,i,invert = F) {
  modulus_points <- results_df[i,] %>%
    select(Y1,P1,Y2,P2) %>%
    pivot_longer(cols = everything()) %>%
    mutate(record = right(name,1),var = left(name,1)) %>%
    select(-name) %>%
    pivot_wider(names_from = var,
                values_from = value)

  if (invert == F) {
  plotlist[[i]] +
    geom_line(data = modulus_points,aes(x=Y,y=P,color='Modulus'),size = 2) +
    geom_point(data = results_df[i,],aes(x=faildefl,y=failload,color='Failure'))
  } else {
    plotlist[[i]] +
      geom_line(data = modulus_points,aes(x=Y,y=P,color='Modulus'),size = 2) +
      geom_point(data = results_df[i,],aes(x=faildefl,y=failload,color='Failure'))+
      coord_flip()
  }
}
