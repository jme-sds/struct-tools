
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

  if (!"dlbf" %in% colnames(data)) {
    data <- data %>%
      mutate(lbf = abs(lbf),
             in. = abs(in.),
             inst.mod = ( (lbf*L^3)/(48*I*in.) )/1000,
             dlbf = c(NA,diff(lbf)),
             din = c(NA,diff(in.))
      )
  }

  data %>%
    filter( (in. <= data[data$lbf==max(data$lbf),]$in.) & (dlbf< -threshold) & (lbf >= after)) %>%
    arrange(sec)

}


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


#' D4476 Analysis
#'
#' Analyze data from an MTS data file to produce data frame of results required by D4476.
#'
#' @param data.file the input file path to the MTS data file.
#' @param meaurements the data frame containing the specimen name, Span Length, original diameter, and thickness measurements for each specimen.
#' @return A data frame of results for the test.
#' @export
d4476.data <- function(data.file,measurements,thresh = 2) {
  name <- readr::read_lines(data.file,n_max = 7)
  sample.name <- unlist(strsplit(name[3],': '))[2]

  data <- read.table(data.file,skip = 7, header = T)

  details <- measurements[measurements$sample==sample.name,]
  L <- details$L
  D <- details$D
  t <- details$t

  R <- D/2
  gam <- t/R
  A <- sqrt(gam*(2-gam))
  B <- 1-gam
  G <- asin(A)
  H <- 2*A*B
  I <- (R^4)*((1/4)*(G-A*B+2*(A^3)*B)-(4/9)*((A^6)/(G-A*B)))
  C <- R*(1-((4*A^3)/(6*G-3*H)))

  cat(I)

  data <- data %>%
    mutate(lbf = abs(lbf),
           in. = abs(in.),
           inst.mod = ( (lbf*L^3)/(48*I*in.) )/1000
    )

  crack <- crack.detect(data,thresh)[1]

  t1 <- data[data$lbf>=0.99*crack$lbf,]$sec[1]
  t2 <- crack$sec

  P1 <- data[data$lbf>=0.99*crack$lbf,]$lbf[1]
  P2 <- crack$lbf

  Y1 <- data[data$lbf>=0.99*crack$lbf,]$in.[1]
  Y2 <- crack$in.

  E_b <- data %>%
    filter( (sec>=t1) & (sec<=t2) ) %>%
    summarize(E_b = mean(inst.mod))

  failload = max(data$lbf)
  faildefl = data[data$lbf==max(data$lbf),1]

  max.stress <- (failload*L*C/(4*I))/1000
  max.strain <- 12*C*faildefl/(L^2)



  data.frame(Sample=sample.name,
             L = L,
             D = D,
             Thickness = t,
             FailLoad=max(data$lbf),
             FailDefl=data[data$lbf==max(data$lbf),1],
             MaxStress=max.stress,
             Modulus=E_b[[1]],
             MaxStrain=max.strain,
             R = R,
             gamma = gam,
             A = A,
             B = B,
             G = G,
             H = H,
             I = I,
             C = C
  )

}


#' D4476 Plots
#'
#' Produce plots of D4476 Testing with highlighted Modulus and Failure points.
#'
#' @param data.file the input file path to the MTS data file.
#' @param results the data frame containing the results from D4476 Analysis.
#' @return A ggplot of the test.
#' @export
d4476plots <- function(data.file,results) {
  name <- read_lines(data.file,n_max = 7)
  sample.name <- unlist(strsplit(name[3],': '))[2]
  data <- read.table(data.file,skip = 7, header = T)
  result <- results[results$Sample==sample.name,]

  dots <- data.frame(row.names = c('Modulus','Failure'),
                     lbf = c(result$P,result$FailLoad),
                     in. = c(result$Y,result$FailDefl)
  )

  data <- data %>%
    mutate(lbf = abs(lbf),
           in. = abs(in.)
    )

  ggplot(data=data,aes(y=in.,x=lbf))+
    geom_point()+
    labs(y='Displacement [in.]',x='Load [lbf]',title = sample.name)+
    geom_point(data=dots,aes(x=in.,y=lbf),color='red')+
    geom_shadowtext(data=dots,label = rownames(dots),nudge_y = 15,check_overlap = T)+
    theme(plot.title = element_text(hjust = 0.5))
}

