## set local paths

library("PWLISPR")
pathToLibs <- "~/RLibs/"

##library("multitaper")
## load dynamic library
filteringLib <- paste(pathToLibs, "filtering", .Platform$dynlib.ext, sep="")
dyn.load(filteringLib)
useFortranLib <- TRUE

## set to false if no lib compiled



## ;;; Note: filter-time-series-direct is primarily intended to be used
## ;;; with short filters for which use of fft's would not be effecient.
##   "given
##    [1] time-series (required)
##        ==> a vector containing a time series
##            x_0, x_1, ..., x_{N-1}
##    [2] the-filter (required)
##        ==> a vector containing the filter coefficients
##            g_0, g_1, ..., x_{K-1}
##    [3] start (keyword; 0)
##        ==> start index of time-series to be used
##    [4] end (keyword; length of time-series)
##        ==> 1 + end index of time-series to be used
##    [5] result (keyword; vector of appropriate length)
##        <== vector to contain filtered time series
##                  K-1
##            y_t = SUM g_k x_{t+K-1-k},  t = 0, ..., N-K+1
##                  k=0
##  returns
##    [1] result, a vector containing the filtered time series
##    [2] the number of values in the filtered time series
## ---
## Note: result can be the same as time-series"

filterTimeSeriesDirectF <- function(timeSeries, theFilter) {

    nFilter <- length(theFilter)
    nTimeSeries <- length(timeSeries)
    nOutput <- nTimeSeries - nFilter + 1

    out <- .Fortran("filterTimeSeriesDirect", as.double(timeSeries),
                    as.double(theFilter), as.integer(nTimeSeries),
                    as.integer(nFilter), nOutput=as.integer(nOutput),
                    result=double(nOutput))
                                              
    return(list(result=out$result, nOutput=out$nOutput))
}

        



filterWfftF <- function(timeSeries, filter1) {
    ##           K-1
    ##     y_t = Sum g_k x_{t+K-1-k},  t = 0, ..., N-K+1
    ##           k=0
    ## http://en.wikipedia.org/wiki/Cyclic_convolution
    
    ndata <- length(timeSeries)
    nfilter <- length(filter1)
    ## 4 * next power of 2
    nfft <- 4*2^(ceiling(log2(nfilter)))
    fPad <- c(filter1, rep(0.0, nfft-nfilter))
    blockLen <- nfft - nfilter +1
    ##nblocks <- floor(ndata / blockLen) +
        ##if(ndata %% blockLen !=0) 1 else 0
    ##npad <- nblocks*blockLen - nfft
    nResult <- ndata - nfilter + 1
    result <- array(NA, nResult)

    out <- .Fortran("filterWfft", as.double(timeSeries), as.double(fPad),
                    as.integer(ndata), as.integer(nfilter), as.integer(nfft),
                    as.integer(nResult), as.integer(blockLen),
                    result=double(nResult), complex(nfft), double(nfft),
                    complex(nfft))

    return(list(result=out$result, nOutput=nResult))
}
