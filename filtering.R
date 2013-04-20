source("~/PWLisp/utilities.R")
dyn.load(paste("~/RLibs/filtering", .Platform$dynlib.ext, sep=""))
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

filterTimeSeriesDirectR <- function(timeSeries, theFilter) {
    nFilter <- length(theFilter)
    nFilterM1 <- nFilter -1
    nOutput <- length(timeSeries) - nFilterM1

    theFilter <- rev(theFilter)
    result <- array(NA, nOutput)
    
    for(i in 1:nOutput) {
        result[i] <- timeSeries[i] * theFilter[1]
        for(j in 1:nFilterM1) {
            result[i] <- result[i] + timeSeries[i+j] * theFilter[j +1]
        }
    }
    return(list(result=result, nOutput=nOutput))
}

filterTimeSeriesDirect <- function(timeSeries, theFilter) {

    nFilter <- length(theFilter)
    nTimeSeries <- length(timeSeries)
    nOutput <- nTimeSeries - nFilter + 1

    out <- .Fortran("filterTimeSeriesDirect", as.double(timeSeries),
                    as.double(theFilter), as.integer(nTimeSeries),
                    as.integer(nFilter), nOutput=as.integer(nOutput),
                    result=double(nOutput))
                                              
    return(list(result=out$result, nOutput=out$nOutput))
}

        



##filterdirect
##filter( c(1,2,3,4,5,-5,-7,-9, 10, 11, 13), c(1,-1))
##filterfft
##convolve( c(1,2,3,4,5,-5,-7,-9), c(1,-1), type="filter")

##perhaps try something like...
##Re(fft(Conj(fft(tspad)*fft(fpad)))/256)[3:(length(ts))]

filterWfft_old<- function(timeSeries, filter) {
    ndata <- length(timeSeries)
    nfilter <- length(filter)
    nfft <- 2^(ceiling(log2(ndata)))
    tsPad <- c(timeSeries, rep(0.0, nfft-ndata))
    fPad <- c(filter, rep(0.0, nfft-nfilter))
    return((Re(fft(Conj(fft(tsPad)*fft(fPad))))/nfft)[nfilter:ndata])
}    



filterWfftR<- function(timeSeries, filter1) {
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
    nblocks <- floor(ndata / blockLen) +
        if(ndata %% blockLen !=0) 1 else 0
    ##npad <- nblocks*blockLen - nfft
    nResult <- ndata - nfilter + 1
    result <- array(NA, nResult)
    fPad <- fft(fPad)
    i <- 0

    while(i*blockLen + nfft < ndata) {
        tsfft <- fft(timeSeries[(i*blockLen +1) : (i*blockLen + nfft)])
        result[(i*blockLen +1) : (i*blockLen + blockLen)] <-
            (Re(fft(tsfft*fPad, inverse=T)/nfft))[nfilter:nfft]
        i <- i +1
    }
    
    lenLeft <- ndata - i*blockLen  
    if(lenLeft != 0) {
        
        tsfft <- fft(c(timeSeries[(i*blockLen +1) : ndata],
                       array(0.0, nfft - lenLeft)))
        result[(i*blockLen +1) : nResult] <-
            (Re(fft(tsfft*fPad,
                    inverse=T)/nfft))[nfilter:(nfilter + lenLeft - nfilter)]
    }   

    return(list(result=result, nOutput=nResult))
}    

filterWfft<- function(timeSeries, filter1) {
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

##filterWfft(  c(1,2,3,4,5,-5,-7,-9), c(1,-1))
