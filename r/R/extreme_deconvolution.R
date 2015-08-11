.fixfix <- function(fix, ngauss) {
    if (is.null(fix)) {
        fix <- rep(0, ngauss)
    } else if (fix == FALSE | fix == TRUE) {
        fix <- rep(as.integer(fix), ngauss)
    } else if (length(fix) != ngauss) {
        warning("Dimension of fix* input does not match data (set all to the first entry)!")
        fix <- rep(as.integer(fix[0]), ngauss)
    } else {
        fix <- as.integer(fix)
    }
    return(fix)
}
extreme_deconvolution <- function(ydata,ycovar,
                          xamp,xmean,xcovar,
                          projection=NULL,
                          weight=NULL,
                          fixamp=NULL,fixmean=NULL,fixcovar=NULL,
                          tol=1.e-6,maxiter=1e9,w=0,logfile=NULL,
                          splitnmerge=0,maxsnm=FALSE,likeonly=FALSE,
                          logweight=FALSE)
{
    ndata <- dim(ydata)[1]
    dataDim <- dim(ydata)[2]
    ngauss <- length(xamp)
    gaussDim <- dim(xmean)[2]
    if (is.null(dim(ycovar)) | dim(ycovar)[2] == 1) {
        diagerrors <- TRUE
    } else {
        diagerrors <- FALSE
    }
    fixamp <- .fixfix(fixamp, ngauss)
    fixmean <- .fixfix(fixmean, ngauss)
    fixcovar <- .fixfix(fixcovar, ngauss)
    avgloglikedata <- 0
    #
    if(is.null(logfile)) {
        clog <- ''
        clog2 <- ''
        n_clog <- 0
        n_clog2 <- 0
    } else {
        clog <- paste(logfile, 'c.log', sep = "_")
        n_clog <- length(clog)
        clog2 <- paste(logfile, 'loglike.log', sep = "_")
        n_clog2 <- length(clog2)
    }
    #
    if(maxsnm)
        splitnmerge <- ngauss*(ngauss-1)*(ngauss-2)/2
    if(is.null(projection)) {
        noprojection <- TRUE
        projection <- list()
    } else {
        noprojection <- FALSE
    }
    if(is.null(weight)) {
        noweight <- TRUE
        logweights <- array(0)
    } else if (!logweight) {
        noweight <- FALSE
        logweights <- log(weight)
    } else {
        noweight <- FALSE
        logweights <- weight
    }
    #
    res <- .C("proj_gauss_mixtures_IDL",
       as.double(as.vector(t(ydata))),
       as.double(as.vector(t(ycovar))), 
			 as.double(as.vector(unlist(lapply(projection,t)))),
       as.double(as.vector(logweights)),
			 as.integer(ndata), as.integer(dataDim), 
			 xamp = as.double(as.vector(t(xamp))),
       xmean = as.double(as.vector(t(xmean))), 
			 xcovar = as.double(as.vector(unlist(lapply(xcovar,t)))),
       as.integer(gaussDim), as.integer(ngauss), 
			 as.integer(as.vector(fixamp)), as.integer(as.vector(fixmean)), 
			 as.integer(as.vector(fixcovar)), 
			 avgloglikedata = as.double(as.vector(avgloglikedata)),
       as.double(tol), as.integer(maxiter), as.integer(likeonly),
       as.double(w), as.character(clog), as.integer(n_clog), as.integer(splitnmerge),
			 as.character(clog2), as.integer(n_clog2),
			 as.integer(noprojection), as.integer(diagerrors),
			 as.integer(noweight),
       PACKAGE = "ExtremeDeconvolution"
       )
    #
    xmean <- matrix(res$xmean, dim(xmean), byrow = TRUE)
    start <- 1
    end <- 0
    for (i in 1:length(xcovar)) {
        end <- end + prod(dim(xcovar[[i]]))
        xcovar[[i]] <- matrix(res$xcovar[start:end], dim(xcovar[[i]]), byrow = TRUE)
        start <- end + 1
    }
    
    return(list(xmean=xmean, xamp=res$xamp, xcovar=xcovar, avgloglikedata=res$avgloglikedata))
}
