### Code for computing Density Recovery Profile (DRP).
### Wed 13 Feb 2002.

autodrp <- function(xs, ys, nbins, r, a=NULL) {
  ## Compute the autoDRP.  This is a simple wrapper around the cross DRP.
  crossdrp(xs, ys, xs, ys, nbins, r, a, auto=TRUE)
}

crossdrp <- function(xs1, ys1, xs2, ys2, nbins, r, a=NULL, auto=FALSE) {
  ## Compute the crossDRP.
  if (is.null(a)) {
    ## If a was not specified, we calculate bounds of dataset from the
    ## array.
    lims <- range(c(xs1, xs2));
    left <- lims[1]; right <-lims[2];

    lims <- range(c(ys1, ys2));
    bottom <- lims[1]; top <-lims[2];
    ##cat(paste(left, right, bottom, top, "\n"))
  }
  else {
    if (is.numeric(a)) {
      ## From the a vector, extract bounds of the area to study.
      left <- a[1]; right <- a[2]; bottom <- a[3]; top <- a[4];
    }
    else {
      stop("a is invalid - should be a vector of 4 numbers\n")
    }
  }
  
  l <- top-bottom; w <- right-left;

  ## now filter out points that are outside the bounding area A.
  subset1.x <- ((xs1 >= left) & (xs1 <= right))
  subset1.y <- ((ys1 >= bottom) & (ys1 <= top))
  subset1 <- which( subset1.x & subset1.y)
  xs1 <- xs1[subset1]; ys1 <- ys1[subset1]
  
  subset2.x <- ((xs2 >= left) & (xs2 <= right))
  subset2.y <- ((ys2 >= bottom) & (ys2 <= top))
  subset2 <- which( subset2.x & subset2.y)
  xs2 <- xs2[subset2]; ys2 <- ys2[subset2]

  if (auto) {
    ## For an autodrp, check subset1 and 2 are equal; this
    ## was not the case Mon 24 Feb 2003, see VC logs.
    stopifnot(identical(all.equal(subset1, subset2), TRUE))
  }

  ns <- binit2(xs1, ys1, xs2, ys2, nbins, r, auto)
  rs <- r * (0:(nbins-1))                 #starting radius of each annulus.
  area <- l * w
  npts1 <- length(xs1); npts2 <- length(xs2);

  ## For cross-correlation, need to use geometrical mean, rather than
  ## arithmetical mean.
  npts <- sqrt(npts1 * npts2)

  fs <- drpcorrections(r, nbins, l, w)

  density <- npts / area
  ## After checking Rodieck code, I now apply the correction factors to
  ## ns, immediately after counting, rather than applying them to areas.
  ## this seems to be the clinch!
  ns <- ns/ fs;
  areas   <- pi * r^2 * ( (2*(1:nbins))-1);# Areas of each annulus (eq 1).

  lambdas <- npts * density * areas     # Expected random count (eq 2).

  ds <- (ns / areas) / npts             # Density measures (eq 3)

  res <- drpeffrad(lambdas, ns,  npts, density)
  effrad <- res$r
  if (effrad < 0) {
    ## have a back-up in case we can't compute effective radius.
    effrad <- 1; p <- 0; maxr <- 0; k <- 0; Dc <- 0;
  } else {
    res <- drppackingfactor(effrad, density);
    p <- res$p; maxr <- res$maxr
    res <- drpreliability(density, area, r);
    k <- res$k; Dc <- res$Dc;
  }

  sd.poisson <- Dc / (sqrt( (2*(1:nbins)) - 1))
  mid.bin <- seq(from=r/2, by=r, length=nbins)
  res <- list(effrad = effrad, p = p, maxr = maxr, k = k, Dc = Dc,
              ds =ds, density=density, n1=npts1, n2=npts2,
              nbins=nbins,r=r,
              ns=ns, fs=fs,
              areas=areas,
              sd.poisson=sd.poisson,
              mid.bin=mid.bin)
  class(res) <- "sjedrp"

  res
}

plot.sjedrp <- function (x, scale=1, title=NULL, mirror=FALSE,
                         show.title=TRUE, ylab='density',
                         xlab='distance') {
  ## Plot the results of the density recovery profile.
  ## The SCALE parameter allows us to change the scale of the y axis.
  ## e.g. when going from um^2 to mm^2 use a scale of 1e6.
  ## MIRROR allows the plot to be mirrored across the y-axis.
  ## SHOW.TITLE:  If false,  do not add any title.
  hts <- (x$ds*scale)
  last.bin <- x$nbins * x$r
  plot.label <- paste(title)
  if (!show.title)
    plot.label <- NULL
  
  ##names(hts) <- x$rs

  if (mirror) {
    ## Make the symmetrical plot.
    hts <- c(rev(hts), hts)
    
    barplot(hts, col="gray",space=0, width=x$r, xlim=c(0,2*last.bin),
            main=plot.label, xlab=xlab, ylab=ylab)

    ## mean density line:
    lines( c(0,2*last.bin), c(x$density, x$density)*scale)

    ## draw x-axis
    axis(1, at=c(0, last.bin, 2*last.bin),
         labels=c(-last.bin, 0, last.bin))

    ## Draw the effective radius.
    lines( c(x$effrad, x$effrad) + last.bin, scale * c(0, x$density))
    
    points( c(x$maxr + last.bin), c(0),pch='|')    #maximum radius.
  } else {
    ## standard, one-sided plot.
    barplot(hts, col="gray",space=0, width=x$r, xlim=c(0,last.bin),
            ylab=ylab, xlab=xlab,
            main=plot.label)
    lines( c(0,last.bin), c(x$density, x$density)*scale)
    axis(1, at=c(0, last.bin))
    lines( c(x$effrad, x$effrad), scale * c(0, x$density))
    points( c(x$maxr), c(0),pch='|')    #maximum radius
  }
}


drppackingfactor <- function (effrad, d) {
  ## Return the packing factor P and the maximum radius MAX_R.
  maxr <-  sqrt(sqrt(4/3)/d)          # (eq 19)
  p <- (effrad/maxr) ^2               # (eq 21)
  if ( (p<0.0) ||(p>1.0))
    warning(paste("packing factor", p,"not in range [0,1]"))
  
  list (p = p, maxr = maxr)
}

drpeffrad <- function (lambdas, ns, n, d) {
  ## Return the effective radius.
  ## If the effective radius cannot be computed, it returns a negative value.
  diffs <- lambdas - ns;
  fracs <-  ns / lambdas;
  
  ##negs <- which(diffs < 0);
  
  ## Rather than finding first value when lambda - n < 0, we instead see
  ## when the value n/lambda > 1.  This should be equivalent, but
  ## sometims rounding errors mean that it never finds such a point.  In
  ## that case, we take the largest value of n/lambda and assume that
  ## bin is the one where the histogram first crosses the density line.

  ## In the 2001+ version of Rodieck's program, his PDF notes that
  ## this problem can occasionly occur in his program too.
  negs <-  which(fracs > 1.0);
  if (length(negs) == 0) {
    v <- max(fracs); firstnegative <- which(fracs == v)
    warning(paste("no negative elements in drpeffrad. Closest"
                  ,v, "at position", firstnegative))
    volume  <- sum( diffs[1:(firstnegative-1)])/n; # (eq 7)
    r <- sqrt(volume/(pi*d));	# (eq 9)
  }
  else {
    ## We have negative values, so we can just sum the difference
    ## values until just before the first negative entry.
    firstnegative <- negs[1];
    if (firstnegative == 1) {
      ## If the firstnegative element is in bin 1, this means that bin1
      ## is above average -- this can happen (and perhaps is indicative of
      ## clustering), in which case I think volume, and hence effective radius
      ## should be zero.
      ##warning("firstnegative should be greater than 1")
      volume <- 0
    } else {
      volume  <- sum( diffs[1:(firstnegative-1)])/n; # (eq 7)
    }
    r <- sqrt(volume/(pi*d));	# (eq 9)
  }
  
  list( r=r, first = firstnegative)
}

drpreliability <- function (density, area, r) {
  ## Return reliability factor K and critical density Dc.
  Dc <- 1.0 / (sqrt( area* pi) * r);    # (eq 11)
  k  <- density/ Dc;                    # (eq 13)
  list (k = k, Dc = Dc)
}

binit2 <- function (xs1, ys1, xs2, ys2, nbins, r, auto) {
  ## Bin the cross-correlation distances.  This is a wrapper around
  ## a C routine that does all the hard work.
  ## The bins are made thus:
  ## [0, r), [r 2r), [2r 3r), ... [(n-1)r r)
  ## which matches Rodieck.
  npts1 <- length(xs1)
  npts2 <- length(xs2)
  z <- .C(C_drp_bin_it_r,
          as.double(xs1), as.double(ys1), as.integer(npts1),
          as.double(xs2), as.double(ys2), as.integer(npts2),
          as.integer(nbins), as.double(r), as.integer(auto),
          ## create memory to store return values.
          ns = integer(nbins) )
  z$ns
}

binit2.check <- function() {
  ## Check whether binit2 is working by using outer-product to calculate
  ## all pairs of distances.  This could be wasteful to use in general.
  n <- 50
  nbins <- 20
  r <- 10
  xs <- 500 * runif(n); ys <- 500 *runif(n)
  b <- binit2(xs, ys, xs, ys, nbins, r, auto=1)

  dist <- function(a, b) {
    ## compute distance between cells A and B.
    dx <- xs[a] - xs[b];
    dy <- ys[a] - ys[b];
    dist <- sqrt((dx**2)+(dy**2))
  }
  
  dists2 <- outer(1:n, 1:n, dist)
  ## Only difference should be that first bin has n pts in it, since
  ## the outer product method includes self-counts.
  h2 <- hist(dists2, breaks=seq(from=0, by=r, to=max(dists2)+r),right=F,plot=F)
  h2$counts[1:nbins] - b
}

drpcorrections <- function(r, nbins, l, w) {       
  ## Return the correction factors.  Note correction factors are applied
  ## to middle of bin, for consistency with Rodieck.  Rodieck had 1/pi=0.312
  ## whereas that should be 0.318, but don't think that was enough alone
  ## to make a difference.
  ## (eq 32)
  rs <- ( (0:(nbins-1)) + 0.5)* r
  fs <- 1 - (2*rs *(l+w)/( pi * l * w)) + ( (rs * rs)/ (pi * l * w))
  fs
}

drp.makestf <- function (x, y, file) {
  ## Make an STF file suitable for reading in with Bob Rodieck's program.
  npts <- length(x)
  stopifnot(npts == length(y))
  type <- numeric(length=npts)+1
  z <- numeric(length=npts)
  d <- data.frame(type=type, x=x, y=y, z=z, tag=z, br=z)
  write.table(d, file=file,
              eol="\r",                 #\r for the macintosh...
              sep="\t", quote=FALSE, row.names=FALSE)
}

drp.makemac <- function(file, bogus=T) {
  ## Tue 25 Feb 2003
  ## Convert  FILE into format suitable for reading into Mac.
  ## See ../data/safe for my output from mac.
  ## Assumes top line has x,y on it.
  data <- read.table(file, header=T)
  x <- data$x
  y <- data$y
  x.range <- range(x)
  y.range <- range(y)
  delta <- (x.range[2] - x.range[1])*0.1

  ## add two data points to extend the lower and upper range.
  
  lo <- c(x.range[1], y.range[1]) - delta
  hi <- c(x.range[2], y.range[2]) + delta
  data2 <- rbind(lo, hi, cbind(x,y))
  newfile <- gsub(".txt", ".macdrp.txt", file)
  write.table(data2, newfile, quote=F,
              eol="\r",                 #get MAC line endings.
              row.names=F, col.names=F)
}


## myauto <- function(xs, ys, nbins=20, dr=10) {
##   ## This was a simple implementation, following closely the
##   ## Rodieck implemenation.  Think I can comment this out for now.
##   bins <- binit2(xs, ys, xs, ys, nbins, dr, TRUE)
##   counts <- bins
##   bw <- dr 
##   mPoints <- length(xs)
##   scale <- 1.0 / (pi * mPoints * bw * bw)
##   k <- scale / (2*(1:nbins) - 1)
##   bins <- bins * k

##   ## now correct for bounds.
##   width <- 700; height <- 700
##   c1 <- 2/pi * ((1/width) + (1/height))
##   c2 <- 0.312 / (width*height)          #should be 0.318
##   corr <- real(nbins)
##   for (i in 1:nbins) {
##     r1 <- ( (i-0.5)*dr)                 #add correction to middle of bin width
##     a <- 1 - (c1*r1) + (c2*r1*r1)
##     corr[i] <- 1/a
##     bins[i] <- (bins[i] * corr[i])

##   }

##   barplot(bins, col="gray")
##   density <- mPoints/(width*height) ##/ (micromilli*micromilli)
##   abline(h=density)

##   ## try to get the effect radius.
##   i <- 1
##   deadVolume <- 0.0
##   while( bins[i] < density) {
##     deadVolume <- deadVolume + ( (density - bins[i])* ((2*i)-1))
##     i <- 1 + i
##   }
##   effrad <- bw * sqrt( deadVolume/density)
##   print(effrad)

##   list( bins=bins, corr=corr, counts=counts)
## }



autocorr <- function(xs, ys, nbins, r, a=NULL) {
  ## Compute the autoDRP.  This is a simple wrapper around the cross DRP.
  crosscorr(xs, ys, xs, ys, nbins, r, a, auto=TRUE)
}

crosscorr <- function(xs1, ys1, xs2, ys2, nbins, r, a=NULL, auto=FALSE) {
  ## Compute the crossDRP.
  if (is.null(a)) {
    ## If a was not specified, we calculate bounds of dataset from the
    ## array.
    lims <- range(c(xs1, xs2));
    left <- lims[1]; right <-lims[2];

    lims <- range(c(ys1, ys2));
    bottom <- lims[1]; top <-lims[2];
    ##cat(paste(left, right, bottom, top, "\n"))
  }
  else {
    if (is.numeric(a)) {
      ## From the a vector, extract bounds of the area to study.
      left <- a[1]; right <- a[2]; bottom <- a[3]; top <- a[4];
    }
    else {
      stop("a is invalid - should be a vector of 4 numbers\n")
    }
  }
  
  l <- top-bottom; w <- right-left;

  ## now filter out points that are outside the bounding area A.
  subset1.x <- ((xs1 >= left) & (xs1 <= right))
  subset1.y <- ((ys1 >= bottom) & (ys1 <= top))
  subset1 <- which( subset1.x & subset1.y)
  xs1 <- xs1[subset1]; ys1 <- ys1[subset1]
  
  subset2.x <- ((xs2 >= left) & (xs2 <= right))
  subset2.y <- ((ys2 >= bottom) & (ys2 <= top))
  subset2 <- which( subset2.x & subset2.y)
  xs2 <- xs2[subset2]; ys2 <- ys2[subset2]

  if (auto) {
    ## For an autodrp, check subset1 and 2 are equal; this
    ## was not the case Mon 24 Feb 2003, see VC logs.
    stopifnot(identical(all.equal(subset1, subset2), TRUE))
  }

  npts1 <- length(xs1)
  npts2 <- length(xs2)
  maxn <- npts1*npts2                   #CONSERVATIVE OVERESTIMATE!!!
  z <- .C(C_cross_corr_r,
          as.double(xs1), as.double(ys1), as.integer(npts1),
          as.double(xs2), as.double(ys2), as.integer(npts2),
          as.integer(nbins), as.double(r), as.integer(auto),
          ## create memory to store return values.
          dx = double(maxn), dy=double(maxn),k=integer(1)
          )

  dx <- z$dx[1:z$k]
  dy <- z$dy[1:z$k]
  res <- list(x = dx, y = dy, n1=npts1, n2=npts2,
              nbins=nbins,r=r)
  class(res) <- "sjecorr"

  res
}

plot.sjecorr <- function(x, pts.cex=0.5, ...) {
  nbins <- x$nbins
  r <- x$r
  maxr <- nbins * r
  plot(NA, xlim=c(-maxr, maxr), ylim=c(-maxr, maxr), yaxt='n',
       xlab='', ylab='', asp=1, bty='n')
  points(x$x, x$y, pch=19, cex=pts.cex)
  symbols(x=rep(0,nbins), y=rep(0,nbins),
          circles=(1:nbins)*r, inch=F, add=T)
}
