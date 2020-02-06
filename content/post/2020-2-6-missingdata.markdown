---
title: Missing data imputation in ctsem -- Kalman filter / smoother.
author: Charles Driver
date: '2020-2-6'
slug: missingdata
categories: []
tags: []
subtitle: ''
summary: ''
authors: []
lastmod: '2019-12-19T18:33:24+01:00'
featured: no
image:
  caption: ''
focal_point: ''
preview_only: no
projects: []
---

  ctsem is R software for statistical modelling using hierarchical state space models, of discrete or continuous time formulations, with possible non-linearities in the parameters. This is a super brief demo of missing data imputation using a Kalman smoother -- for a more complete quick start see https://cdriver.netlify.com/post/ctsem-quick-start/ , and for even more details see the current manual at https://github.com/cdriveraus/ctsem/raw/master/vignettes/hierarchicalmanual.pdf




# Data
Lets load ctsem (if you haven't installed it see the quick start post!) and pull in some data:


```r
library(ctsem)
ssdat <- data.frame(id=1,
  time=do.call(seq,as.list(attributes(sunspot.year)$tsp))[1:40],
  ss=sunspot.year[1:40])

missings <- c(6,26,35) #a few random observations...

ssdatmissings=ssdat #save the missings to check later
ssdatmissings$ss[-missings]<-NA

ssdat$ss[missings]<-NA #remove the selected obs from our analysis data set

head(ssdat)
##   id time ss
## 1  1 1700  5
## 2  1 1701 11
## 3  1 1702 16
## 4  1 1703 23
## 5  1 1704 36
## 6  1 1705 NA
```


# Model
If we're going to impute missing data, we need some kind of model for the imputation. The default model in ctsem is a first order auto / cross regressive style model, in either discrete or continuous time. There is correlated system noise, and measurement error. When multiple subjects are specified, the default is to have random (subject specific) initial states and measurement intercepts (with correlation between the two). For our purposes here, we're going to rely on the defaults, and use a discrete time, difference equation format. 

```r
ssmodel <- ctModel(type='stanct',
  manifestNames='ss',
  latentNames=c('lss'),
  T0VAR=1e-3, #only one subject, must fix starting variance
  LAMBDA=1)

ctModelLatex(ssmodel) 
```


<img src="/post/2020-2-6-missingdata_files/figure-html/TEX-1.png" width="672" />

Fit using optimization and maximum likelihood:

```r
ssfit<- ctStanFit(ssdat, ssmodel, optimize=TRUE, nopriors=TRUE, cores=2)
```

Then we can use summary and plotting functions:

```r
summary(ssfit)

kp<-ctKalman(ssfit,plot=TRUE, #predicted (conditioned on past time points) predictions.
  kalmanvec=c('y','yprior'),timestep=.01)
```

<img src="/post/2020-2-6-missingdata_files/figure-html/unnamed-chunk-4-1.png" width="672" />

```r

ks<-ctKalman(ssfit,plot=TRUE, #smoothed (conditioned on all time points) latent states.
  kalmanvec=c('y','ysmooth'),timestep=.01 )
```

<img src="/post/2020-2-6-missingdata_files/figure-html/unnamed-chunk-4-2.png" width="672" />


```r
library(ggplot2)
kp= kp + geom_point(data = data.frame(Variable='ss', Value=ssdatmissings$ss,
    Element='y',Time=ssdatmissings$time), col='black')
plot(kp)
```

<img src="/post/2020-2-6-missingdata_files/figure-html/unnamed-chunk-5-1.png" width="672" />

```r

ks= ks + geom_point(data = data.frame(Variable='ss', Value=ssdatmissings$ss,
    Element='y',Time=ssdatmissings$time), col='black')
plot(ks)
```

<img src="/post/2020-2-6-missingdata_files/figure-html/unnamed-chunk-5-2.png" width="672" />

