---
title: Kalman filter vs smoother -- Missing data imputation in ctsem.
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

  ctsem is R software for statistical modelling using hierarchical state space models, of discrete or continuous time formulations, with possible non-linearities (ie state / time dependence) in the parameters. This is a super brief demo to show the basic intuition for Kalman filtering / smoother, and missing data imputation -- for a general quick start see https://cdriver.netlify.com/post/ctsem-quick-start/ , and for more details see the current manual at https://github.com/cdriveraus/ctsem/raw/master/vignettes/hierarchicalmanual.pdf




# Data
Lets load ctsem (if you haven't installed it see the quick start post!), generate some data from a simple discrete time model, and create some artificial missingness:


```r
set.seed(3)
library(ctsem)

y <- 6 #start value
n <- 40 #number of obs
for(i in 2:n){ 
  y[i] = .9 * y[i-1] + .1 + rnorm(1,0,.5) # ar * y + intercept + system noise
}
y=y+rnorm(n,0,.2) #add measurement error
y=data.frame(id=1,y=y,time=1:n) #create a data.frame for ctsem

missings <- c(5,17:19,26,35) #a few random observations...

ymissings=y #save the missings to check later
ymissings$y[-missings]<-NA

y$y[missings]<-NA #remove the selected obs from our analysis data set

head(y)
##   id        y time
## 1  1 6.158752    1
## 2  1 5.176335    2
## 3  1 4.408774    3
## 4  1 4.592951    4
## 5  1       NA    5
## 6  1 3.284191    6
```


# Model
If we're going to impute missing data, we need some kind of model for the imputation. The default model in ctsem is a first order auto / cross regressive style model, in either discrete or continuous time. There is correlated system noise, and uncorrelated measurement error. When multiple subjects are specified, the default is to have random (subject specific) initial states and measurement intercepts (with correlation between the two). For our purposes here, we're going to rely mostly on the defaults, and use a continuous time, differential equation format -- the discrete time form would also work fine for these purposes, but visualising the difference between the filter and smoother is much more obvious with the continuous time approach. Besides the defaults, we also fix the initial latent variance to a very low value (because we only have one subject), and the measurement error variance (because having system and measurement noise makes for better visuals but we can't easily estimate both with so little data). 

```r
model <- ctModel(type='stanct', # use 'standt' for a discrete time setup
  manifestNames='y',latentNames='ly',
  T0VAR=.01, #only one subject, must fix starting variance
  MANIFESTVAR=.2, #fixed measurement error, not easy to estimate with such limited data
  LAMBDA=1) #Factor loading fixed to 1

ctModelLatex(model) #requires latex install -- will prompt with instructions in any case
```


<img src="/post/2020-2-6-missingdata_files/figure-html/TEX-1.png" width="672" />

Fit using optimization and maximum likelihood (Possibly a couple of spurious warnings while estimating Hessian -- fixed on github):

```r
fit<- ctStanFit(y, model, optimize=TRUE, nopriors=TRUE, cores=2)
```

Then we can use summary and plotting functions. Note the differences between the first plot, using the Kalman *filter* predictions for each point, where the model simply extrapolates forwards in time and has to make sudden updates as new information arrives, and the smoothed estimates, which are conditional on *all* time points in the data -- past, present, and future. These plots are based on the maximum likelihood estimate / posterior mean of the parameters. 

```r
summary(fit)

kp<-ctKalman(fit,plot=TRUE, #predicted (conditioned on past time points) predictions.
  kalmanvec=c('y','yprior'),timestep=.1)
```

<img src="/post/2020-2-6-missingdata_files/figure-html/unnamed-chunk-4-1.png" width="672" />

```r

ks<-ctKalman(fit,plot=TRUE, #smoothed (conditioned on all time points) latent states.
  kalmanvec=c('y','ysmooth'),timestep=.1 )
```

<img src="/post/2020-2-6-missingdata_files/figure-html/unnamed-chunk-4-2.png" width="672" />

We can modify the ggplot objects we created to include the data we dropped:

```r
library(ggplot2)
kp= kp + geom_point(data = data.frame(Variable='y', Value=ymissings$y,
    Element='y',Time=ymissings$time), col='black')
plot(kp)
```

<img src="/post/2020-2-6-missingdata_files/figure-html/unnamed-chunk-5-1.png" width="672" />

```r

ks= ks + geom_point(data = data.frame(Variable='y', Value=ymissings$y,
    Element='y',Time=ymissings$time), col='black')
plot(ks)
```

<img src="/post/2020-2-6-missingdata_files/figure-html/unnamed-chunk-5-2.png" width="672" />

To access the imputed, smoothed values without plotting directly, we use the mean of our parameter samples to calculate various state and observation expectations using the ctStanKalman function. (In this case the samples are based on the Hessian at the max likelihood, but they could potentially come via importance sampling or Stan's dynamic HMC.) If we wanted more or different values to be available from such an approach, the easiest way is to include them as missing observations in the dataset used for fitting.

```r
k<-ctStanKalman(fit,collapsefunc = mean) 
str(k$ysmooth) 
k$ysmooth[1,missings,] 
```


<!-- We can also use the estimated mean and uncertainty to generate new data. For imputing purposes, we would generally want to account for uncertainty about the model parameters also, so we begin by creating mean and uncertainty estimates for random parameter vectors from our parameter distribution.  -->
<!-- ```{r } -->

<!-- imputed <- ctStanGenerateFromFit(fit,fullposterior=TRUE,imputemissings = TRUE,nsamples = 50) -->

<!-- matplot(imputed$generated$Y[,,1],type='l') -->

<!-- nsamples = 50 -->
<!-- k <- ctStanKalman(fit = fit, nsamples = nsamples) -->

<!-- matplot( t(k$ysmooth[,,1]), #selecting first manifest variable -->
<!--   type='l',lty=1,col=rgb(1,0,0,.3)) -->
<!-- ``` -->

<!-- As expected from such little data, there is quite some variability in the expected trajectory due to parameter uncertainty. We would also expect variability in the uncertainty (ysmoothcov) though this is not visualised. -->

<!-- For the next step, we randomly generate samples from our randomly drawn trajectories plus uncertainties. -->

<!-- ```{r } -->
<!-- newdat <-k$yprior -->
<!-- for(sampi in 1:nsamples){ -->
<!--   for(rowi in 1:length(k$time)){ -->
<!--     newdat[sampi,rowi,] <- k$yprior[sampi,rowi,] + #mean estimate -->
<!--       t(chol(k$ypriorcov[sampi,rowi,,])) %*% #plus Cholesky factor of cov matrix -->
<!--       rnorm(dim(newdat)[3]) #multiplied by some standard normal samples -->
<!--   } -->
<!-- } -->
<!-- matplot( t(k$ysmooth[,,1]), #plot expected trajectories -->
<!--   type='l',lty=1,col=rgb(1,0,0,.3)) -->

<!-- matplot(t(newdat[,,1]), #and sampled trajectories -->
<!--   type='p',pch=1,col=rgb(0,0,1,.3),add=TRUE) -->

<!-- ``` -->

<!-- With this, we could select just our missing time points, leaving us with an nsamples by nmissingobs by nvariables array: -->

<!-- ```{r } -->
<!-- imputed <- newdat[,missings,,drop=FALSE] -->
<!-- ``` -->

<!-- For a sanity check, lets compare to the original plots: -->

<!-- ```{r } -->
<!-- library(plyr) -->
<!-- newdf=adply(newdat,c(1,2,3)) #convert to data.frame -->
<!-- newdf[,2] = as.numeric(newdf[,2]) -->
<!-- colnames(newdf)[c(2,3,4)] <- c('Time','Element','Value') -->

<!-- ks= ks + geom_point(data = newdf[newdf$Time %in% missings,], col='black',alpha=.3) -->
<!-- plot(ks) -->

<!-- ``` -->
