Generalised Ensemble Kalman Inversion - Multivariate Normal Model
================
Owen Jackson
2025-01-15

# Results For Multivariate Normal Model

## Import Functions

### External Packages

``` r
pacman::p_load(pacman, purrr, mvtnorm, mcmc, MASS, glue)
```

### Algorithms

``` r
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/eki.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/mcmc_normal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/mcmc_lognormal.R')
```

### Models

``` r
# MVN Normal
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/mvn_normal/eki_normal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/mvn_normal/eki_normal_known_var.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/mvn_normal/eki_normal_known_mean.R')
# MVN Lognormal
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/mvn_lognormal/eki_lognormal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/mvn_lognormal/eki_lognormal_var_only.R')
```

### Sampling and PDFs

``` r
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/pdfs/pdfs_normal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/pdfs/pdfs_lognormal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples/samples_normal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples/samples_lognormal.R')
```

### Utils

``` r
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/utils/tempering.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/utils/eki_helper.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/results/plots_normal.R')
```

``` r
generate_results <- function(num_particles, true_parameters) {
  
  true_data <- likelihood_normal(true_parameters)
  
  # Model 1
  eki_result_known_var <- eki_normal_known_var(num_particles, true_data, true_parameters, 
                                               prior_params, adaptive = adaptive)
  plot_eki_normal_known_var(eki_result_known_var, true_data, 
                            true_parameters, prior_params)
  
  # Model 2
  eki_result_known_mean <- eki_normal_known_mean(num_particles, true_data, true_parameters, 
                                                 prior_params, adaptive = adaptive)
  plot_eki_normal_known_mean(eki_result_known_mean, true_parameters, prior_params)
  
  # Model 3
  eki_result <- eki_normal(num_particles, true_data, true_parameters, 
                           prior_params, adaptive = adaptive)
  plot_eki_normal(eki_result, true_parameters, prior_params)
  
  # MCMC Comparison for Model 3
  run.1 <- normal_mcmc(true_data, true_parameters, prior_params)
  chain.1 <- run.1$batch
  run.2 <- normal_mcmc(true_data, true_parameters, prior_params)
  chain.2 <- run.2$batch
  plot_mcmc_trace_plots(chain.1, chain.2, burnin = 500)
  plot_mcmc_histogram(chain.1, chain.2, true_parameters, prior_params, 
                      burnin = 500)
  
  plot_overlay(eki_result, chain.1, chain.2, true_parameters, prior_params)
  
}
```

``` r
generate_lognormal_results <- function(num_particles, true_parameters, prior_params) {
  
  true_data <- likelihood_lognormal(true_parameters)

  # Model 2
  eki_result_known_mean <- eki_lognormal_var_only(num_particles, true_data, 
                                                 true_parameters,
                                                 prior_params, adaptive = adaptive)
  plot_eki_normal_known_mean(eki_result_known_mean, true_parameters, prior_params)
  
  # Model 3
  eki_result <- eki_lognormal(num_particles, true_data, true_parameters, 
                           prior_params, adaptive = adaptive)
  plot_eki_normal(eki_result, true_parameters, prior_params)
  
  # Correlation between particles
  alpha_particles <- eki_result$particles[, 1]
  sigma2_particles <- exp(eki_result$particles[, 2])
  
  plot(alpha_particles, sigma2_particles)
  
  # MCMC Comparison for Model 3
  run.1 <- lognormal_mcmc(true_data, true_parameters, prior_params)
  chain.1 <- run.1$batch
  run.2 <- lognormal_mcmc(true_data, true_parameters, prior_params)
  chain.2 <- run.2$batch
  plot_mcmc_trace_plots(chain.1, chain.2, burnin = 500)
  plot_mcmc_histogram(chain.1, chain.2, true_parameters, prior_params,
                      burnin = 500)
  
}
```

## The Model

We assume that we have access to a single observation from a
multivariate normal distribution $y \sim N(\alpha x, \sigma^2I)$ where
$\alpha$ is a scalar and $x$ is a known vector.

I look at three variations of this model

- The variance is known, but the mean is not $\theta = \alpha$

- The mean is known, but the variance is not $\theta = \sigma^2$

- Both parameters are unknown $\theta = (\alpha, \sigma^2)$.

I look at a range of combinations of $\alpha, x$ and $\sigma^2$,
comparing against MCMC sampling. Of particular interest is how well GEKI
is able to estimate the noise parameter, since previous EKI algorithms
require the noise parameter to be known.

### Priors

In all combinations we draw from the same priors (if that parameter is
considered unknown):

- $\alpha \sim ~ N(0, 5^2)$

- $\log(\sigma^2) \sim N(2, 1^2)$

``` r
prior_params <- list(alpha.mean = 0, alpha.sd = 5, 
                     sigma2.mean = 2, sigma2.sd = 1)
```

### Adaptive Tempering

``` r
adaptive = TRUE
```

### Output

For each variation of the model, I plot the EKI histograms of the
particles against their prior. In the known mean, unknown variance case,
I also plot the posterior since there is an analytical posterior in this
special case. This serves as a check that the algorithm has been
implemented correctly.

## Changing Number of Dimensions

In this experiment I keep $\alpha$ and $\sigma^2$ fixed at 2 and 4
respectively and let $x$ be a vector of 1s corresponding to 10, 50 and
100 dimensions. I keep the number of particles fixed at 400. The purpose
of this is to see how the EKI algorithm performs when the number of
dimensions increases.

``` r
num_dimensions <- c(10, 50, 100)
num_particles <- 400
for (dimension in num_dimensions) {
  true_parameters <- list(alpha = 2, sigma = 2, x = rep(1, dimension))
  generate_results(num_particles, true_parameters)
}
```

    ## Next temp is 0.08687011232338
    ## Next temp is 0.857868273815884
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

    ## Next temp is 0.153642666180727
    ## Next temp is 0.50829918722836
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-5.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-6.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-7.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-8.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-9.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-10.png)<!-- -->

    ## Next temp is 0.0194308207375097
    ## Next temp is 0.194321477484249
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-11.png)<!-- -->

    ## Next temp is 0.208092317357542
    ## Next temp is 0.416285603745412
    ## Next temp is 0.614351866171276
    ## Next temp is 0.829981167352962
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-12.png)<!-- -->

    ## Next temp is 0.0249817952942393
    ## Next temp is 0.0957354148649936
    ## Next temp is 0.202387931379878
    ## Next temp is 0.338687866536726
    ## Next temp is 0.483577987401768
    ## Next temp is 0.647518274156959
    ## Next temp is 0.82954635118449
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-13.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-14.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-15.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-16.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-17.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-18.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-19.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-20.png)<!-- -->

    ## Next temp is 0.00995985875274052
    ## Next temp is 0.100067859677949
    ## Next temp is 0.969700771806105
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-21.png)<!-- -->

    ## Next temp is 0.105925604212918
    ## Next temp is 0.214492818483211
    ## Next temp is 0.317379252084323
    ## Next temp is 0.425111636999605
    ## Next temp is 0.526364672571806
    ## Next temp is 0.630329019982785
    ## Next temp is 0.732567515720513
    ## Next temp is 0.840037288656411
    ## Next temp is 0.967453157536834
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-22.png)<!-- -->

    ## Next temp is 0.0146104359485384
    ## Next temp is 0.050975826479809
    ## Next temp is 0.104908457981391
    ## Next temp is 0.172878442728953
    ## Next temp is 0.249838282341316
    ## Next temp is 0.332588164352904
    ## Next temp is 0.419182574866912
    ## Next temp is 0.509040831837881
    ## Next temp is 0.598157551681333
    ## Next temp is 0.686961288855672
    ## Next temp is 0.774072483104258
    ## Next temp is 0.865289306425813
    ## Next temp is 0.961098506988116
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-23.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-24.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-25.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-26.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-27.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-28.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-29.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-30.png)<!-- -->

Overall, EKI covers both parameters quite well in this scenario.
However, it is worth noting that the prior for $\sigma^2$ already covers
the true parameter well and that the $\sigma^2$ particles have not
shifted significantly from their prior distribution. By contrast MCMC
concentrates well for both parameters and we can see a noticeable shift
from the prior distribution.

## Changing $\alpha$

In this experiment, I keep $\sigma^2$ fixed at 4, make $x$ a
50-dimensional vector of 1s and change $\alpha$ from 0 to 10 in
increments of 2. The purpose of this is to make sure that the EKI
algorithm is accurate across of range of $\alpha$ and if it changes how
the noise parameter is estimated.

``` r
alpha_seq <- seq(0, 10, length.out = 6)
num_dimensions <- 50 
num_particles <- 400
for (alpha in alpha_seq) {    
  true_parameters <- list(alpha = alpha, sigma = 2, x = rep(1, num_dimensions))   
  generate_results(num_particles, true_parameters)
}
```

    ## Next temp is 0.0204202489372129
    ## Next temp is 0.175823905042393
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

    ## Next temp is 0.18542020625145
    ## Next temp is 0.369377541149523
    ## Next temp is 0.555419603706681
    ## Next temp is 0.749102935069128
    ## Next temp is 0.929831455756061
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

    ## Next temp is 0.0345684054265834
    ## Next temp is 0.120788196135665
    ## Next temp is 0.251563413973423
    ## Next temp is 0.409843809840606
    ## Next temp is 0.585249154090942
    ## Next temp is 0.771237939893208
    ## Next temp is 0.959987650782649
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-5.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-6.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-7.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-8.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-9.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-10.png)<!-- -->

    ## Next temp is 0.0155514064540065
    ## Next temp is 0.146745343849284
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-11.png)<!-- -->

    ## Next temp is 0.163043389027133
    ## Next temp is 0.327040595745879
    ## Next temp is 0.485141281714457
    ## Next temp is 0.644932715897876
    ## Next temp is 0.79299094087241
    ## Next temp is 0.944536773832195
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-12.png)<!-- -->

    ## Next temp is 0.0245548680346588
    ## Next temp is 0.0926489433012542
    ## Next temp is 0.197552259811653
    ## Next temp is 0.32791246221421
    ## Next temp is 0.488453747160589
    ## Next temp is 0.67022001011305
    ## Next temp is 0.85967601794879
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-13.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-14.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-15.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-16.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-17.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-18.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-19.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-20.png)<!-- -->

    ## Next temp is 0.0109859742093441
    ## Next temp is 0.102609623406089
    ## Next temp is 0.851521290978515
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-21.png)<!-- -->

    ## Next temp is 0.248364338718228
    ## Next temp is 0.499062966322259
    ## Next temp is 0.786122323449135
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-22.png)<!-- -->

    ## Next temp is 0.0181101215074195
    ## Next temp is 0.0818007987459065
    ## Next temp is 0.200258825582731
    ## Next temp is 0.35846348329935
    ## Next temp is 0.549080190834945
    ## Next temp is 0.751800090215404
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-23.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-24.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-25.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-26.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-27.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-28.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-29.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-30.png)<!-- -->

    ## Next temp is 0.00723757735188418
    ## Next temp is 0.051257552676826
    ## Next temp is 0.414238333434658
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-31.png)<!-- -->

    ## Next temp is 0.269282005211813
    ## Next temp is 0.526484154669525
    ## Next temp is 0.789685538508991
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-32.png)<!-- -->

    ## Next temp is 0.0118228119979911
    ## Next temp is 0.0631036370702188
    ## Next temp is 0.158170791806545
    ## Next temp is 0.326542055335237
    ## Next temp is 0.54719505994247
    ## Next temp is 0.813894765167895
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-33.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-34.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-35.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-36.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-37.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-38.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-39.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-40.png)<!-- -->

    ## Next temp is 0.00374257292983237
    ## Next temp is 0.0198563961841701
    ## Next temp is 0.126093348333334
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-41.png)<!-- -->

    ## Next temp is 0.194923649584554
    ## Next temp is 0.400840500437064
    ## Next temp is 0.621616431039261
    ## Next temp is 0.829073894236793
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-42.png)<!-- -->

    ## Next temp is 0.00884258766004701
    ## Next temp is 0.0402609565489777
    ## Next temp is 0.111576311936736
    ## Next temp is 0.22726613667642
    ## Next temp is 0.382910839867645
    ## Next temp is 0.560218029343285
    ## Next temp is 0.769628326285667
    ## Next temp is 0.976595407414083
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-43.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-44.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-45.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-46.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-47.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-48.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-49.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-50.png)<!-- -->

    ## Next temp is 0.00276602785039855
    ## Next temp is 0.0136040865270129
    ## Next temp is 0.0942725742673419
    ## Next temp is 0.828153511127036
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-51.png)<!-- -->

    ## Next temp is 0.173004232238072
    ## Next temp is 0.353718069447931
    ## Next temp is 0.539266456552848
    ## Next temp is 0.748215888258394
    ## Next temp is 0.976629523608135
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-52.png)<!-- -->

    ## Next temp is 0.00529261256949334
    ## Next temp is 0.0204978099727375
    ## Next temp is 0.0693871945271121
    ## Next temp is 0.168520114849203
    ## Next temp is 0.297296113584747
    ## Next temp is 0.443868318133355
    ## Next temp is 0.622460143999576
    ## Next temp is 0.807237860517126
    ## Next temp is 0.989254948875243
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-53.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-54.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-55.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-56.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-57.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-58.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-59.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-60.png)<!-- -->

All EKI models and MCMC do well at estimating the mean parameter. There
is no change in how the noise parameter is being estimated across a
range of $\alpha$ values, which is not surprising since the covariance
matrices used to update the particles should not depend on $\alpha$.

## Changing $\sigma^2$

In this experiment I change the true value of $\sigma$ from $2$ to $10$
in increments of $2$. The purpose of this is to determine whether the
EKI algorithm can detect values of $\sigma$ which are further away from
the prior. Additionally, we would expect the $\alpha$ particles to
become more diffuse as the true value of $\sigma$ increases.

``` r
sigma_seq <- seq(2, 10, length.out = 2)
num_dimensions <- 50
num_particles <- 400
for (sigma in sigma_seq) {

  true_parameters <- list(alpha = 2, sigma = sigma, x = rep(1, num_dimensions))
  generate_results(num_particles, true_parameters)
}
```

    ## Next temp is 0.018074697550398
    ## Next temp is 0.153698786148632
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

    ## Next temp is 0.253538950686758
    ## Next temp is 0.509779004560826
    ## Next temp is 0.78276754801435
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

    ## Next temp is 0.0301539874997372
    ## Next temp is 0.113751452984853
    ## Next temp is 0.267878401018794
    ## Next temp is 0.480838235267461
    ## Next temp is 0.71916695264373
    ## Next temp is 0.985965688954379
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-5.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-6.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-7.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-8.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-9.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-10.png)<!-- -->

    ## Next temp is 0.428838691808616
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-11.png)<!-- -->

    ## Next temp is 0.0053303396023119
    ## Next temp is 0.0106132907950226
    ## Next temp is 0.0158987074006318
    ## Next temp is 0.0210482692681921
    ## Next temp is 0.0262238784982584
    ## Next temp is 0.0314290564696028
    ## Next temp is 0.0365798108457294
    ## Next temp is 0.0418453619984084
    ## Next temp is 0.0471680967062769
    ## Next temp is 0.0524243032912217
    ## Next temp is 0.0576383423371467
    ## Next temp is 0.0628568103833526
    ## Next temp is 0.0679531504835789
    ## Next temp is 0.0730177143163032
    ## Next temp is 0.0780592826548269
    ## Next temp is 0.0831529156713002
    ## Next temp is 0.088187074682006
    ## Next temp is 0.0933158305274567
    ## Next temp is 0.0984502095521412
    ## Next temp is 0.103609025466036
    ## Next temp is 0.108801869637954
    ## Next temp is 0.113992352543617
    ## Next temp is 0.119173700865116
    ## Next temp is 0.12434297071709
    ## Next temp is 0.129549852521595
    ## Next temp is 0.134815154777953
    ## Next temp is 0.1400526946199
    ## Next temp is 0.145257436326973
    ## Next temp is 0.150518516077045
    ## Next temp is 0.155709773965398
    ## Next temp is 0.161000528789561
    ## Next temp is 0.166335175935211
    ## Next temp is 0.171748601658212
    ## Next temp is 0.177157078577245
    ## Next temp is 0.182673802905295
    ## Next temp is 0.188150989284456
    ## Next temp is 0.193611443803134
    ## Next temp is 0.199100558531086
    ## Next temp is 0.204595727336292
    ## Next temp is 0.210039342187414
    ## Next temp is 0.215472730949893
    ## Next temp is 0.22089904244963
    ## Next temp is 0.226306663443495
    ## Next temp is 0.231651039132701
    ## Next temp is 0.236947368022866
    ## Next temp is 0.242185999013421
    ## Next temp is 0.247494634435809
    ## Next temp is 0.252837855099676
    ## Next temp is 0.258187534435996
    ## Next temp is 0.263529822199567
    ## Next temp is 0.268824275117075
    ## Next temp is 0.274123513588768
    ## Next temp is 0.279424529438331
    ## Next temp is 0.284779472178755
    ## Next temp is 0.290106999057895
    ## Next temp is 0.295468569658208
    ## Next temp is 0.300813615063802
    ## Next temp is 0.306156979659383
    ## Next temp is 0.311505054469156
    ## Next temp is 0.316856701003045
    ## Next temp is 0.322165206226455
    ## Next temp is 0.327488767050997
    ## Next temp is 0.332838761116779
    ## Next temp is 0.338238034703352
    ## Next temp is 0.343563940079075
    ## Next temp is 0.34886993890733
    ## Next temp is 0.354248027988625
    ## Next temp is 0.359636903841771
    ## Next temp is 0.364962879815816
    ## Next temp is 0.370194196691935
    ## Next temp is 0.375429383296382
    ## Next temp is 0.380664453963553
    ## Next temp is 0.386013801137576
    ## Next temp is 0.391481519229583
    ## Next temp is 0.39694448165989
    ## Next temp is 0.402395425490229
    ## Next temp is 0.407818979076868
    ## Next temp is 0.413245073197571
    ## Next temp is 0.418732439413895
    ## Next temp is 0.42423616671665
    ## Next temp is 0.429782638173024
    ## Next temp is 0.435301668914177
    ## Next temp is 0.440819806541006
    ## Next temp is 0.446417436600274
    ## Next temp is 0.452010676755025
    ## Next temp is 0.457640341158795
    ## Next temp is 0.463321060377687
    ## Next temp is 0.468968895495823
    ## Next temp is 0.474557986736685
    ## Next temp is 0.480189768087741
    ## Next temp is 0.485785479512251
    ## Next temp is 0.491355700759797
    ## Next temp is 0.496938630840027
    ## Next temp is 0.502556193835519
    ## Next temp is 0.508170082956678
    ## Next temp is 0.513765897960156
    ## Next temp is 0.519294383890469
    ## Next temp is 0.524816930900038
    ## Next temp is 0.530235043535726
    ## Next temp is 0.53564621926171
    ## Next temp is 0.54110626337548
    ## Next temp is 0.546528114152394
    ## Next temp is 0.551866216229484
    ## Next temp is 0.5571741166764
    ## Next temp is 0.562503391793885
    ## Next temp is 0.567854538253026
    ## Next temp is 0.5732547974044
    ## Next temp is 0.578704523020123
    ## Next temp is 0.584161900948626
    ## Next temp is 0.58956950174252
    ## Next temp is 0.595017925052978
    ## Next temp is 0.600447321963513
    ## Next temp is 0.605869431283725
    ## Next temp is 0.6113006057599
    ## Next temp is 0.61677752180075
    ## Next temp is 0.622310892076507
    ## Next temp is 0.627967610923185
    ## Next temp is 0.633701886642208
    ## Next temp is 0.639360236295889
    ## Next temp is 0.644981785115446
    ## Next temp is 0.650619952773377
    ## Next temp is 0.656200090530194
    ## Next temp is 0.661796920577093
    ## Next temp is 0.667438591873142
    ## Next temp is 0.673114714654511
    ## Next temp is 0.678902392953873
    ## Next temp is 0.684634505088309
    ## Next temp is 0.690357557863697
    ## Next temp is 0.696121416240936
    ## Next temp is 0.701875336905562
    ## Next temp is 0.707621100535038
    ## Next temp is 0.713301196097776
    ## Next temp is 0.718872062429613
    ## Next temp is 0.724458691244341
    ## Next temp is 0.730008831685376
    ## Next temp is 0.735514859352769
    ## Next temp is 0.741011561717887
    ## Next temp is 0.746554796173287
    ## Next temp is 0.752098802886693
    ## Next temp is 0.757597112680645
    ## Next temp is 0.763161493069577
    ## Next temp is 0.768727570040843
    ## Next temp is 0.774230575546879
    ## Next temp is 0.779827092970895
    ## Next temp is 0.785467856779534
    ## Next temp is 0.791135007017349
    ## Next temp is 0.796732268955458
    ## Next temp is 0.802355374653086
    ## Next temp is 0.807964456736252
    ## Next temp is 0.813635496050408
    ## Next temp is 0.819438216822466
    ## Next temp is 0.825231290561792
    ## Next temp is 0.831034135969401
    ## Next temp is 0.836855770122367
    ## Next temp is 0.842717542314865
    ## Next temp is 0.848531738761833
    ## Next temp is 0.85424899621929
    ## Next temp is 0.859948268315559
    ## Next temp is 0.865576095566439
    ## Next temp is 0.871281416686324
    ## Next temp is 0.877069973498734
    ## Next temp is 0.882817244804842
    ## Next temp is 0.888639129867279
    ## Next temp is 0.894457572814708
    ## Next temp is 0.900215777002622
    ## Next temp is 0.906031811496029
    ## Next temp is 0.91179346070259
    ## Next temp is 0.917564869266399
    ## Next temp is 0.923342308114261
    ## Next temp is 0.929119136655984
    ## Next temp is 0.934926285595032
    ## Next temp is 0.940672176659443
    ## Next temp is 0.946408348971278
    ## Next temp is 0.952134832094445
    ## Next temp is 0.95779462521981
    ## Next temp is 0.963376733973834
    ## Next temp is 0.968931835050972
    ## Next temp is 0.974495688256572
    ## Next temp is 0.980073535343877
    ## Next temp is 0.985655458548857
    ## Next temp is 0.991250419747613
    ## Next temp is 0.996876871721706
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-12.png)<!-- -->

    ## Next temp is 0.00425008955849161
    ## Next temp is 0.00879994992711752
    ## Next temp is 0.0135912604181889
    ## Next temp is 0.0185115235663068
    ## Next temp is 0.0235302951044751
    ## Next temp is 0.0286179192031597
    ## Next temp is 0.0337769271295703
    ## Next temp is 0.038948354926249
    ## Next temp is 0.0440791338347297
    ## Next temp is 0.0493175298181347
    ## Next temp is 0.054639539101463
    ## Next temp is 0.0599712548350425
    ## Next temp is 0.0653310129530087
    ## Next temp is 0.0707162799289747
    ## Next temp is 0.0763178977617192
    ## Next temp is 0.0819560804158868
    ## Next temp is 0.0876316766485527
    ## Next temp is 0.0931638961811253
    ## Next temp is 0.0986661837395562
    ## Next temp is 0.104259282076634
    ## Next temp is 0.109919101134658
    ## Next temp is 0.115602020304633
    ## Next temp is 0.121300392298417
    ## Next temp is 0.127069634106533
    ## Next temp is 0.132767299406192
    ## Next temp is 0.13849112952757
    ## Next temp is 0.144302244536967
    ## Next temp is 0.150089873541179
    ## Next temp is 0.155859757165017
    ## Next temp is 0.161647922976246
    ## Next temp is 0.167386916915232
    ## Next temp is 0.173310681253567
    ## Next temp is 0.179329536064696
    ## Next temp is 0.1853647370877
    ## Next temp is 0.191223431788732
    ## Next temp is 0.196952916104527
    ## Next temp is 0.202684584457345
    ## Next temp is 0.208425461350571
    ## Next temp is 0.214275604734055
    ## Next temp is 0.220270140292877
    ## Next temp is 0.226244664088166
    ## Next temp is 0.232198182164257
    ## Next temp is 0.238126753836431
    ## Next temp is 0.243985050606232
    ## Next temp is 0.249822624287248
    ## Next temp is 0.255672064556624
    ## Next temp is 0.261626445357864
    ## Next temp is 0.267528822187097
    ## Next temp is 0.273470175165111
    ## Next temp is 0.279305480946284
    ## Next temp is 0.285120144448242
    ## Next temp is 0.290952951105256
    ## Next temp is 0.296782378268376
    ## Next temp is 0.302486760446965
    ## Next temp is 0.308149893126897
    ## Next temp is 0.31381229457139
    ## Next temp is 0.319492913259348
    ## Next temp is 0.325221071704927
    ## Next temp is 0.330823549085456
    ## Next temp is 0.336338358327769
    ## Next temp is 0.34184175610584
    ## Next temp is 0.347383894340236
    ## Next temp is 0.352935308281654
    ## Next temp is 0.358446195618837
    ## Next temp is 0.364092565700834
    ## Next temp is 0.369805934888457
    ## Next temp is 0.375592480561119
    ## Next temp is 0.381427537269835
    ## Next temp is 0.387221663870087
    ## Next temp is 0.393041201074597
    ## Next temp is 0.398831657931916
    ## Next temp is 0.404727260996446
    ## Next temp is 0.410673056709386
    ## Next temp is 0.416609612884796
    ## Next temp is 0.422544053939414
    ## Next temp is 0.42854825236062
    ## Next temp is 0.434623985867702
    ## Next temp is 0.440738323029581
    ## Next temp is 0.446779916003798
    ## Next temp is 0.452809535483096
    ## Next temp is 0.45884762846789
    ## Next temp is 0.464874162267392
    ## Next temp is 0.471005152150637
    ## Next temp is 0.477165528200514
    ## Next temp is 0.483298196657015
    ## Next temp is 0.489479661910179
    ## Next temp is 0.49569824943855
    ## Next temp is 0.501933622753622
    ## Next temp is 0.508126372313014
    ## Next temp is 0.514371351870601
    ## Next temp is 0.520548057811629
    ## Next temp is 0.526705674699907
    ## Next temp is 0.532910981124429
    ## Next temp is 0.539046442271134
    ## Next temp is 0.54525470048257
    ## Next temp is 0.5514310343
    ## Next temp is 0.557568594403946
    ## Next temp is 0.56375348840048
    ## Next temp is 0.569993578883179
    ## Next temp is 0.576189431839605
    ## Next temp is 0.582479237947686
    ## Next temp is 0.588716059998128
    ## Next temp is 0.594896794328265
    ## Next temp is 0.601078178856362
    ## Next temp is 0.607311470780431
    ## Next temp is 0.613464668570116
    ## Next temp is 0.619505838189572
    ## Next temp is 0.625626301963046
    ## Next temp is 0.631788502709293
    ## Next temp is 0.638037871323399
    ## Next temp is 0.644345029785364
    ## Next temp is 0.650703878587705
    ## Next temp is 0.657017669640215
    ## Next temp is 0.66333289665038
    ## Next temp is 0.669625898443885
    ## Next temp is 0.675926974795427
    ## Next temp is 0.682100368036118
    ## Next temp is 0.688231648436049
    ## Next temp is 0.694361822067494
    ## Next temp is 0.700627312796146
    ## Next temp is 0.706905500042848
    ## Next temp is 0.713177733887621
    ## Next temp is 0.719492564359834
    ## Next temp is 0.725761411397448
    ## Next temp is 0.732127229709252
    ## Next temp is 0.738496569330338
    ## Next temp is 0.744993492541933
    ## Next temp is 0.751404013967568
    ## Next temp is 0.757869534635992
    ## Next temp is 0.764387227451732
    ## Next temp is 0.770990029272298
    ## Next temp is 0.777696975920153
    ## Next temp is 0.784422140313744
    ## Next temp is 0.791153914155219
    ## Next temp is 0.797824313294271
    ## Next temp is 0.804424966804915
    ## Next temp is 0.811169943641201
    ## Next temp is 0.81791461218963
    ## Next temp is 0.824531902825963
    ## Next temp is 0.831130744941228
    ## Next temp is 0.837694174677598
    ## Next temp is 0.844215331796346
    ## Next temp is 0.850851645874904
    ## Next temp is 0.857580219972261
    ## Next temp is 0.864366592015495
    ## Next temp is 0.871175593008172
    ## Next temp is 0.878055426642179
    ## Next temp is 0.884831133723877
    ## Next temp is 0.891661230342365
    ## Next temp is 0.898553859156836
    ## Next temp is 0.905561305905564
    ## Next temp is 0.912579502026605
    ## Next temp is 0.919582393500809
    ## Next temp is 0.926593626531315
    ## Next temp is 0.933588999903902
    ## Next temp is 0.940717655743129
    ## Next temp is 0.947821540349686
    ## Next temp is 0.954918735276001
    ## Next temp is 0.961954555805828
    ## Next temp is 0.968974707424676
    ## Next temp is 0.976133186542264
    ## Next temp is 0.983223052495174
    ## Next temp is 0.990382520873174
    ## Next temp is 0.997577993160633
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-13.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-14.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-15.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-16.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-17.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-18.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-19.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-20.png)<!-- -->

In the known variance case, the results for the mean parameter are what
we expect. That is, the particles are centered around the true value and
they become more diffuse as $\sigma^2$ increases. However, the ability
of the other EKI models to estimate the true noise parameter is poor.
Regardless of the true value of $\sigma^2$, the particles remain
anchored to their prior distribution. In the cases where the prior is
lower than the true value of $\sigma^2$, the $\alpha$ values become too
concentrated. This is problematic if we don’t have strong information
about the noise parameter to specify a good prior distribution.

## Changing Number of Particles

Finally I change the number of particles used for the EKI algorithm. The
main purpose of this is to check that the EKI model with known variance
concentrates around the true posterior distribution.

``` r
num_dimensions <- 50
particle_seq <- c(100, 500, 1000) 
true_parameters <- list(alpha = 2, sigma = 2, x = rep(1, num_dimensions))

for (num_particles in particle_seq) {
  generate_results(num_particles, true_parameters)
}
```

    ## Next temp is 0.0168005819181129
    ## Next temp is 0.313966746401731
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

    ## Next temp is 0.160407189656696
    ## Next temp is 0.413308349267595
    ## Next temp is 0.710883241727487
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-13-2.png)<!-- -->

    ## Next temp is 0.0312629153870564
    ## Next temp is 0.156404582866844
    ## Next temp is 0.373391611019594
    ## Next temp is 0.620192795624161
    ## Next temp is 0.85190394143846
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-13-3.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-4.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-5.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-6.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-7.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-8.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-9.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-10.png)<!-- -->

    ## Next temp is 0.0175424215807934
    ## Next temp is 0.162074418513314
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-13-11.png)<!-- -->

    ## Next temp is 0.276061527822876
    ## Next temp is 0.526233022392729
    ## Next temp is 0.779019508100977
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-13-12.png)<!-- -->

    ## Next temp is 0.0336527279792067
    ## Next temp is 0.11686055920205
    ## Next temp is 0.24660142514884
    ## Next temp is 0.409143142062176
    ## Next temp is 0.588772739011596
    ## Next temp is 0.776875105371088
    ## Next temp is 0.971029661775292
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-13-13.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-14.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-15.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-16.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-17.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-18.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-19.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-20.png)<!-- -->

    ## Next temp is 0.0167203212074068
    ## Next temp is 0.153368295575384
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-13-21.png)<!-- -->

    ## Next temp is 0.16919797202442
    ## Next temp is 0.336901792702724
    ## Next temp is 0.507565652463643
    ## Next temp is 0.672931616445919
    ## Next temp is 0.836655171666221
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-13-22.png)<!-- -->

    ## Next temp is 0.0290224881886644
    ## Next temp is 0.0984484842828443
    ## Next temp is 0.195205607463718
    ## Next temp is 0.310100516809191
    ## Next temp is 0.437780492936273
    ## Next temp is 0.573777026795985
    ## Next temp is 0.721434942701475
    ## Next temp is 0.875822109279098
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-13-23.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-24.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-25.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-26.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-27.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-28.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-29.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-13-30.png)<!-- -->

Indeed we find that this is the case.

## Lognormal Model

Considering that the normal model does poorly when looking at the noise
parameter, we try an alternative model specification. Specifically, we
look at the lognormal distribution. Unlike the normal distribution, the
noise parameter effects the mean of the observations, so we should
expect a signal with the data. I only look at changing the noise
parameter in this case.

``` r
prior_params <- list(alpha.mean = 0, alpha.sd = 1, 
                     sigma2.mean = -0.5, sigma2.sd = 1)
```

``` r
sigma_seq <- seq(0.5, 2.5, length.out = 5)
num_dimensions <- 50
num_particles <- 400
for (sigma in sigma_seq) {

  true_parameters <- list(alpha = 2, sigma = sigma, x = rep(1, num_dimensions))
  generate_lognormal_results(num_particles, true_parameters, prior_params)
}
```

    ## Next temp is 0.19916654331749
    ## Next temp is 0.461807502698022
    ## Next temp is 0.822569780373638
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

    ## Next temp is 0.0112116861912029
    ## Next temp is 0.0238528798866366
    ## Next temp is 0.0374043024408499
    ## Next temp is 0.0521147327733864
    ## Next temp is 0.0681289022600505
    ## Next temp is 0.0859733203866159
    ## Next temp is 0.105380250835933
    ## Next temp is 0.12714710685059
    ## Next temp is 0.15140493940014
    ## Next temp is 0.178916293533801
    ## Next temp is 0.209928273995097
    ## Next temp is 0.244346210710543
    ## Next temp is 0.284175895544516
    ## Next temp is 0.329538552371244
    ## Next temp is 0.378737914565093
    ## Next temp is 0.432488876157138
    ## Next temp is 0.49300638672172
    ## Next temp is 0.556599081394025
    ## Next temp is 0.623146391703158
    ## Next temp is 0.696435757976091
    ## Next temp is 0.774080050698594
    ## Next temp is 0.85552627367261
    ## Next temp is 0.938830622643806
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-15-2.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-3.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-4.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-5.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-6.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-7.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-8.png)<!-- -->

    ## Next temp is 0.130715905332
    ## Next temp is 0.259443755009631
    ## Next temp is 0.414244759183395
    ## Next temp is 0.751883702917038
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-15-9.png)<!-- -->

    ## Next temp is 0.0091118836097642
    ## Next temp is 0.0204045386427113
    ## Next temp is 0.0331417859892501
    ## Next temp is 0.0499365906301914
    ## Next temp is 0.0694451077385043
    ## Next temp is 0.0924738120280769
    ## Next temp is 0.121498591677243
    ## Next temp is 0.157684929382882
    ## Next temp is 0.200196428298955
    ## Next temp is 0.251876018831406
    ## Next temp is 0.309438196188045
    ## Next temp is 0.375513986806074
    ## Next temp is 0.448452173266884
    ## Next temp is 0.550560128334593
    ## Next temp is 0.654984610637353
    ## Next temp is 0.773420199535261
    ## Next temp is 0.916393415933545
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-15-10.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-11.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-12.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-13.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-14.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-15.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-16.png)<!-- -->

    ## Next temp is 0.038719837528134
    ## Next temp is 0.0773217452684478
    ## Next temp is 0.113426990653692
    ## Next temp is 0.152830695757509
    ## Next temp is 0.20525267104835
    ## Next temp is 0.27715397877297
    ## Next temp is 0.369126159083577
    ## Next temp is 0.480076334587506
    ## Next temp is 0.640052475824755
    ## Next temp is 0.833808860022667
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-15-17.png)<!-- -->

    ## Next temp is 0.00726036116396764
    ## Next temp is 0.0167275847392116
    ## Next temp is 0.0283764414670398
    ## Next temp is 0.0425769239136372
    ## Next temp is 0.0591612669348966
    ## Next temp is 0.0792520291726114
    ## Next temp is 0.103304907381016
    ## Next temp is 0.129507183345708
    ## Next temp is 0.157396393893557
    ## Next temp is 0.191835952665203
    ## Next temp is 0.231991143276082
    ## Next temp is 0.284164573274131
    ## Next temp is 0.34348978100547
    ## Next temp is 0.41275525252468
    ## Next temp is 0.50348022340493
    ## Next temp is 0.610031740877824
    ## Next temp is 0.716728996519583
    ## Next temp is 0.836868661739435
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-15-18.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-19.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-20.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-21.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-22.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-23.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-24.png)<!-- -->

    ## Next temp is 0.0167808079763405
    ## Next temp is 0.0359515395626651
    ## Next temp is 0.0638742409874259
    ## Next temp is 0.130976062329949
    ## Next temp is 0.241292262417142
    ## Next temp is 0.330048831029776
    ## Next temp is 0.457213360889159
    ## Next temp is 0.691917132876133
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-15-25.png)<!-- -->

    ## Next temp is 0.00565303993428789
    ## Next temp is 0.0129282358236919
    ## Next temp is 0.0219973062071527
    ## Next temp is 0.0321447388848473
    ## Next temp is 0.0506138006281967
    ## Next temp is 0.0729646678452675
    ## Next temp is 0.10002989625664
    ## Next temp is 0.128571389327323
    ## Next temp is 0.158538223696624
    ## Next temp is 0.195586214156795
    ## Next temp is 0.231403978138345
    ## Next temp is 0.267201205968365
    ## Next temp is 0.305098363896664
    ## Next temp is 0.343208738838065
    ## Next temp is 0.38029693157603
    ## Next temp is 0.42036102547836
    ## Next temp is 0.462698569804139
    ## Next temp is 0.506824434773416
    ## Next temp is 0.549665925642895
    ## Next temp is 0.590209930042035
    ## Next temp is 0.632399967765749
    ## Next temp is 0.671944621855146
    ## Next temp is 0.71103935582321
    ## Next temp is 0.753934485624172
    ## Next temp is 0.798334084584441
    ## Next temp is 0.844371949470087
    ## Next temp is 0.881455596856146
    ## Next temp is 0.915700923989179
    ## Next temp is 0.94862597327758
    ## Next temp is 0.979393578072594
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-15-26.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-27.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-28.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-29.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-30.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-31.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-32.png)<!-- -->

    ## Next temp is 0.0120691326918458
    ## Next temp is 0.0255172145003302
    ## Next temp is 0.0442573912444094
    ## Next temp is 0.0632480411901581
    ## Next temp is 0.0838068635384067
    ## Next temp is 0.111326735561598
    ## Next temp is 0.154737431616794
    ## Next temp is 0.19773032874503
    ## Next temp is 0.252593669665505
    ## Next temp is 0.551896241458079
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-15-33.png)<!-- -->

    ## Next temp is 0.00573951706763352
    ## Next temp is 0.0141608181694062
    ## Next temp is 0.0284955821903802
    ## Next temp is 0.0514748953322122
    ## Next temp is 0.074540032421508
    ## Next temp is 0.0999470155701103
    ## Next temp is 0.124163473800316
    ## Next temp is 0.153231831382424
    ## Next temp is 0.19757095478733
    ## Next temp is 0.265151848000563
    ## Next temp is 0.339041358983855
    ## Next temp is 0.421482175849051
    ## Next temp is 0.500534069350087
    ## Next temp is 0.635184752555274
    ## Next temp is 0.845212117429116
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-15-34.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-35.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-36.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-37.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-38.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-39.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-15-40.png)<!-- -->

## Summary Statistics

In this section, we further summarize the data into two summary
statistics $(\bar{y}, s_y)$ which are jointly sufficient for the
parameters. This is because we expect a correlation between $\sigma^2$
and $s_y$ which should allow us to estimate the parameters. Because they
are sufficient, the posterior is still exact in the linear Gaussian
case.

We can see that GEKI does well at estimating both parameters in this
model when we use summary statistics. This is interesting considering it
doesn’t do well when using the full data.

## Conclusion

To summarize, GEKI does well at estimating the mean parameter, but does
poorly at estimating the noise parameter, remaining anchored to its
prior distribution. By contrast, MCMC does well at estimating both
parameters. One explanation for this phenomenon is that the GEKI
algorithm relies on a cross covariance between the parameters and the
data, but $\sigma^2$ and $y$ are theoretically uncorrelated. See
`results_covariances.md` for further detail.
