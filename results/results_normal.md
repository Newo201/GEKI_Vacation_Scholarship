Generalised Ensemble Kalman Inversion - Multivariate Normal Model
================
Owen Jackson
2025-01-15

# Results For Multivariate Normal Model

## Import Functions

### External Packages

``` r
pacman::p_load(pacman, purrr, mvtnorm, mcmc, MASS)
```

### Algorithms

``` r
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/eki.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/mcmc_normal.R')
```

### Models

``` r
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_normal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_normal_known_var.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_normal_known_mean.R')
```

### Sampling and PDFs

``` r
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/pdfs/pdfs_normal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples/samples_normal.R')
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
  plot_mcmc_trace_plots(chain.1, chain.2, 500)
  plot_mcmc_histogram(chain.1, chain.2, true_parameters, 500)
  
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

![](results_normal_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-5.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-6.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-7.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-8.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-9.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-10.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-11.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-12.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-13.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-14.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-15.png)<!-- -->

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

![](results_normal_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-5.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-6.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-7.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-8.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-9.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-10.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-11.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-12.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-13.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-14.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-15.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-16.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-17.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-18.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-19.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-20.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-21.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-22.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-23.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-24.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-25.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-26.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-27.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-28.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-29.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-30.png)<!-- -->

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
sigma_seq <- seq(2, 10, length.out = 5)
num_dimensions <- 50
num_particles <- 400
for (sigma in sigma_seq) {

  true_parameters <- list(alpha = 2, sigma = sigma, x = rep(1, num_dimensions))
  generate_results(num_particles, true_parameters)
}
```

![](results_normal_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-5.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-6.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-7.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-8.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-9.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-10.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-11.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-12.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-13.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-14.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-15.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-16.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-17.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-18.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-19.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-20.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-21.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-22.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-23.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-24.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-25.png)<!-- -->

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

![](results_normal_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-5.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-6.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-7.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-8.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-9.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-10.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-11.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-12.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-13.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-14.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-15.png)<!-- -->

Indeed we find that this is the case.

## Conclusion

To summarize, GEKI does well at estimating the mean parameter, but does
poorly at estimating the noise parameter, remaining anchored to its
prior distribution. By contrast, MCMC does well at estimating both
parameters. One explanation for this phenomenon is that the GEKI
algorithm relies on a cross covariance between the parameters and the
data, but $\sigma^2$ and $y$ are theoretically uncorrelated. See
`results_covariances.md` for further detail.

While beyond the scope of this project, one would expect this phenomenon
to be observed in similar parameters that affect variance and kurtosis
only (e.g. the degrees of freedom in a t-distribution). Improvements to
the algorithm which may be able to address this issue are also beyond
the scope of this project.
