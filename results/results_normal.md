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
  
  true_sample <- likelihood_normal(true_parameters)
  true_data <- c(true_sample, sd(true_sample))
  
  # Model 1
  eki_result_known_var <- eki_normal_known_var(num_particles, true_data, true_parameters, 
                                               prior_params, adaptive = adaptive)
  plot_eki_normal_known_var(eki_result_known_var, true_sample, 
                            true_parameters, prior_params)
  
  # Model 2
  eki_result_known_mean <- eki_normal_known_mean(num_particles, true_data, true_parameters, 
                                                 prior_params, adaptive = adaptive)
  plot_eki_normal_known_mean(eki_result_known_mean, true_parameters, prior_params)
  
  # Model 3
  eki_result <- eki_normal(num_particles, true_data, true_parameters, 
                           prior_params, adaptive = adaptive)
  plot_eki_normal(eki_result, true_parameters, prior_params)
  
  # # MCMC Comparison for Model 3
  # run.1 <- normal_mcmc(true_data, true_parameters, prior_params)
  # chain.1 <- run.1$batch
  # run.2 <- normal_mcmc(true_data, true_parameters, prior_params)
  # chain.2 <- run.2$batch
  # plot_mcmc_trace_plots(chain.1, chain.2, burnin = 500)
  # plot_mcmc_histogram(chain.1, chain.2, true_parameters, prior_params, 
  #                     burnin = 500)
  
}
```

## The Model

We assume that we have access to a single observation from a
multivariate normal distribution $y \sim N(\alpha x, \sigma^2I)$ where
$\alpha$ is a scalar and $x$ is a known vector. We further augment this
data with the sample standard deviation $s_y$.

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
  # true_sample <- likelihood_normal(true_parameters)
  # true_data <- c(mean(true_sample), sd(true_sample))
  # eki_result <- eki_normal(num_particles, true_data, true_parameters, 
  #                        prior_params, adaptive = adaptive)
  # plot_eki_normal(eki_result, true_parameters, prior_params)
}
```

    ## Next temp is 0.0663722386824061
    ## Next temp is 0.286601703114163
    ## Next temp is 0.530821185893534
    ## Next temp is 0.808986256483745
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

    ## Next temp is 0.414124616862147
    ## Next temp is 0.767358666552512
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->

    ## Next temp is 0.156400109408777
    ## Next temp is 0.423104485047443
    ## Next temp is 0.717282545114891
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->

    ## Next temp is 0.0159232349938719
    ## Next temp is 0.080579855714708
    ## Next temp is 0.175422814916486
    ## Next temp is 0.27761693205788
    ## Next temp is 0.381176236713693
    ## Next temp is 0.491093549500933
    ## Next temp is 0.591715202033683
    ## Next temp is 0.693850990167666
    ## Next temp is 0.8048261926058
    ## Next temp is 0.914840472964458
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-9-5.png)<!-- -->

    ## Next temp is 0.0968266388962507
    ## Next temp is 0.190531148428993
    ## Next temp is 0.28794365908653
    ## Next temp is 0.385190824626928
    ## Next temp is 0.485503524831705
    ## Next temp is 0.592020187819491
    ## Next temp is 0.690314548311411
    ## Next temp is 0.782775674482756
    ## Next temp is 0.887720394481156
    ## Next temp is 0.988027010231292
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-9-6.png)<!-- -->

    ## Next temp is 0.0246128909500824
    ## Next temp is 0.0743082900738266
    ## Next temp is 0.145950213484078
    ## Next temp is 0.230318385361513
    ## Next temp is 0.319903155552539
    ## Next temp is 0.422254237447387
    ## Next temp is 0.508978557441594
    ## Next temp is 0.594324756262281
    ## Next temp is 0.696754762008449
    ## Next temp is 0.809825059992278
    ## Next temp is 0.920820003024172
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-9-7.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-8.png)<!-- -->

    ## Next temp is 0.00626976698748959
    ## Next temp is 0.0404509024306742
    ## Next temp is 0.0944836201658137
    ## Next temp is 0.165370584941028
    ## Next temp is 0.237780938811783
    ## Next temp is 0.307811618182495
    ## Next temp is 0.372429766305375
    ## Next temp is 0.438279006952206
    ## Next temp is 0.510630615474525
    ## Next temp is 0.583429990463426
    ## Next temp is 0.656489562661078
    ## Next temp is 0.730428010430272
    ## Next temp is 0.806623533400516
    ## Next temp is 0.87742445781236
    ## Next temp is 0.951392567008668
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-9-9.png)<!-- -->

    ## Next temp is 0.0443979133950869
    ## Next temp is 0.0954856694503441
    ## Next temp is 0.151086475262517
    ## Next temp is 0.210393710250987
    ## Next temp is 0.27578516588109
    ## Next temp is 0.342892724082407
    ## Next temp is 0.40859500487348
    ## Next temp is 0.466762654576915
    ## Next temp is 0.529791130434685
    ## Next temp is 0.593595649749437
    ## Next temp is 0.663745957616508
    ## Next temp is 0.742832774760197
    ## Next temp is 0.806234895296754
    ## Next temp is 0.878878323179921
    ## Next temp is 0.93854146692227
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-9-10.png)<!-- -->

    ## Next temp is 0.0118899534163036
    ## Next temp is 0.0345269658904381
    ## Next temp is 0.0675126791926049
    ## Next temp is 0.101906740796542
    ## Next temp is 0.152192848615069
    ## Next temp is 0.203649578126522
    ## Next temp is 0.26085318280886
    ## Next temp is 0.318717750141703
    ## Next temp is 0.376760614598095
    ## Next temp is 0.434599200663659
    ## Next temp is 0.496824882563265
    ## Next temp is 0.57351297257659
    ## Next temp is 0.637909414798081
    ## Next temp is 0.6970026010055
    ## Next temp is 0.765731689734511
    ## Next temp is 0.827681408051776
    ## Next temp is 0.894073260604136
    ## Next temp is 0.970023987625907
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-9-11.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-9-12.png)<!-- -->

Overall, EKI is able to estimate both parameters quite well in this
scenario. We can see a significant movement away from the prior
distribution for both parameters, which is promising. Additionally, in
the known variance scenario, the posterior is still exact, which is what
we expect with sufficient statistics.

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

    ## Next temp is 0.0144809873248642
    ## Next temp is 0.0772813787150981
    ## Next temp is 0.170798358199618
    ## Next temp is 0.272486697552656
    ## Next temp is 0.376356027848421
    ## Next temp is 0.477663908985177
    ## Next temp is 0.584235705638273
    ## Next temp is 0.694684922390558
    ## Next temp is 0.821399608594613
    ## Next temp is 0.928624415870984
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

    ## Next temp is 0.0999789045418
    ## Next temp is 0.189110803690641
    ## Next temp is 0.294068398076003
    ## Next temp is 0.3944171821208
    ## Next temp is 0.497819247526897
    ## Next temp is 0.591821230803007
    ## Next temp is 0.694711881711052
    ## Next temp is 0.802708025929122
    ## Next temp is 0.901078555825336
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

    ## Next temp is 0.0257251063297429
    ## Next temp is 0.0798640791726288
    ## Next temp is 0.150845191368852
    ## Next temp is 0.229536511723272
    ## Next temp is 0.315273027566377
    ## Next temp is 0.41381262738128
    ## Next temp is 0.520739273484486
    ## Next temp is 0.611120374224664
    ## Next temp is 0.711482472482838
    ## Next temp is 0.801405589265275
    ## Next temp is 0.902944039541763
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-3.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-4.png)<!-- -->

    ## Next temp is 0.0130303813972454
    ## Next temp is 0.0655411906043328
    ## Next temp is 0.137775863016858
    ## Next temp is 0.233383729126749
    ## Next temp is 0.330487608812669
    ## Next temp is 0.435845244843768
    ## Next temp is 0.546324427949597
    ## Next temp is 0.646502934474777
    ## Next temp is 0.762224596691969
    ## Next temp is 0.861104426288302
    ## Next temp is 0.965579650420128
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-5.png)<!-- -->

    ## Next temp is 0.100680767859887
    ## Next temp is 0.21489035186147
    ## Next temp is 0.31815742447507
    ## Next temp is 0.413920105427979
    ## Next temp is 0.513290230629973
    ## Next temp is 0.616683456629398
    ## Next temp is 0.734772584669165
    ## Next temp is 0.83861017005229
    ## Next temp is 0.936693695595008
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-6.png)<!-- -->

    ## Next temp is 0.0237271668756541
    ## Next temp is 0.0760119745471268
    ## Next temp is 0.144752845903668
    ## Next temp is 0.227199129784127
    ## Next temp is 0.326178485434939
    ## Next temp is 0.434177339155251
    ## Next temp is 0.534956677866143
    ## Next temp is 0.623441375506481
    ## Next temp is 0.724835433234697
    ## Next temp is 0.845816427994876
    ## Next temp is 0.952116216120036
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-7.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-8.png)<!-- -->

    ## Next temp is 0.008855295169334
    ## Next temp is 0.0546168948551043
    ## Next temp is 0.13270382698229
    ## Next temp is 0.229946188065067
    ## Next temp is 0.331320335671745
    ## Next temp is 0.438096534334498
    ## Next temp is 0.550755320961714
    ## Next temp is 0.654995601607082
    ## Next temp is 0.77295131673697
    ## Next temp is 0.872036509396466
    ## Next temp is 0.980999758704467
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-9.png)<!-- -->

    ## Next temp is 0.0833920297331841
    ## Next temp is 0.178324509765058
    ## Next temp is 0.270462922709667
    ## Next temp is 0.370538471171746
    ## Next temp is 0.479371410680246
    ## Next temp is 0.59275404982171
    ## Next temp is 0.692319964861436
    ## Next temp is 0.785999125004146
    ## Next temp is 0.903151837853462
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-10.png)<!-- -->

    ## Next temp is 0.0174865421769199
    ## Next temp is 0.0633331320705049
    ## Next temp is 0.118870006592336
    ## Next temp is 0.200257999011128
    ## Next temp is 0.281618573973706
    ## Next temp is 0.373420560104373
    ## Next temp is 0.474132735937428
    ## Next temp is 0.57607455240095
    ## Next temp is 0.685445907235655
    ## Next temp is 0.783959595885767
    ## Next temp is 0.889856312917872
    ## Next temp is 0.992864539389064
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-11.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-12.png)<!-- -->

    ## Next temp is 0.006211685191458
    ## Next temp is 0.0440625533881076
    ## Next temp is 0.128390302661966
    ## Next temp is 0.228126754087217
    ## Next temp is 0.328743487001627
    ## Next temp is 0.441718075088194
    ## Next temp is 0.544034847672454
    ## Next temp is 0.648193686735425
    ## Next temp is 0.751459687773993
    ## Next temp is 0.858683916146792
    ## Next temp is 0.964493732400772
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-13.png)<!-- -->

    ## Next temp is 0.104793565451647
    ## Next temp is 0.195976856377863
    ## Next temp is 0.29331494257051
    ## Next temp is 0.389634775486219
    ## Next temp is 0.492295263237608
    ## Next temp is 0.588737519692811
    ## Next temp is 0.690560399967952
    ## Next temp is 0.797541263185095
    ## Next temp is 0.894680826600142
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-14.png)<!-- -->

    ## Next temp is 0.0129726741366896
    ## Next temp is 0.0497613473590634
    ## Next temp is 0.114485756347468
    ## Next temp is 0.194388559824385
    ## Next temp is 0.283477515145029
    ## Next temp is 0.370103345958855
    ## Next temp is 0.4736235764864
    ## Next temp is 0.566486224552591
    ## Next temp is 0.674118697911552
    ## Next temp is 0.77972767006694
    ## Next temp is 0.884387486990513
    ## Next temp is 0.99526376834102
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-15.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-16.png)<!-- -->

    ## Next temp is 0.00407809283558522
    ## Next temp is 0.0230746763272435
    ## Next temp is 0.0772284581382804
    ## Next temp is 0.166800106000885
    ## Next temp is 0.26014110632299
    ## Next temp is 0.35728257854952
    ## Next temp is 0.464687235148139
    ## Next temp is 0.558149071944147
    ## Next temp is 0.655003466887107
    ## Next temp is 0.753844757268642
    ## Next temp is 0.842246841680275
    ## Next temp is 0.939248316806412
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-17.png)<!-- -->

    ## Next temp is 0.0867579812421335
    ## Next temp is 0.1913892526915
    ## Next temp is 0.286232350531644
    ## Next temp is 0.385791615145375
    ## Next temp is 0.470342830474749
    ## Next temp is 0.568459777860514
    ## Next temp is 0.675103533053063
    ## Next temp is 0.77891022226881
    ## Next temp is 0.884412644022732
    ## Next temp is 0.987408812613212
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-18.png)<!-- -->

    ## Next temp is 0.00985754312981875
    ## Next temp is 0.0326733118685458
    ## Next temp is 0.0820562897760152
    ## Next temp is 0.1556170719954
    ## Next temp is 0.241182570656025
    ## Next temp is 0.329322671653149
    ## Next temp is 0.427891267239759
    ## Next temp is 0.520771471700151
    ## Next temp is 0.612537543375436
    ## Next temp is 0.71554127394506
    ## Next temp is 0.820549552599721
    ## Next temp is 0.932601196819375
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-19.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-20.png)<!-- -->

    ## Next temp is 0.00215564912824267
    ## Next temp is 0.00806369699989309
    ## Next temp is 0.0345663888083769
    ## Next temp is 0.10527902415581
    ## Next temp is 0.190815060525455
    ## Next temp is 0.301676732928779
    ## Next temp is 0.408019381564265
    ## Next temp is 0.511757977500771
    ## Next temp is 0.612930373127768
    ## Next temp is 0.713926179083176
    ## Next temp is 0.833277615375933
    ## Next temp is 0.946654490600338
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-21.png)<!-- -->

    ## Next temp is 0.0773547295654063
    ## Next temp is 0.156951751093035
    ## Next temp is 0.246743934816546
    ## Next temp is 0.350351655257319
    ## Next temp is 0.451284331487749
    ## Next temp is 0.565519893096129
    ## Next temp is 0.666482110226188
    ## Next temp is 0.785222400559906
    ## Next temp is 0.884347152052067
    ## Next temp is 0.983438817982916
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-22.png)<!-- -->

    ## Next temp is 0.00529407775399103
    ## Next temp is 0.0183886888820913
    ## Next temp is 0.0492195488052411
    ## Next temp is 0.102632099289143
    ## Next temp is 0.18020388013673
    ## Next temp is 0.266213243138171
    ## Next temp is 0.364199671433694
    ## Next temp is 0.46222858605499
    ## Next temp is 0.545177424450994
    ## Next temp is 0.652147767767227
    ## Next temp is 0.753366116505842
    ## Next temp is 0.859840541375873
    ## Next temp is 0.962342480907915
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-10-23.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-10-24.png)<!-- -->

All EKI models do well at estimating the mean parameter. There is no
change in how the noise parameter is being estimated across a range of
$\alpha$ values, which is not surprising since the covariance matrices
used to update the particles should not depend on $\alpha$.

## Changing $\sigma^2$

In this experiment I change the true value of $\sigma$ from $2$ to $10$
in increments of $2$. The purpose of this is to determine whether the
EKI algorithm can detect values of $\sigma$ which are further away from
the prior. Additionally, we would expect the $\alpha$ particles to
become more diffuse as the true value of $\sigma$ increases.

``` r
sigma_seq <- seq(0.5, 8.5, length.out = 5)
num_dimensions <- 50
num_particles <- 400
for (sigma in sigma_seq) {

  true_parameters <- list(alpha = 2, sigma = sigma, x = rep(1, num_dimensions))
  generate_results(num_particles, true_parameters)
}
```

    ## Next temp is 0.00082867134976852
    ## Next temp is 0.00764028562359403
    ## Next temp is 0.042719298847722
    ## Next temp is 0.122308715924309
    ## Next temp is 0.215272935032492
    ## Next temp is 0.316814509215657
    ## Next temp is 0.43453893849406
    ## Next temp is 0.538274838465655
    ## Next temp is 0.651779990773881
    ## Next temp is 0.756603724395652
    ## Next temp is 0.860964370083573
    ## Next temp is 0.979500033902429
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

    ## Next temp is 0.0694174118545819
    ## Next temp is 0.126527039137796
    ## Next temp is 0.177456776979639
    ## Next temp is 0.225111073762815
    ## Next temp is 0.28600244965704
    ## Next temp is 0.349816726803175
    ## Next temp is 0.43059541578018
    ## Next temp is 0.509594883167636
    ## Next temp is 0.590284253354129
    ## Next temp is 0.673587046577448
    ## Next temp is 0.768536238784417
    ## Next temp is 0.867375566765372
    ## Next temp is 0.967278287887686
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

    ## Next temp is 0.0247427789708974
    ## Next temp is 0.0680083448930907
    ## Next temp is 0.109587650467597
    ## Next temp is 0.151707733184753
    ## Next temp is 0.199998747524172
    ## Next temp is 0.252258951465052
    ## Next temp is 0.310574932267493
    ## Next temp is 0.369678209388307
    ## Next temp is 0.437257119295541
    ## Next temp is 0.50966371206254
    ## Next temp is 0.59113093543701
    ## Next temp is 0.66880837029259
    ## Next temp is 0.757462261908549
    ## Next temp is 0.840642788820035
    ## Next temp is 0.939374699572812
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-3.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-4.png)<!-- -->

    ## Next temp is 0.020554223027983
    ## Next temp is 0.089533123710322
    ## Next temp is 0.18316709646175
    ## Next temp is 0.28732808278078
    ## Next temp is 0.401971668395659
    ## Next temp is 0.512702989235894
    ## Next temp is 0.617815827456568
    ## Next temp is 0.720078270363608
    ## Next temp is 0.836748322896018
    ## Next temp is 0.953776730898824
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-5.png)<!-- -->

    ## Next temp is 0.100841487810705
    ## Next temp is 0.195253658250371
    ## Next temp is 0.286898509565689
    ## Next temp is 0.383610711938137
    ## Next temp is 0.484770277407867
    ## Next temp is 0.58909261844995
    ## Next temp is 0.699746435633655
    ## Next temp is 0.801341388641088
    ## Next temp is 0.907382706112557
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-6.png)<!-- -->

    ## Next temp is 0.0254754343337749
    ## Next temp is 0.0840281243057291
    ## Next temp is 0.152640814889358
    ## Next temp is 0.233300376129105
    ## Next temp is 0.327943663847907
    ## Next temp is 0.427422249580689
    ## Next temp is 0.528910827651617
    ## Next temp is 0.626850397533231
    ## Next temp is 0.724707895762919
    ## Next temp is 0.841194021733056
    ## Next temp is 0.952767262172663
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-7.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-8.png)<!-- -->

    ## Next temp is 0.0496334084261725
    ## Next temp is 0.132815927169721
    ## Next temp is 0.221013943623984
    ## Next temp is 0.329958465518383
    ## Next temp is 0.431116491460217
    ## Next temp is 0.544099907819745
    ## Next temp is 0.653904454972175
    ## Next temp is 0.760542839900727
    ## Next temp is 0.856979605299311
    ## Next temp is 0.964191724569994
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-9.png)<!-- -->

    ## Next temp is 0.0935976869337699
    ## Next temp is 0.185222567694457
    ## Next temp is 0.283922169137732
    ## Next temp is 0.373976564486796
    ## Next temp is 0.478429655256771
    ## Next temp is 0.581463934075035
    ## Next temp is 0.680216021253695
    ## Next temp is 0.777057313659798
    ## Next temp is 0.884114706071031
    ## Next temp is 0.983924095925776
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-10.png)<!-- -->

    ## Next temp is 0.0250884761537722
    ## Next temp is 0.0882161334715146
    ## Next temp is 0.169645794816616
    ## Next temp is 0.258607176060574
    ## Next temp is 0.358763903415275
    ## Next temp is 0.457159113032379
    ## Next temp is 0.562966288295024
    ## Next temp is 0.656930038821103
    ## Next temp is 0.762444670864912
    ## Next temp is 0.854766164269854
    ## Next temp is 0.953867029821043
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-11.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-12.png)<!-- -->

    ## Next temp is 0.0620537615640808
    ## Next temp is 0.164650578926143
    ## Next temp is 0.257927823437086
    ## Next temp is 0.355526602975778
    ## Next temp is 0.457429677154315
    ## Next temp is 0.573711872542556
    ## Next temp is 0.68244558418426
    ## Next temp is 0.793868302870657
    ## Next temp is 0.89743901919454
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-13.png)<!-- -->

    ## Next temp is 0.0450259539570422
    ## Next temp is 0.118060875265514
    ## Next temp is 0.197827535100775
    ## Next temp is 0.282085907966012
    ## Next temp is 0.379448298158458
    ## Next temp is 0.484263613804344
    ## Next temp is 0.596517253534706
    ## Next temp is 0.708583414654557
    ## Next temp is 0.813116957085503
    ## Next temp is 0.910004682528542
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-14.png)<!-- -->

    ## Next temp is 0.0214713177956659
    ## Next temp is 0.0683700150943348
    ## Next temp is 0.128491013773246
    ## Next temp is 0.212363268503772
    ## Next temp is 0.292363078360792
    ## Next temp is 0.373791698255395
    ## Next temp is 0.474381992283124
    ## Next temp is 0.572346418611092
    ## Next temp is 0.678087670807492
    ## Next temp is 0.781691102763251
    ## Next temp is 0.893229550132672
    ## Next temp is 0.992231011000768
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-15.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-16.png)<!-- -->

    ## Next temp is 0.0713859939242031
    ## Next temp is 0.167774857710307
    ## Next temp is 0.261411117845757
    ## Next temp is 0.363494464973717
    ## Next temp is 0.474108014136815
    ## Next temp is 0.579555633182481
    ## Next temp is 0.679523767725672
    ## Next temp is 0.773390446481057
    ## Next temp is 0.868901417071101
    ## Next temp is 0.96460076016485
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-17.png)<!-- -->

    ## Next temp is 0.0233117258784618
    ## Next temp is 0.0767275215921766
    ## Next temp is 0.14142287150722
    ## Next temp is 0.228026471511724
    ## Next temp is 0.325155675116057
    ## Next temp is 0.411233300035386
    ## Next temp is 0.513257633033009
    ## Next temp is 0.619586796288997
    ## Next temp is 0.724093707176862
    ## Next temp is 0.830387899636195
    ## Next temp is 0.947494020496031
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-18.png)<!-- -->

    ## Next temp is 0.0179885915440214
    ## Next temp is 0.0487213323906699
    ## Next temp is 0.0993065604863536
    ## Next temp is 0.176075062656317
    ## Next temp is 0.256250272663043
    ## Next temp is 0.352992762498744
    ## Next temp is 0.454272573269791
    ## Next temp is 0.55039910517751
    ## Next temp is 0.651422344437804
    ## Next temp is 0.759710313839266
    ## Next temp is 0.859212654700572
    ## Next temp is 0.960662173025128
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-11-19.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-11-20.png)<!-- -->

We see that the EKI method is able to estimate $\sigma^2$ quite well,
even as it gets larger. The particles for $\alpha$ are more diffuse as
$\sigma^2$ increases, which is what we naturally expect.

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

    ## Next temp is 0.00928181088956578
    ## Next temp is 0.0409538762330974
    ## Next temp is 0.12468943051649
    ## Next temp is 0.196457640538734
    ## Next temp is 0.308787149713333
    ## Next temp is 0.38795129510286
    ## Next temp is 0.473059357766192
    ## Next temp is 0.567967455930458
    ## Next temp is 0.655520910498937
    ## Next temp is 0.759775062357178
    ## Next temp is 0.836328519462736
    ## Next temp is 0.904887638173547
    ## Next temp is 0.990423319692113
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

    ## Next temp is 0.0834310385295366
    ## Next temp is 0.153974229931592
    ## Next temp is 0.225181805660791
    ## Next temp is 0.302787983669729
    ## Next temp is 0.405396336157892
    ## Next temp is 0.500346675662361
    ## Next temp is 0.616026422456751
    ## Next temp is 0.699937938314315
    ## Next temp is 0.789030825911357
    ## Next temp is 0.878389619521034
    ## Next temp is 0.941121895158802
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-2.png)<!-- -->

    ## Next temp is 0.0116467677499403
    ## Next temp is 0.0491279587617451
    ## Next temp is 0.092349785300649
    ## Next temp is 0.147562436176495
    ## Next temp is 0.210368723358359
    ## Next temp is 0.282655388450917
    ## Next temp is 0.35631973247712
    ## Next temp is 0.433195297292213
    ## Next temp is 0.498945659609137
    ## Next temp is 0.583509645855736
    ## Next temp is 0.653555421196374
    ## Next temp is 0.739342070530022
    ## Next temp is 0.822968683126579
    ## Next temp is 0.904592210849247
    ## Next temp is 0.98076005166997
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-3.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-4.png)<!-- -->

    ## Next temp is 0.0153224769125928
    ## Next temp is 0.072253923733641
    ## Next temp is 0.158942449409779
    ## Next temp is 0.250115801968081
    ## Next temp is 0.357642933943958
    ## Next temp is 0.467128091165393
    ## Next temp is 0.576174716977655
    ## Next temp is 0.691956395706913
    ## Next temp is 0.793464184854212
    ## Next temp is 0.899251218356913
    ## Next temp is 0.99523596867421
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-5.png)<!-- -->

    ## Next temp is 0.106442777495298
    ## Next temp is 0.201645619293658
    ## Next temp is 0.298364994468934
    ## Next temp is 0.391955468034737
    ## Next temp is 0.500067734408486
    ## Next temp is 0.612923977485685
    ## Next temp is 0.714691845811126
    ## Next temp is 0.816024124015301
    ## Next temp is 0.918632091344816
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-6.png)<!-- -->

    ## Next temp is 0.025304526540023
    ## Next temp is 0.0840728129971337
    ## Next temp is 0.151318940895822
    ## Next temp is 0.225133876932338
    ## Next temp is 0.316648213549853
    ## Next temp is 0.417591314697598
    ## Next temp is 0.513495306791508
    ## Next temp is 0.620268997024166
    ## Next temp is 0.717381710928101
    ## Next temp is 0.815209043423934
    ## Next temp is 0.915991351356513
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-7.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-8.png)<!-- -->

    ## Next temp is 0.0141808178992063
    ## Next temp is 0.0740231824025985
    ## Next temp is 0.168083910384471
    ## Next temp is 0.268021893065288
    ## Next temp is 0.379799541882456
    ## Next temp is 0.478149139026803
    ## Next temp is 0.578232622347829
    ## Next temp is 0.683434073302406
    ## Next temp is 0.791164497246618
    ## Next temp is 0.896935457975982
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-9.png)<!-- -->

    ## Next temp is 0.114383491256156
    ## Next temp is 0.215115147817456
    ## Next temp is 0.318685284054495
    ## Next temp is 0.417996702067478
    ## Next temp is 0.517525373843158
    ## Next temp is 0.618969612831789
    ## Next temp is 0.731821288425834
    ## Next temp is 0.835321332330358
    ## Next temp is 0.938483547484928
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-10.png)<!-- -->

    ## Next temp is 0.0279387548651342
    ## Next temp is 0.0862216846315664
    ## Next temp is 0.16291681577622
    ## Next temp is 0.248918349986221
    ## Next temp is 0.341975955413097
    ## Next temp is 0.441563223151467
    ## Next temp is 0.548314692836748
    ## Next temp is 0.652320046996621
    ## Next temp is 0.761551374186337
    ## Next temp is 0.860561979820929
    ## Next temp is 0.970435674410568
    ## Next temp is 1

![](results_normal_files/figure-gfm/unnamed-chunk-12-11.png)<!-- -->![](results_normal_files/figure-gfm/unnamed-chunk-12-12.png)<!-- -->

Indeed we find that this is the case.

## Conclusion

To summarize, GEKI does well at estimating both parameters when we use
summary statistics. This is interesting considering it doesn’t do well
when using the full data.
