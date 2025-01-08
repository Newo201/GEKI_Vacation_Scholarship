# Generalised Ensemble Kalman Inversion

## By Owen Jackson

### Supervised by Dr Imke Botha and Professor Jennifer Flegg

### University of Melbourne Vacation Scholarship

## Background

Ensemble Kalman techniques are a class of Monte Carlo algorithms that can be used in high dimensional settings with the restriction of additive Gaussian likelihoods where the noise parameter is required to be known.

For example, Ensemble Kalman Filters (EnKF), as the name suggests uses an ensemble of Kalman filters where the Kalman gain is calculated based on empirical covariance matrices. Ensemble Kalman Inversion adds a tempering step to this process. **Probably needs more explanation, some references and some equations.**

On the other hand, Approximate Bayesian Computation is a collection of algorithms which can be used on general likelihoods that cannot be evaluated pointwise. However, due to the likelihoods being intractable, they have only been effective in low dimensional settings.

Duffield and Singh (2022) proposed a Generalised Ensemble Kalman Inversion algorithm which makes the best of both worlds. It is able to operate in high dimensional settings, but it can also operate on general likelihoods.

## Purpose

The purpose of this research is to implement the GEKI algorithm (Duffield and Singh, 2022) and to test its performance in two situations:

-   A simulated multivariate normal model where we seek to estimate the mean and variance

-   A malaria transmission model

## GEKI Algorithm

**Inputs:** a sequence of inverse temperatures $\{\lambda_L\}^L_{l=0}$, the prior distributions for each of the parameters, and a single draw from the likelihood using the true parameters

-   The observation may be mutlidimensional (i.e. a vector, or a realization of a time series) but we cannot have more than one observation. Mathematically the dimensions of $y$ are $d_y \times 1$ (a column vector)

-   In simulated settings we know the true parameters, but in real world settings the true parameters are unknown and we are only given the observation

1.  Sample a pre-determined number of particles $N$ from the prior distribution.

    -   This represents our initial ensemble where each particle represents a sample of the parameters

2.  For $l = 1 \ldots n$ do:

3.  For each particle, make a single draw from the likelihood using the sampled parameters. $y_{l-1}^{(i)} \sim p(\cdot|x_{l-1}^{(i)})$

4.  Form sample covariance and cross covariance matrices

    -   Note the following dimensions of the covariance matrices

        $C^{yy}: d_y \times d_y$

        $C^{xx}: d_x \times d_x$

        $C^{xy}: d_x \times d_y$

        $C^{y|x}: d_y \times d_y$

5.  Set step size $h_l = \lambda_l - \lambda_{l-1}$

6.  Generate observation perturbations $\eta_l^{(i)} \sim N(0, (h_l^{-1}-1)C_{l-1}^{y|x})$

7.  Move particles according to the following update:

    $x_l^{(i)} = x_{l-1}^{(i)} + C_{l-1}^{xy}(C_{l-1}^{yy} + (h_l^{-1}-1)C_{l-1}^{y|x})^{-1}(y-y_{l-1}^{(i)}-\eta_l^{(i)})$

8.  Return the final ensemble $\{x_L^{(i)}\}$

### Selecting Temperature Adaptively

Instead of providing a sequence of inverse temperatures as the input, we can instead choose to select the temperature adaptively.

The idea is to use pseudo weights $w_l^{(i)} \propto exp(-\frac{1}{2}(\lambda_l-\lambda_{l-1})(y-y_{l-1}^{(i)})^T(C_{l-1}^{y|x})^{-1}(y-y_{l-1}^{(i)}))$.

Given $\lambda_{l-1}$ we want to choose $\lambda_l$ such that the $ESS = \rho N$ where $\rho$ is a hyper-parameter (set to $\frac{1}{2}$). We can find this using any standard root finding algorithm; I use the `uniroot` function in R. In cases where the $\lambda_l > 1$ we set $\lambda_l = 1$ and make this our final iteration.

## References
