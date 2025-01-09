## Sampling

### Implementation

-   Sampling parameters from normal prior distributions; $\alpha \sim N(0, 5^2)$, $log(\sigma^2) \sim N(0, 2^2)$

-   Sampling the normal likelihood from a multivariate normal distribution $N(\alpha x, \sigma^2I)$ where $x$ represents a known vector

### Diagnostics

-   Plotting a histogram of prior and likelihood distributions to check that it closely resembles theoretical distribution

### Unit Tests

**Prior Samples:** `tests/test_priors.R`

-   That the dimensions of the prior samples are as expected

-   That the mean and standard deviation of the prior samples are close to the true mean and standard deviation

-   That different samples from the prior distribution produce a consistent mean

**Likelihood Samples:** `tests/test_likelihoods.R`

-   Test that the dimensions of the likelihood samples are correct (single and multiple dimensions)

-   Test that the dimensions of the likelihood density are correct (single and multiple dimensions)

## Densities

### Implementation

-   Calculate the densities for normal prior distributions

-   Calculate the densities for normal likelihood distribution

### **Diagnostics**

### **Unit Tests**

## MCMC For Normal Distribution

### Implementation

-   Use the `mcmc` package

### Diagnostics

-   Trace plots

-   Histograms

### Unit Tests

## **GEKI Algorithm For Normal Distribution**

### **Implementation**

**Standard Algorithm**

-   Draw a single sample from the likelihood using the true parameters

-   Initialize particles by sampling from prior distributions

-   For each particle draw from the likelihood

-   Calculate covariance matrices

-   Generate perturbations

-   Update particles

**Adaptive Temperature**

-   Helper functions to calculate ESS and find the next temperature

-   Integrate into the EKI algorithm

### Diagnostics

**Standard Algorithm**

-   Plot histograms for marginal posterior distributions of parameters

    -   Across different number of dimensions

    -   Using different parameter values

    -   Using adaptive and non-adaptive temperature sequence

**Adaptive Temperature**

-   Plot the temperature sequence to determine if it is increasing

### Unit Tests

**Dimensions:** `tests/test_eki_normal.R`

-   Dimensions of generated particles are correct

-   Dimensions of updated particles are correct

-   Dimensions of perturbations are correct

-   Dimension of particle pseudo-weights are correct

**Covariance Matrices:** `tests/test_eki_normal.R`

-   Dimensions of covariance matrices are correct

-   Covariance matrices are symmetric and positive definite

**Adaptive Temperature:** `tests/test_eki_normal.R`

-   Effective sample size is between 1 and N

-   Next selected temperature is between current temp and 1

-   Next selected temperature is chosen to match the target ESS

### Results

-   Change $\alpha$

-   Change number of dimensions

-   Change $\sigma^2$

-   Change $x$
