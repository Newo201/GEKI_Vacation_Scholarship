source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/samples_normal.R')
source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/pdfs_normal.R')

true_params_1d = list(alpha = 2, sigma = 2, x = 1)
true_params_2d = list(alpha = 2, sigma = 2, x = c(1, 2))

# Test the dimensions of the likelihood samples

test_dimensions <- function(true_params, num_samples) {
  
  l_samples <- likelihood_sample(true_params, num_samples)
  return(as.numeric(dim(l_samples)))
}

test_that('Dimensions of likelihood samples are correct', 
          {
            expect_identical(test_dimensions(true_params_1d, 100), c(100, 1))
            expect_identical(test_dimensions(true_params_2d, 100), c(100, 2))
          })

# Test the dimension of the likelihood density
test_density_dimensions <- function(true_params, num_samples) {
  
  simulated_data <- likelihood_sample(true_params, num_samples)
  density <- loglike_pdf(simulated_data, true_params)
  return(length(density))
}

test_that('Dimensions of likelihood density are correct', 
          {
            expect_equal(test_density_dimensions(true_params_1d, 100), 100)
            expect_equal(test_density_dimensions(true_params_2d, 100), 100)
          })


