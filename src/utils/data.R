source('C:/Users/owenj/OneDrive/Uni/Vacation Scholarship/GEKI_Vacation_Scholarship/src/models/eki_normal.R')

set.seed(2025)

true_parameters <- list(alpha = 2, sigma = 2, x = rep(1, 50))
test_data <- likelihood_normal(true_parameters)