library("brms")
library("dplyr")
library("tidyr")
library("cmdstanr")

als = read.csv("../Data/PROACT_ALSFRS.csv")
muscle = read.csv("../Data/PROACT_MUSCLESTRENGTH.csv")

als_t0 = als %>%
  mutate(
    Q4_Handwriting = as.integer(Q4_Handwriting),
    Q5a_Cutting_without_Gastrostomy = as.integer(Q5a_Cutting_without_Gastrostomy),
    Q6_Dressing_and_Hygiene = as.integer(Q6_Dressing_and_Hygiene)
  ) %>%
  filter(
    ALSFRS_Delta == 0
  )

muscle_t0 = muscle %>%
  filter(
    Test_Location == "WRIST JOINT",
    Test_trial == 1,
    Test_Laterality == "RIGHT",
    MS_Delta == 0
  )

data_t0 = als_t0 %>%
  left_join(
    muscle_t0,
    by=c("subject_id")
  ) %>%
  select(
    subject_id, Q4_Handwriting, Q5a_Cutting_without_Gastrostomy,
    Q6_Dressing_and_Hygiene, Test_Result, ALSFRS_Delta
  ) %>%
  filter(
    !is.na(Q4_Handwriting),
    !is.na(Q5a_Cutting_without_Gastrostomy),
    !is.na(Q6_Dressing_and_Hygiene),
    !is.na(Test_Result)
  )

Y = cbind(
  as.integer(data_t0$Q4_Handwriting),
  as.integer(data_t0$Q5a_Cutting_without_Gastrostomy),
  as.integer(data_t0$Q6_Dressing_and_Hygiene)
) + 1L

X = data_t0$Test_Result

stan_data = list(
  N = nrow(data_t0),
  D = 3,
  K = 5,
  P = 1,
  Y = Y,
  X = matrix(X, ncol=1)
)

mod = cmdstan_model("ghk_probit_model.stan")

fit = mod$sample(
  data = stan_data,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  init = 0
)

print(
  fit$summary(
    variables = c(
      "beta",
      "tau",
      "Sigma[1,2]",
      "Sigma[1,3]",
      "Sigma[2,3]"),
    mean, ~quantile(.x, probs = c(0.025, 0.975))
  ),
  n = 18
)
