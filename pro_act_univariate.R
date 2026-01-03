library("brms")
library("dplyr")
library("cmdstanr")

als = read.csv("../Data/PROACT_ALSFRS.csv")
muscle = read.csv("../Data/PROACT_MUSCLESTRENGTH.csv")

als = als %>%
  mutate(
    Q4_Handwriting = as.integer(Q4_Handwriting)
  )

als_t0 = als %>%
  filter(ALSFRS_Delta == 0)

muscle_t0 = muscle %>%
  filter(
    MS_Delta == 0,
    Test_Location == "WRIST JOINT",
    Test_trial == 1,
    Test_Laterality == "RIGHT"
  )

data_t0 = als_t0 %>%
  left_join(
    muscle_t0, by="subject_id"
  ) %>%
  select(
    subject_id, Q4_Handwriting, Test_Result
  ) %>%
  filter(
    !is.na(Q4_Handwriting),
    !is.na(Test_Result)
  ) %>% mutate(
    Q4_Handwriting = ordered(Q4_Handwriting)
  )

stan_data = list(
  N = nrow(data_t0),
  D = 1,
  K = 5,
  P = 1,
  Y = matrix(as.integer(data_t0$Q4_Handwriting), ncol=1),
  X = matrix(data_t0$Test_Result, ncol=1)
)

ghk_mod = cmdstan_model("ghk_probit_model.stan")

ghk_fit = ghk_mod$sample(
  data = stan_data,
  seed = 42,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

print(
  ghk_fit$summary(
    variables = c(
      "beta",
      "tau"),
    mean, ~quantile(.x, probs = c(0.025, 0.975))
  ),
  n = 15
)
