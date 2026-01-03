library("brms")
library("dplyr")
library("tidyr")
library("cmdstanr")

als = read.csv("../Data/PROACT_ALSFRS.csv")
muscle = read.csv("../Data/PROACT_MUSCLESTRENGTH.csv")

als = als %>%
  mutate(
    Q4_Handwriting = as.integer(Q4_Handwriting)
  )

muscle = muscle %>%
  filter(
    Test_Location == "WRIST JOINT",
    Test_trial == 1,
    Test_Laterality == "RIGHT"
  )

data = als %>%
  left_join(
    muscle,
    by=c("subject_id", "ALSFRS_Delta" = "MS_Delta")
  ) %>%
  select(
    subject_id, Q4_Handwriting, Test_Result, ALSFRS_Delta
  ) %>%
  filter(
    !is.na(Q4_Handwriting),
    !is.na(Test_Result),
    !is.na(subject_id),
    !is.na(ALSFRS_Delta)
  ) %>% mutate(
    Q4_Handwriting = ordered(Q4_Handwriting),
    Time_yrs = ALSFRS_Delta / 365
  )


longitudinal_model = brm(
  formula = Q4_Handwriting ~ Test_Result + Time_yrs + (1 | subject_id),
  data = data,
  seed = 42,
  family = cumulative("probit"),
  prior = prior(exponential(1), class = "sd"),
  backend = "cmdstanr",
  cores = parallel::detectCores()
)

summary(longitudinal_model)

ps = posterior_summary(longitudinal_model, digits = 5)
ps[c("Intercept[1]", "Intercept[2]", "Intercept[3]", "Intercept[4]",
     "b_Test_Result", "b_Time_yrs", "sd_subject_id__Intercept"), ]
