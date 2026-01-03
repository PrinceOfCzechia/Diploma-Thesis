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


###
# model using brms
###

brms_model = brm(
  formula = Q4_Handwriting ~ Test_Result,
  data = data_t0,
  seed = 42,
  family = cumulative("probit"),
  backend = "cmdstanr",
  cores = parallel::detectCores()
)
summary(brms_model)
