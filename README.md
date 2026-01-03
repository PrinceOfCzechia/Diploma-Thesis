This repository stores the (mostly R) code accompanying my master's thesis.
The prerequisites are the following R packages:
+ dplyr
+ cmdstanr
+ brms
+ mvtnorm

and the data from the PRO-ACT database.
These cannot be shared by me, but anyone can request access at the database
[website](https://ncri1.partners.org/ProACT/Home/Index).

A short description of the files in this repository follows.
The files are listed in that order, in which they appear in the thesis body.

#### pro_act_brms_univariate.R
A first example of an actual cumulative model fit, using cross-sectional data
from PRO-ACT.
There is a single univariate covariate, a univariate response, and the fit is
obtained by the *brms* library.

The result is shown in section 2.5.

#### pro_act_brms_longitudinal.R
In sections 1.4 and 2.4, we discuss mixed-effects models and their use in the
cumulative model framework.

This code produces the example fit we show in section 2.6.

#### ghk_probit_model.stan
Our implementation extends an existing one which can be found at
[Stan's development team github](https://github.com/stan-dev/example-models/blob/master/misc/multivariate-probit/probit-multi-good.stan).
For detailed description, see the thesis text, especially sections
3.2, 3.3.3, 4.1 and 4.2.

#### ghk_probit_synthetic.R
We test the model from *ghk_probit_model.stan* on 3D synthetic data where we
specify the slope which should be estimated by the linear model and the
correlation between components of the latent variable, also to be estimated.

The results are shown in section 4.3.

#### pro_act_univariate.R
This file uses the GHK-style model to mirror the univariate fit from
*pro_act_brms_univariate.R*.

We discuss the model summary and comparison to the *brms* fit from earlier in
section 4.4. 

#### pro_act_multivariate.R
Here, we fit the GHK-style in the following setting:
+ multivariate response
+ the linear predictor can have different effect on each of the response components
+ there is an unknown correlation structure between the response components

This is the true use case of the GHK-style implementation.
The results are discussed in section 4.5.

#### ghk_integration.R
We made a toy implementation of the GHK algorithm as means of Monte Carlo
integration, too.
It can be found in this file, and the result is compared to *mvtnorm* library
result of the same integral.
The code is discussed in Appendix C.1.