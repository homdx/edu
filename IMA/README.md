# IMA-M2IBC

## Statistics
https://www.statlect.com/

The Rmd files in the Statistics folder are incomplete. Take a look at the GitHub repo for complete files.

When we are talking about confidence intervals, we are modeling the mean.
When we are talking about prediction intervals, we are modeling the individual.


## Logistic regression

Why do we use log(odds) of an outcome instead of probability?

Because p in [0, 1], so eventually the linear model will move out of the bounds.
prob in (0, 1] => odds in (0, infty] => log(odds) in (-infty, infty)

The regression coeeficients are on a log(odds) scale, but the fitted values are on a probability scale.
