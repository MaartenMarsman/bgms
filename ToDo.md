# sampleMRF
We should check the new version of sampleMRF

# bgm.panel.r
## Implement na.action
We should create a function to compute and simulate from the full-conditional

## Implement threshold_model = "constant" (stationarity assumption)
We have a first implementation of the c++ code, but we need to adjust the input
from R in bgm.panel.r. We also still need to check for the consistency with the
reformat_data function.

## Count first time point as one, not zero.
I have made the change in the analysis input, but I need to adjust some
additional functionalities in bgm.panel.r, and reformat_data. It would also be a
good idea to check for the consistency with the sampleMRF function.

## Name output
We should incorporate the changes we made to the analysis output of bgm() to the
panel model output. It would also be good to offer the output of the thresholds
in case of a collapsed category in a more meaningful manner.

# C++ functions
## Restmatrix
We need to implement a restmatrix approach similar to the c++ code for bgm().
