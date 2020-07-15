# Vittadello2020
MATLAB code for numerical simulations of the functional differential equation model in Vittadello et al. 2020.

In NumericalSoln_v1.m, large values of the shape parameter for the Erlang distribution may cause MATLAB to overflow to infinity. This is due to the direct calculation of the Erlang distribution.

In NumericalSoln_v2.m we remedy this situation by first calculating the logarithm of the Erlang distribution and then taking the exponential of the result. This allows for large values of the shape parameter, however the code is considerably slower to run than NumericalSoln_v1, so it is recommended to use the latter when possible. Note that NumericalSoln_v2.m requires the function file NumericalSoln_v2_RTErlang.m, which contains the right-truncated Erlang distribution.
