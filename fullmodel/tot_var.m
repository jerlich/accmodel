function [sigma2_tot var_a var_s] = tot_var(sigma2_a, sigma2_s, nclicks, dt)

var_a = sigma2_a*dt;
var_s = nclicks*sigma2_s/41;

sigma2_tot = var_a + var_s;