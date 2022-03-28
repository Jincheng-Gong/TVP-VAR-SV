%%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%
%%
%%  dv = fGeweke(vx, iBm)
%%
%%  "fGeweke" computes Geweke statistics
%%  for convergence check of MCMC iteration
%%
%%  [input]
%%    vx:   data series (ns*1 vector)
%%    iBm:  bandwidth (scalar)
%%
%%  [output]
%%    dpv:  P-value for convergence
%%

function dpv = fGeweke(vx, iBm)

ns = length(vx);
n1 = floor(ns * 0.1);
n2 = floor(ns * 0.5);
vx1 = vx(1:n1);
vx2 = vx(n2+1:end);

dm1 = mean(vx1);
dm2 = mean(vx2);
dv1 = ftsvar(vx1, iBm);
dv2 = ftsvar(vx2, iBm);

dz = (dm1 - dm2) / sqrt(dv1 / n1 + dv2 / n2);

dpv = 2 * (1 - normcdf(abs(dz)));
