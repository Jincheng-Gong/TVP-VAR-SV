%%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%
%%
%%  dv = ftsvar(vx, iBm)
%%
%%  "ftsvar" computes time-series variance
%%  for certain bandwidth with Parzen window
%%
%%  [input]
%%    vx:   data series (ns*1 vector)
%%    iBm:  bandwidth (scalar)
%%
%%  [output]
%%    dvar:   estimated variance
%%

function dvar = ftsvar(vx, iBm)

ns = length(vx);
vi = (1 : iBm)' / iBm;
vK1 = 1 - 6 * vi.^2 .* (1 - vi);
vK2 = 2 * (1 - vi).^3;

vxm = vx - mean(vx);

dvar = 0;
for i = 1 : iBm
    davar = vxm(1:end-i)' * vxm(i+1:end) / ns;
    if vi(i) < 0.5
        dvar = dvar + vK1(i) * davar;
    else
        dvar = dvar + vK2(i) * davar;
    end
end

dvar = sum(vxm.^2) / ns + 2 * ns / (ns - 1) * dvar;

