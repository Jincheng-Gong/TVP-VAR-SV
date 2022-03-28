%%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%
%%
%%  mAi = finvm(mA)
%%
%%  "finvm" outputs inv(A)
%%  when singular, Moore-Penrose pseudoinverse of matrix
%%
%%  [input]
%%    mA:      nk*nk matrix
%%
%%  [output]
%%    mAi:     nk*nk matrix
%%

function mAi = finvm(mA)

drA = rcond(mA);

if isnan(drA) || (drA < eps*10^2)
    mAi = inv(diag(diag(mA)));
else
    mAi = inv(mA);
end