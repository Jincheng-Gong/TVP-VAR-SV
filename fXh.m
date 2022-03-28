%%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%
%%
%%  mXh = fXh(vyh, nk, na)
%%
%%  "fXh" sets matrix "X_hat" from 1*nk vector "y_hat"
%%
%%  [input]
%%    vyh:     1*nk vector
%%    nk, na:  scalar
%%
%%  [output]
%%    mXh:     nk*na matrix
%%

function mXh = fXh(vyh, nk, na)

mXh = zeros(nk, na);

for i = 2 : nk
   mXh(i, (i-1)*(i-2)/2+1 : i*(i-1)/2) = -vyh(1:i-1);
end
