%%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%
%%
%%  mXt = fXt(my, fli)
%%
%%  "fXt" sets matrix "X_t" from y(t-1)...y(t-p)
%%
%%  [input]
%%    myi:  data (nl*nk matrix)
%%    fli:  intercept flag
%%
%%  [output]
%%    mXt:  nk*nb matrix
%%

function mXt = fXt(myi, fli)

myif = flipud(myi)';
if fli == 1
    vyi = [1 myif(:)'];
else
    vyi = myif(:)';
end

mXt = kron(eye(size(myi, 2)), vyi);
