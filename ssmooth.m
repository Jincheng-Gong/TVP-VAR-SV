%%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%
%%
%%  malpha = ssmooth(my, amZ, amG2, mH2, va0, mH02)
%%
%%  "ssmooth" implements simulation smoother
%%  by de Jong & Shephard (1995)
%%
%%  [model]
%%      y_t = Z_t*alpha_t + G_t*e_t
%%      alpha_{t+1} = alpha_t + H_t*e_t
%%      e_t ~ N(0, I)
%%      G_t*H_t' = O
%%
%%      y_t:     nk*1 vector
%%      Z_t:     nk*np matrix
%%      alpha_t: np*1 vector
%%      G_t:     nk*(nk+np) matrix
%%      H_t:     np*(nk+np) matrix
%%
%%  [input]
%%      my:     response (ns*nk vector)
%%      amZ:    independent variable (nk*np*ns array)
%%      amG2:   G_t*G_t' (nk*nk*ns array)
%%      mH2:    H_t*H_t' = H*H' (np*np matrix)
%%      va0:    alpha_0 (np*1 vector)
%%      mH02:   H_0*H_0' (np*np matrix)
%%
%%  [output]
%%      malpha:  sampled state variable (np*ns matrix)
%%

function malpha = ssmooth(my, amZ, amG2, mH2, va0, mH02)

%%--- set variables ---%%

ns = size(my, 1);    % # of time periods
nk = size(my, 2);    % # of series
np = size(mH2, 1);    % # of state

va = va0;
mP = mH02;
vr = zeros(np, 1);
mU = zeros(np);

me = zeros(nk, ns);
amDinv = zeros(nk, nk, ns);
amL = zeros(np, np, ns);
meta = zeros(np, ns);

%--- Kalman filter ---%%

for i = 1 : ns
    me(:, i) = my(i, :)' - amZ(:, :, i) * va;

    mD = amZ(:,:,i) * mP * amZ(:,:,i)' + amG2(:,:,i);
    drD = rcond(mD);
    if isnan(drD) || (drD < eps*10^2)
        mDinv = eye(nk) * 10;
    else
        mDinv = inv(mD);
    end
    amDinv(:, :, i) = mDinv;

    mK = mP * amZ(:, :, i)' * mDinv;
    amL(:, :, i) = eye(np) - mK * amZ(:, :, i);
    
    va = va + mK * me(:, i);
    mP = mP * amL(:, :, i)' + mH2;
end


%%--- simulation smoother ---%%

i = ns;
while i >= 1
    mC = mH2 - mH2 * mU * mH2;
    mC = (mC + mC')/2;
    [mCc, fl] = chol(mC, 'lower');
    drC = rcond(mC);
    if fl > 0
        mCc = eye(np) * 0.01;
        mCinv = eye(np) * 10^4;
    elseif isnan(drC) || (drC < eps*10^2)
        mCinv = eye(np) * 10^4;
    else
        mCinv = inv(mC);
    end
    
    veps = mCc * randn(np, 1);
    meta(:, i) = mH2 * vr + veps;
    mV = mH2 * mU * amL(:, :, i);

    vr = amZ(:,:,i)' * amDinv(:,:,i) * me(:, i) ...
       + amL(:,:,i)' * vr - mV' * mCinv * veps;
    mU = amZ(:,:,i)' * amDinv(:,:,i) * amZ(:,:,i) ...
       + amL(:,:,i)' * mU * amL(:,:,i) + mV' * mCinv * mV;

    i = i - 1;
end

mC = mH02 - mH02 * mU * mH02;
mC = (mC + mC')/2;
[mCc, fl] = chol(mC, 'lower');
if fl > 0
    mCc = eye(np) * 0.07;
end
veta0 = mH02 * vr + mCc *  randn(np, 1);

malpha = zeros(np, ns);
malpha(:, 1) = va0 + veta0;
for i = 1 : ns-1
    malpha(:, i+1) = malpha(:, i) + meta(:, i);
end
