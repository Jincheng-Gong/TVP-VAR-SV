%%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%
%%
%%  [] = setvar(...)
%%
%%  "setvar" sets variables or options as global variable
%%
%%  [input]
%%   (stype, ...)
%%
%%   ('data', my, nlag) -> set data and lags
%%
%%   ('intercept', fli) -> set time-varying intercept
%%      0: off, demean (default)
%%      1: on
%%
%%   ('SigB', flSb) -> set non-diagonal matrix for Sigma_beta
%%      0: diagonal (default)
%%      1: non-digaonl
%%
%%   ('impulse', nimp) -> set length of impulse response
%%                        (default: 12)
%%
%%   ('fastimp', flfi) -> fast computing of response
%%                        (default: 1)
%%
%%   ('ranseed', iseed) -> set ranseed (default: 1)
%%
%%   ('prior', spar, ., .) -> set prior
%%
%%  [global variables]
%%      m_my:     variable y_t (n*k matrix)
%%      m_asvar:  variable names
%%      m_nl:     # of lags (l, scalar)
%%      m_ns:     # of time periods (n, scalar)
%%      m_nk:     # of series (k, scalar)
%%      m_fli:    intercept type (flag)
%%      m_flSb:   non-diagonal Sig_b (flag)
%%      m_nimp:   # of impulse response
%%      m_flfi:   fast computing of response (flag)
%%      m_iseed:  ranseed
%%
%%      m_dvb0 m_dVb0 m_dva0 m_dVa0 m_dvh0 m_dVh0: priors
%%

function setvar(stype, arg1, arg2, arg3)

global m_my m_asvar m_nl m_ns m_nk m_fli m_flSb m_nimp m_flfi ...
       m_iseed m_dvb0 m_dVb0 m_dva0 m_dVa0 m_dvh0 m_dVh0

switch stype
    case 'data'
        m_my = arg1;
        m_asvar = arg2;
        m_nl = arg3;
        m_ns = size(m_my, 1);
        m_nk = size(m_my, 2);
    case 'intercept'
        m_fli = arg1;
    case 'SigB'
        m_flSb = arg1;
    case 'impulse'
        m_nimp = arg1 + 1;
    case 'fastimp'
        m_flfi = arg1;
    case 'ranseed'
        m_iseed = arg1;
    case 'prior'
      switch arg1
          case 'b'
            m_dvb0 = arg2;    % SetSigB(1): Sigb ~ IW(vb0, I*Vb0) 
            m_dVb0 = arg3;    % SetSigB(0): sigb_i^2 ~ IG(vb0, Vb0)
          case 'a'
            m_dva0 = arg2*2;  % siga_i^2 ~ IG(va0, Va0)
            m_dVa0 = arg3*2;
          case 'h'
            m_dvh0 = arg2*2;  % sigh_i^2 ~ IG(vh0, Vh0)
            m_dVh0 = arg3*2;
      end
end
