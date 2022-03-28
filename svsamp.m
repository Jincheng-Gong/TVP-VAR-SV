%%--------------------------------------------------------%%
%%                    TVP-VAR package                     %%
%%--------------------------------------------------------%%
%%
%%  vhs = svsamp(vy, vh, dsig2, dh00, dsig02, nK)
%%
%%  "svsamp" implements multi-move sampler for SV model
%%  by Shephard & Pitt (1997) and Omori & Watanabe (2004)
%%
%%  [model]
%%    y_t = exp(h_t/2)*eps_t,  eps_t ~ N(0, 1)
%%    h_{t+1} = h_t + ets_t,  eta_t ~ N(0, sig^2)
%%    h_0 = h00,  eta_0 ~ N(0, sig0^2)
%%
%%  [input]
%%      vy:     response (ns*1 vector)
%%      vh:     current point of h (ns*1 vector)
%%      dsig2, dsig02:
%%              parameter (scalar)
%%      nK:     # of blocks for multi-move sampler
%%
%%  [output]
%%      vhs:  sampled stochastic volatility (ns*1 vector)
%%

function vhs = svsamp(vy, vh, dsig2, dh00, dsig02, nK)

global m_iachar m_iachmh m_k


%%--- set variables ---%%

ns = size(vy, 1);
nite = 5;       % # of iteration for h_hat

vhs = vh;       % current point

vk = [1 2];
while sum(diff(vk)<2) > 0
  vk = [1; floor(ns * ((1:nK)+rand(1, nK)) / (nK+2))'...
         ; ns+1];    % stochastic knots                
end
                
%%--- sampling start ---%%

for i = 1 : nK+1
  ir  = vk(i);
  id  = vk(i+1) - vk(i);
  ird = vk(i+1) - 1;

  vyi = vy(ir:ird);
  vho = vh(ir:ird);   % current (old) point
  vhn = zeros(id, 1); % new point

  if i <= nK
      dhrd1 = vh(ird+1);  % h(r+d+1)
  end
    
 %%--- finding model & draw candidate ---%
    
  for j = 1 : nite+1
        
      if j == 1
        vhh = vho;  % h_hat
      else
        vhh = vhn;
      end
        
      vgder2 = -0.5 * vyi.^2 ./ exp(vhh); % g''(h)
      vgder1 = -0.5 - vgder2;             % g'(h)
      vsiga2 = -1 ./ vgder2;              % sig2_ast
      vha = vhh + vsiga2 .* vgder1;       % h_ast
        
      if i <= nK
        vsiga2(id) = 1 / (-vgder2(id) + 1/dsig2);
        vha(id) = vsiga2(id) ...
                * (vgder1(id) - vgder2(id)*vhh(id) ...
                   + dhrd1/dsig2);             
      end

     %%--- simulation smoother ---%%

      if i == 1
        dh0 = dh00;
        dH20 = dsig02;
      else
        dh0 = vhs(ir-1);
        dH20 = dsig2;
      end
          
      da = dh0;
      dP = dH20;
      dH2 = dsig2;
      ve = zeros(id, 1);
      vDinv = zeros(id, 1);
      vK = zeros(id, 1);
      vL = zeros(id, 1);
      vu = zeros(id, 1);

      for t = 1 : id
        ve(t) = vha(t) - da;           % Kalman filter 
        vDinv(t) = 1 / (dP + vsiga2(t));
        vK(t) = dP * vDinv(t);
        vL(t) = 1 - vK(t);
        
        da = da + vK(t) * ve(t);
        dP = dP * vL(t) + dH2;
      end

      if j <= nite
                    % finding mode
        dr = 0;       
        dU = 0;
        
        t = id;                        % simulation smoother
        while t >= 1
          dC = dH2 * (1 - dU * dH2);
          deps = 0;
          vu(t) = dH2 * dr + deps;
          dV = dH2 * dU * vL(t);
          
          dCinv = 1 / dC;
          dr = vDinv(t) * ve(t) + vL(t) * dr ...
             - dV * dCinv * deps;
          dU = vDinv(t) + vL(t)^2 * dU + dV^2 * dCinv;
          t = t - 1;
        end
        
        du0 = dH20 * dr;
        vhn(1) = dh0 + du0;
        for t = 1 : id-1
          vhn(t+1) = vhn(t) + vu(t);
        end
      
      else
                    % draw candidate
        fl = 0;
        icyc = 0;
        while (fl == 0) && (icyc < 100)
            
          dr = 0;
          dU = 0;
        
          t = id;                      % simulation smoother
          while t >= 1
            dC = dH2 * (1 - dU * dH2);   
            deps = sqrt(dC) * randn();
            vu(t) = dH2 * dr + deps;
            dV = dH2 * dU * vL(t);
          
            dCinv = 1 / dC;
            dr = vDinv(t) * ve(t) + vL(t) * dr ...
               - dV * dCinv * deps;
            dU = vDinv(t) + vL(t)^2 * dU + dV^2 * dCinv;
            t = t - 1;
          end

          dC = dH20 * (1 - dU * dH20);
          du0 = dH20 * dr + sqrt(dC) * randn();
          vhn(1) = dh0 + du0;
          for t = 1 : id-1
            vhn(t+1) = vhn(t) + vu(t);
          end

         %%--- AR step ---%%
         
          dpron ...
           = sum( ...
             -0.5 * (vhh + vyi.^2 ./ exp(vhh)) ...% g(h_hat)
             + vgder1 .* (vhn - vhh) ...         % g'(h_hat)
             + 0.5 * vgder2 .* (vhn - vhh).^2); % g''(h_hat)
          
          dposn = sum(-0.5 * (vhn + vyi.^2 ./ exp(vhn)));
          
          if m_k >= 0
              m_iachar = m_iachar + 1;
          end
          
          if rand() < exp(dposn - dpron)
              fl = 1;
          end
          
          icyc = icyc + 1;

        end
                    
      end

  end
  
  if icyc < 100

     %%--- MH step ---%%
     
      dproo ...
       = sum( ...
          -0.5 * (vhh + vyi.^2 ./ exp(vhh)) ... % g(h_hat)
          + vgder1 .* (vho - vhh) ...          % g'(h_hat)
          + 0.5 * vgder2 .* (vho - vhh).^2);  % g''(h_hat)

      dposo = sum(-0.5 * (vho + vyi.^2 ./ exp(vho)));
      
      dfrac = exp(  dposn + min([dposo; dproo]) ...
                  - dposo - min([dposn; dpron]) );

      if rand() < dfrac
          vhs(ir:ird) = vhn;
          if m_k >= 0
              m_iachmh = m_iachmh + 1;
          end
      end
  end
  
end  