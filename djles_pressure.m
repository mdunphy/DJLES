% djles_pressure.m - Pressure computations. Specifically, we compute
%  a) the hydrostatic background pressure
%  b) the hydrostatic pressure due to the internal solitary wave
%  c) the non-hydrostatic pressure due to the internal wave
%  d) residuals from substituting the results into the governing equations
%
% The pressure solve part is courtesy of Derek Steinmoeller.
%
% Note: This is only tested/verified for cases without a background current.
% ToDo: Make it work with background current cases.


%%% Hydrostatic background pressure %%%

% We break rho(z) into two parts: a constant-N background state and an odd
% function of z:  rho(z) = rholin(z) + rhoodd(z)
rhomax = rho(-H); rhomin=rho(0); drho = rhomax-rhomin;
rholin = @(z) rhomin + (z/-H)*drho;
rhoodd = @(z) rho(z) - rholin(z);

% We integrate dplin/dz = -g*rholin(z) analytically to get plin(z)
plin = @(z) -g*rhomin*z + g*z.^2/(2*H)*drho;

% Now we use a sine transform to integrate dpodd/dz = -g*rhoodd(z)
temp = djles_extend(-g*rhoodd(Z(:,1)), [], 'odd', targetgrid); % extend
temp = real(ifft(1i*msi.*fft(temp)));                          % integrate
podd = temp(1:length(Z(:,1)));                                 % extract
poddze = djles_shift_grid(podd, NX, NZ, [], 'even'); % get podd on z endpoint grid
podd = podd - poddze(end);     % enforce podd(0)=0

% Total background pressure is plin+podd
pbg = plin(Z(:,1)) + podd;


%%% Wave pressure %%%

% First find the internal solitary wave hydrostatic pressure
rhowave = rho(Z-eta) - rho(Z);
% Integrate dph/dz = -g*rhowave vertically (odd, odd symmetry in x,z)
temp = djles_extend(-g*rhowave, 'odd', 'odd', targetgrid);    % extend
temp = real(ifft(bsxfun(@times,1i*msi,fft(temp,[],1)),[],1)); % integrate
ph = temp(1:size(eta,1), 1:size(eta,2));                      % extract
phze = djles_shift_grid(ph, NX, NZ, [], 'even'); % get ph on z endpoint grid
ph = bsxfun(@minus, ph, phze(end,:));  % Enforce ph(x,0) = 0

% Now solve for total wave pressure
T1 = -((u-c).*ux + w.*uz);                    % odd, even symmetry in x,z
T2 = -((u-c).*wx + w.*wz + (g/rho0)*rhowave); % even, odd symmetry in x,z
% Get divergence of T vector
[Tx, ~] = djles_gradient(T1, ks, ms, 'odd', 'even', targetgrid);
[~, Tz] = djles_gradient(T2, ks, ms, 'even', 'odd', targetgrid);
% The Poisson problem has even, even symmetry in x,z
temp = djles_extend(rho0*(Tx+Tz), 'even', 'even', targetgrid);   % extend
temp = real(ifft2(INVLAP.*fft2(temp)));                          % invert
p = temp(1:size(eta,1), 1:size(eta,2));                          % extract
p = p - p(end,end);  % Use p(inf, 0) = 0 (no wave pressure in the farfield)
pnh = p - ph;  % non-hydrostatic part

clear T1 T2 Tx Tz phze poddze temp


%%% Residuals %%%
% Now we verify that the results are a steady state solution to
% the Boussinesq equations (trend terms and divergence should be zero)
[px ,pz  ] = djles_gradient(p      , ks, ms, 'even', 'even', targetgrid);
[~  ,pnhz] = djles_gradient(pnh    , ks, ms, 'even', 'even', targetgrid);
[~  ,phz ] = djles_gradient(ph     , ks, ms, 'even', 'even', targetgrid);
[rwx,rwz ] = djles_gradient(rhowave, ks, ms, 'odd' , 'odd' , targetgrid);

% Boussinesq equations
ut = -px/rho0 -(u-c).*ux -w.*uz;
wt = -pz/rho0 -(u-c).*wx -w.*wz - (g/rho0)*rhowave;
div = ux+wz;
rt = -(u-c).*rwx - w.*(rhoz(Z) + rwz);

% Also verify that the split pressure satisfies the vertical
% momentum eqn and hydrostatic balance
wtnh = -pnhz/rho0 -(u-c).*wx -w.*wz;
hb = -phz/rho0 -(g/rho0)*rhowave;

% Compute relative errors for each
ue=max(abs(ut (:))) / max(abs(-px(:)/rho0));
we=max(abs(wt (:))) / max(abs(-pz(:)/rho0));
de=max(abs(div(:))) / max(abs(ux(:)));
re=max(abs(rt (:))) / max(abs(-(u(:)-c).*rwx(:)));

wnhe=max(abs(wtnh(:))) / max(abs(-pnhz(:)/rho0));
hbe =max(abs(hb  (:))) / max(abs(phz(:)));

fprintf('ut, wt, div, rho_t residuals: %.1e, %.1e, %.1e, %.1e\n',ue,we,de,re);
fprintf('wnht, hb, residuals: %.1e, %.1e\n',wnhe,hbe);
