% Computes a variety of diagnostics from the solved eta and c

djles_common

%%% Target grid selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
targetgrid='interior';
% targetgrid='endpoint';

if isequal(targetgrid, 'endpoint')
    % DJLES computes eta on the interior grid; if we want diagnostics on
    % the endpoint grid then we need to do a shift.
    eta=djles_shift_grid(eta,NX,NZ,'odd','odd');
    Z=ZE; z=ze; x=xe;
else
    Z=ZC; z=zc; x=xc;
end

%%% Compute the diagnostics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute velocities (Via SL2002 Eq 27)
[etax, etaz] = djles_gradient(eta, ks, ms, 'odd', 'odd', targetgrid);
u = Ubg(Z-eta).*(1-etaz) + c*etaz;
w = -Ubg(Z-eta).*(-etax) - c*etax;
uwave = u - Ubg(Z);

% Kinetic energy perturbation (m^2/s^2) - first and second order terms
kepert = Ubg(Z).*uwave + 0.5*(uwave.^2 + w.^2);

% Wave kinteic energy density (m^2/s^2) - second order terms only
kewave = 0.5*(uwave.^2 + w.^2);
%[kewavex, kewavez] = djles_gradient(kewave, ks, ms, 'even', targetgrid);

% APE density (m^2/s^2)
apedens = djles_compute_apedens(rho, eta, Z, g, wl, zl);

% Get gradient of u and w
[ux, uz] = djles_gradient(u, ks, ms, 'even', 'even', targetgrid);
[wx, wz] = djles_gradient(w, ks, ms, 'even', 'odd', targetgrid);

% Surface strain rate
uxze=djles_shift_grid(ux,NX,NZ,[],'even'); % shift ux to z endpoints
surf_strain = -uxze(end,:); % = -du/dx(z=0)

% Vorticity, density and Richardson number
vorticity = uz - wx;
density = rho(Z-eta);
ri = N2(Z-eta)./(uz.*uz);

% Wavelength (currently works only on interior grid)
wavelength = djles_wavelength(eta,L);

% Residual in DJL equation
[residual, LHS, RHS] = djles_residual(ks, ms, eta, Ubg, Ubgz, N2, Z, c, targetgrid);
fprintf('Relative residual %e\n',max(abs(residual(:))) / max(abs(LHS(:))));
