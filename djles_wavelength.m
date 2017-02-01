function [wavelength] = djles_wavelength(eta,L)
% L_w following via Eq 3.6 in
% Aghsaee, P., Boegman, L., and K. G. Lamb. 2010. "Breaking of shoaling internal
% solitary waves". J. Fluid Mech. 659 289-317. doi:10.1017/S002211201000248X.
% We return 2*L_w as the wavelength
[maxeta, idx] = max(abs(eta(:)));
[iz, ~] = ind2sub(size(eta),idx);
etaL = abs(eta(iz,:));
[~,NX] = size(eta);
w = djles_quadweights(NX)*(L/pi);
Lw = sum(w.*etaL) / maxeta;
% We take wavelength as twice Lw
wavelength = 2*Lw;
end
