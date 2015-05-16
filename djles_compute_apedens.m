function apedens = djles_compute_apedens(rho, eta, z, g, wl, zl)
% Use Gauss quadrature to find ape density
apedens = rho(z-eta);
for ii=1:numel(wl)
    apedens = apedens - wl(ii).*rho(z-zl(ii)*eta);
end
apedens = g*apedens.*eta;
