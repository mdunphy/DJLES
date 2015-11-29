function [residual, LHS, RHS] = djles_residual(ks, ms, eta, Ubg, Ubgz, N2, z, c, gridtype)
% DJL residual using Eq 2.32 in (Stastna, 2001)

% Odd extend the function in both directions
etaextended = djles_extend(eta, 'odd', 'odd', gridtype);

[SZ,SX] = size(eta);

% Compute left hand side
LAP = -bsxfun(@plus,ms.^2,ks.^2);           % Laplacian operator
LHS = real(ifft2(LAP.*fft2(etaextended)));  % Laplacian of extended eta
LHS = LHS(1:SZ, 1:SX);                      % Trim for eta on gridtype

% Compute right hand side
[etax, etaz] = djles_gradient(eta, ks, ms, 'odd', 'odd', gridtype);
Umc = Ubg(z-eta)-c;
RHS = -(Ubgz(z-eta)./Umc).*(1 - (etax.^2 + (1-etaz).^2)) - N2(z-eta).*eta./(Umc.^2);

% Residual
residual = LHS-RHS;
end