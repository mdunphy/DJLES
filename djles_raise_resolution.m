function eta = djles_raise_resolution(eta0,NX,NZ)
% Increases the resolution of eta to NZxNX
% This function is idempotent, that is, it does nothing if the
% requested size matches the old size.
% Size should be an integer multiple of old size (2x, 4x, etc).

[NZ0,NX0] = size(eta0);  % old size

% use FFT interpolation
if ~isequal(NX0,NX)
    RX=NX/NX0;
    
    % Change resolution in X
    eta0in = [eta0 -fliplr(eta0)];
    eta0out = real(interpft(eta0in, 4*NX, 2));
    eta0out = circshift(eta0out, [0 (RX-1)]);
    eta0 = eta0out(:, 1:2:2*NX);
end

if ~isequal(NZ0,NZ)
    RZ=NZ/NZ0;
    
    % Change resolution in Z
    eta0in = [eta0; -flipud(eta0)];
    eta0out = real(interpft(eta0in, 4*NZ, 1));
    eta0out = circshift(eta0out, [(RZ-1) 0]);
    eta0 = eta0out(1:2:2*NZ, :);
end

eta = eta0;
end
