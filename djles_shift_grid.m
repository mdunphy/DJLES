function fe = djles_shift_grid(fc, NX, NZ, symmx, symmz)
% Shifts the data fc from the interior grid to the endpoint grid
% Specify symmetry as [] to skip shifting a dimension

[SZ,SX] = size(fc); % Size of incoming grid

if ~isequal(SX, NX+1) && ~isempty(symmx)
    % Shift grid in X using FFT interpolation
    eta0in = djles_extend(fc, symmx, [], 'interior');
    eta0out = real(interpft(eta0in, 4*NX, 2));
    eta0out = circshift(eta0out, [0 1]);
    fc = eta0out(:, 1:2:2*NX+1);
end

if ~isequal(SZ, NZ+1) && ~isempty(symmz)
    % Shift grid in Z using FFT interpolation
    eta0in = djles_extend(fc, [], symmz, 'interior');
    eta0out = real(interpft(eta0in, 4*NZ, 1));
    eta0out = circshift(eta0out, [1 0]);
    fc = eta0out(1:2:2*NZ+1, :);
end

fe=fc;
end
