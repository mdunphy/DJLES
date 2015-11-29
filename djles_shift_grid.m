function fe = djles_shift_grid(fc, NX, NZ, symmx, symmz)
% Shifts the data fc from the interior grid to the endpoint grid

[SZ,SX] = size(fc); % Size of incoming grid

if isequal([SZ SX], [NZ NX]+1)
    % We're already on the endpoint grid (ie, this got called twice)
    fe=fc;
else
    % We're on the interior grid, need to move to endpoint grid
    % Shift grid in X using FFT interpolation
    eta0in = djles_extend(fc, symmx, [], 'interior');
    eta0out = real(interpft(eta0in, 4*NX, 2));
    eta0out = circshift(eta0out, [0 1]);
    temp = eta0out(:, 1:2:2*NX+1);
    
    % Shift grid in Z using FFT interpolation
    eta0in = djles_extend(temp, [], symmz, 'interior');
    eta0out = real(interpft(eta0in, 4*NZ, 1));
    eta0out = circshift(eta0out, [1 0]);
    fe = eta0out(1:2:2*NZ+1, :);
end
end
