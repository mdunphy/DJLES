function fe = djles_shift_grid(fc, NX, NZ, symmetry)
% Shifts the data fc from the interior grid to the endpoint grid

switch symmetry
    case 'odd'
        s=-1;
    case 'even'
        s=1;
    otherwise
        warning('You must specify even or odd symmetry here')
end

[SZ,SX] = size(fc); % Size of incoming grid

if isequal([SZ SX], [NZ NX]+1)
    % We're already on the endpoint grid (ie, this got called twice)
    fe=fc;
else
    % We're on the interior grid, need to move to endpoint grid
    % Shift grid in X using FFT interpolation
    eta0in = [fc s*flipdim(fc,2)];
    eta0out = real(interpft(eta0in, 4*NX, 2));
    eta0out = circshift(eta0out, [0 1]);
    temp = eta0out(:, 1:2:2*NX+1);
    
    % Shift grid in Z using FFT interpolation
    eta0in = [temp; s*flipdim(temp,1)];
    eta0out = real(interpft(eta0in, 4*NZ, 1));
    eta0out = circshift(eta0out, [1 0]);
    fe = eta0out(1:2:2*NZ+1, :);
end
end
