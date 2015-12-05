function eta = djles_change_resolution(eta0,NX,NZ)
% Changes the resolution of eta to NZxNX
% This function is idempotent, that is, it does nothing if the
% requested size matches the old size.

[NZ0,NX0] = size(eta0);
if ~isequal(NX0,NX)
    % Change resolution in X dimension
    
    RX=NX/NX0;
    if RX>1 && double(int32(RX))==RX
        % Raise resolution in X using interpft
        eta0in = djles_extend(eta0, 'odd', [], 'interior');
        eta0out = real(interpft(eta0in, 4*NX, 2));
        eta0out = circshift(eta0out, [0 (RX-1)]);
        eta0 = eta0out(:, 1:2:2*NX);
    else
        % Change resolution in X using linear interpolation
        x0 = (  0:NX0)*(1/NX0);  % old x endpoints
        x  = (0.5:NX )*(1/NX );  % new x centres
        tmp = djles_shift_grid(eta0, NX0, NZ0, 'odd', []); % endpoint grid in x
        eta0 = interp2(x0,(1:NZ0).',tmp,x,(1:NZ0).');
    end
end

[NZ0,NX0] = size(eta0);
if ~isequal(NZ0,NZ)
    % Change resolution in Z dimension
    
    RZ=NZ/NZ0;
    if RZ>1 && double(int32(RZ))==RZ
        % Raise resolution in Z using interpft
        eta0in = djles_extend(eta0, [], 'odd', 'interior');
        eta0out = real(interpft(eta0in, 4*NZ, 1));
        eta0out = circshift(eta0out, [(RZ-1) 0]);
        eta0 = eta0out(1:2:2*NZ, :);
    else
        % Change resolution in Z using linear interpolation
        z0 = (  0:NZ0)*(1/NZ0);  % old z endpoints
        z  = (0.5:NZ )*(1/NZ );  % new z centres
        tmp = djles_shift_grid(eta0, NX0, NZ0, [], 'odd');  % endpoint grid in z
        eta0 = interp2(1:NX0,z0.',tmp,1:NX0,z.');
    end
end

eta = eta0;
end
