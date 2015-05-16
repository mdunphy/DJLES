% Constructs a centered differentiation matrix with up/down wind at the ends
% - order can be 1 (1st deriv) or 2 (2nd deriv)
function D = djles_diffmatrix(dx, N, order, ends)
a=[0;0;0]; a(order+1)=1; dx2 = dx*dx;

% Upwind 2nd order coefs
matu=[ 1      1     1  ; ...
      -2*dx  -dx    0  ; ...
       2*dx2  dx2/2 0 ];
coefsu=matu \ a;

% Centered 2nd order coefs
matc=[ 1      1  1      ; ...
      -dx     0  dx     ; ...
       dx2/2  0  dx2/2 ];
coefsc=matc \ a;

% Downwind 2nd order coefs
matd=[ 1  1     1     ; ...
       0  dx    2*dx  ; ...
       0  dx2/2 2*dx2 ];
coefsd=matd \ a;

D = zeros(N,N);

% Centered in the interior
for ii=2:N-1,
    D (ii,ii+[-1 0 1])=coefsc;
end

switch ends
    case 'periodic'
        % use centered at the ends (wrap around)
        D(1, [N   1  2]) = coefsc;
        D(N, [N-1 N  1]) = coefsc;
    otherwise
        % Downwind at the first point
        D(1, 1+[0 1 2]) = coefsd;
        % Upwind at the last point
        D(N, N+[-2 -1 0]) = coefsu;
end
end
