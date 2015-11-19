function w = djles_sinequadrature(NX,NZ,H,L)
% Gets the weights for 2D quadrature over the domain
wtx=get1dweights(NX)*(L/pi);
wtz=get1dweights(NZ)*(H/pi);
w=bsxfun(@times,wtx,wtz.');
end

function w=get1dweights(N)
% Quadrature weights for integrating an odd periodic function over a
% half period using the interior grid. From John Boyd's Chebyshev and
% Fourier Spectral Methods, 2nd Ed, Pg 568, Eq F.32.
w=zeros(1,N);
for j=1:N
    xj=(2*j-1)*pi/(2*N);
    m=1:N-1;
    s=sum(sin(m*xj).*(sin(m*pi/2).^2)./m);
    w(j) = (2/N/N)*sin(N*xj)*(sin(N*pi/2).^2) + (4/N)*s;
end
end
