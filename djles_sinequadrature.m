function w = djles_sinequadrature(NX,NZ,H,L)
% Gets the weights for 2D quadrature over the domain
wtx=djles_quadweights(NX)*(L/pi);
wtz=djles_quadweights(NZ)*(H/pi);
w=bsxfun(@times,wtx,wtz.');
end
