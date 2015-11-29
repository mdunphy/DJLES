function [fx, fz] = djles_gradient(f, ks, ms, symmx, symmz, gridtype)
% Computes gradient using the specified symmetry and grid type
[SZ,SX] = size(f);
fx = real(ifft(bsxfun(@times,1i*ks(:).',fft(djles_extend(f, symmx, [], gridtype),[],2)),[],2));
fx = fx(1:SZ,1:SX);
fz = real(ifft(bsxfun(@times,1i*ms(:)  ,fft(djles_extend(f, [], symmz, gridtype),[],1)),[],1));
fz = fz(1:SZ,1:SX);
