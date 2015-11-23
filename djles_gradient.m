function [fx, fz] = djles_gradient(f, ks, ms, symmx, symmz, gridtype)
% Computes gradient using the specified symmetry and grid type

sx=getsymm(symmx);
sz=getsymm(symmz);

switch gridtype
    case 'interior'
        % Cell centred or 'interior' grid
        [SZ,SX] = size(f);
        fx = real(ifft(bsxfun(@times,1i*ks(:).',fft(cat(2,f,sx*compat_flip(f,2)),[],2)),[],2));
        fx = fx(1:SZ,1:SX);
        fz = real(ifft(bsxfun(@times,1i*ms(:)  ,fft(cat(1,f,sz*compat_flip(f,1)),[],1)),[],1));
        fz = fz(1:SZ,1:SX);
    case 'endpoint'
        % Cell edges grid or 'endpoint' grid, which is one larger in
        % each dimension than the interior grid
        [SZ,SX] = size(f);
        fx = real(ifft(bsxfun(@times,1i*ks(:).',fft(cat(2,f,sx*compat_flip(f(:,2:end-1),2)),[],2)),[],2));
        fx = fx(1:SZ,1:SX);
        fz = real(ifft(bsxfun(@times,1i*ms(:)  ,fft(cat(1,f,sz*compat_flip(f(2:end-1,:),1)),[],1)),[],1));
        fz = fz(1:SZ,1:SX);
    otherwise
        warning('You must specify interior or endpoint grid here')
end
end

function s = getsymm(symm)
switch symm
    case 'odd'
        s=-1;
    case 'even'
        s=1;
    otherwise
        warning('You must specify even or odd symmetry here')
end
end
