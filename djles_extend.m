function fe = djles_extend(f, symmx, symmz, gridtype)
% Extends the function f using symmetry symmx, symmz in x and z
% Specify symmetry as [] to skip extending a dimension

fe = extendx(f , symmx, gridtype);
fe = extendz(fe, symmz, gridtype);
end

function fe = extendx(f, symm, gridtype)
% Extend the function in the x dimension using symmetry symm
if ~isempty(symm)
    s=getsymm(symm);
    switch gridtype
        case 'interior'
            fe = cat(2, f, s*compat_flip(f,2));
        case 'endpoint'
            fe = cat(2, f, s*compat_flip(f(:,2:end-1),2));
    end
else
    fe = f;
end
end

function fe = extendz(f, symm, gridtype)
% Extend the function in the z dimension using symmetry symm
if ~isempty(symm)
    s=getsymm(symm);
    switch gridtype
        case 'interior'
            fe = cat(1, f, s*compat_flip(f,1));
        case 'endpoint'
            fe = cat(1, f, s*compat_flip(f(2:end-1,:),1));
    end
else
    fe = f;
end
end

function s = getsymm(symm)
% Converts the symmetry to a sign
switch symm
    case 'odd'
        s=-1;
    case 'even'
        s=1;
    otherwise
        s=0;
end
end

function f = compat_flip(f,dim)
% Wrapper function for MATLAB's flipdim -> flip replacement
try
    f=flip(f,dim);    % R2013a onward
catch
    f=flipdim(f,dim); % Prior to R2013a
end
end
