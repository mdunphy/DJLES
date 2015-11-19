function f = compat_flip(f,dim)
% Wrapper function for MATLAB's flipdim -> flip replacement
try
    f=flip(f,dim);    % R2013a onward
catch
    f=flipdim(f,dim); % Prior to R2013a
end

