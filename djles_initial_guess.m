% djles_initial_guess.m
% - Finds an initial guess for eta and c via weakly nonlinear theory

djles_common

Dz   = djles_diffmatrix(dz, NZ, 1, 'not periodic');
Dzz  = djles_diffmatrix(dz, NZ, 2, 'not periodic');
Dzzc = Dzz(2:end-1, 2:end-1); % cut the endpoints for Dirichlet conditions

% get n2, u and uzz data
tmpz   = zc(2:end-1)';
n2vec  = N2(tmpz);
uvec   = Ubg(tmpz);
uzzvec = Ubgzz(tmpz);

% create diagonal matrices
N2d  = diag(n2vec );
Ud   = diag(uvec  );
Uzzd = diag(uzzvec);

% setup quadratic eigenvalue problem
B0 = sparse(N2d + Ud.*Ud*Dzzc - Ud*Uzzd);
B1 = sparse(-2*Ud*Dzzc + Uzzd);
B2 = sparse(Dzzc);

% Solve eigenvalue problem; extract first eigenvalue & eigenmode
if ~exist('OCTAVE_VERSION', 'builtin')
    % We're in MATLAB, proceed normally
    [V,cc]=polyeig(B0,B1,B2);
else
    % We're in Octave
    if exist('polyeig','file')
        % This version of Octave has polyeig, use it
        [V,cc]=polyeig(B0,B1,B2);
    else
        error(mfilename, 'polyeig function not found! see comments here for workaround');
        % This version of Octave is ancient and does not have polyeig.m
        % You can grab polyeig.m from Octave's source code repository and
        % save it in the DJLES directory to proceed. This version works:
        %  http://hg.savannah.gnu.org/hgweb/octave/file/ca1648b2e673/scripts/polynomial/polyeig.m
    end
    % Convert cc from diagonal_matrix to vector if needed
    if ~isvector(cc), cc=diag(cc); end
end

% Sort eigenvalues, largest is clw
[c,ii]=sort(cc,'descend'); clw=c(1);

% Add boundary conditions
phi1=[0; V(:,ii(1)); 0];
uvec=[Ubg(-H); uvec; Ubg(0)];

% Compute E1, normalise
E1=clw*phi1./(clw-uvec);
E1=abs(E1)/max(abs(E1));

% Compute r10 and r01
E1p=Dz*E1; E1p2=E1p.^2; E1p3=E1p.^3;
bot=sum((clw-uvec).*E1p2);
r10=(-0.75/clw)*sum((clw-uvec).*(clw-uvec).*E1p3)/bot;
r01=-0.5*sum((clw-uvec).*(clw-uvec).*E1.*E1)/bot;
if (verbose)
    fprintf('WNL gives: c_lw = %f, r10 = %f, r01 = %f\n\n', clw, r10, r01);
end

% Now optimise the b0, lambda parameters
E = repmat(E1, 1, NX); E(:,1)=0; E(:,end)=0;
b0=sign(r10)*0.05*H;  % Start b0 as 5% of domain height
lambda = sqrt( -6*r01 / (clw * r10 * b0) );
c = (1+(2/3)*r10*b0)*clw;
if (verbose)
    fprintf('init b0 = %f, lambda = %f, V = %f\n', b0, lambda, c);
end

flag = 1;
while (flag>0)
    flag=flag+1;
    eta0 = -b0*E.*sech((XC-L/2)/lambda).^2;
    eta0(1,:)=0; eta0(:,1)=0; eta0(:,end)=0; eta0(end,:)=0;
    
    % Find the APE (DSS2011 Eq 21 & 22)
    apedens = djles_compute_apedens(rho, eta0, ZC, g, wl, zl);
    F = sum(wsine(:).*apedens(:));
    
    % New b0 by rescaling
    afact = max(min(A/F, 1.05), 0.95);
    b0 = b0 * afact;
    
    % New c, lambda
    c = (1+(2/3)*r10*b0)*clw;
    lambda = sqrt( -6*r01 / (clw * r10 * b0) );
    if ~isreal(lambda)
        disp('problem finding new lambda !!')
    end
    
    if (verbose >= 1)
        fprintf('F=%e, desired = %e, rescaling b0 by factor of %f...\n',F,A,afact);
        fprintf('new b0 = %f, lambda = %f, V = %f\n\n',b0,lambda,c);
    end
    
    % Stop conditions: the wave gets too big, or we get matching APE
    if abs(b0) > 0.75*H || abs(afact-1) < 0.01
        eta=eta0;
        flag=0;
    end
end

% Cleanup unneeded variables (comment these lines for debugging)
clear Dz Dzz Dzzc tmpz n2vec uvec uzzvec N2d Ud Uzzd B0 B1 B2 V cc ii flag F clw
clear phi1 uvec E1 E1p E1p2 E1p3 bot r10 r01 E b0 lambda apedens afact eta0
