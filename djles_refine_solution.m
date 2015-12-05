% djles_refine_solution.m
% - Iterates eta to convergence following the procedure
%   described in Stastna and Lamb, 2002 (SL2002) and
%   also in Dunphy, Subich and Stastna 2011 (DSS2011).

djles_common

% Initialise t for recording timing data
t.start = clock; t.solve = 0; t.int = 0;

% Get initial guess from WNL theory if needed
if ~exist('eta','var'), djles_initial_guess; end

% Check for nonzero velocity profile, save time below if it's zero
uflag = any(Ubg(zc));

% Change the resolution if needed
eta = djles_change_resolution(eta,NX,NZ);

if (verbose >= 1)
    [~,idx]   = max(abs(eta(:)));
    wave_ampl = eta(idx);
    fprintf('Initial guess:\n wave ampl = %+.10e,   c = %+.10e\n\n',wave_ampl,c);
end

flag = 1; iteration = 0;
while (flag)
    % Iteration shift
    iteration = iteration + 1;
    eta0      = eta;
    c0        = c;
    lambda0   = g*H/(c*c);

    S = N2(ZC-eta0).*eta0/(g*H); % compute S (DSS2011 Eq 19)

    % Compute R, assemble RHS
    if (uflag)
        [eta0x, eta0z] = djles_gradient(eta0, ks, ms, 'odd', 'odd', 'interior');
        uhat  = Ubg(ZC-eta0)/c0;
        uhatz = Ubgz(ZC-eta0)/c0;
        R = (uhatz ./ (uhat -1)) .* ( 1 - (eta0x.^2 + (1-eta0z).^2));
        rhs = - ( lambda0 * (S./((uhat-1).^2))  + R );
    else
        rhs = - lambda0 * S;  % RHS of (DSS2011 Eq 18)
    end

    % Solve the linear Poisson problem (DSS2011 Eq 18)
    t0 = clock;
    temp = djles_extend(rhs, 'odd', 'odd', 'interior');
    temp = real(ifft2(INVLAP.*fft2(temp)));
    nu = temp(1:NZ, 1:NX);
    t.solve = t.solve + etime(clock, t0);

    % Find the APE (DSS2011 Eq 21 & 22)
    t0 = clock;
    apedens = djles_compute_apedens(rho, eta0, ZC, g, wl, zl);
    F = sum(wsine(:).*apedens(:));
    % Compute S1, S2 (components of DSS2011 Eq 20)
    S1 = g*H*rho0*sum( wsine(:).*S(:).*nu(:)   );
    S2 = g*H*rho0*sum( wsine(:).*S(:).*eta0(:) );
    t.int = t.int + etime(clock, t0);

    % Find new lambda (DSS2011 Eq 20)
    lambda = lambda0*(A-F+S2)/S1;

    % check if lambda is OK
    if (lambda < 0)
        fprintf('new lambda has wrong sign --> nonconvergence of iterative procedure\n');
        fprintf('new lambda = %1.6e\n',lambda);
        fprintf('   A = %1.10e, F = %1.10e\n',A,F);
        fprintf('   S1 = %1.8e, S2 = %1.8e, S2/S1 = %1.8e\n',S1,S2,S2/S1);
        break;
    end

    % Compute new c, eta
    c   = sqrt(g*H/lambda); % DSS2011 Eq 24
    eta = (lambda/lambda0)*nu; % (DSS2011 Eq 23)

    % Apply underrelaxation factor
    eta = (1-relax)*eta0 + relax*eta;

    % Find wave amplitude
    [~,idx]   = max(abs(eta(:)));
    wave_ampl = eta(idx);

    % Compute relative difference between present and previous iteration
    reldiff = max(abs(eta(:)-eta0(:))) / abs(wave_ampl);

    % Report on state of the operation
    if (verbose >=1)
        fprintf('Iteration %4d:\n',iteration);
        fprintf(' A       = %+.10e, wave ampl = %+16.10f m\n',A,wave_ampl);
        fprintf(' F       = %+.10e, c         = %+16.10f m/s\n',F,c);
        fprintf(' reldiff = %+.10e\n\n',reldiff);
    end

    % Stop conditions
    if (iteration >= min_iteration) && (reldiff < epsilon)
        flag = 0;
    end
    if (iteration >= max_iteration)
        flag = 0;
        fprintf('Reached maximum number of iterations (%d >= %d)\n',iteration,max_iteration);
    end
end

t.stop = clock; t.total = etime(t.stop, t.start);
% Report the timing data
if (verbose >= 2)
    fprintf('Poisson solve time: %6.2f seconds\n', t.solve);
    fprintf('Integration time:   %6.2f seconds\n', t.int);
    fprintf('Other time:         %6.2f seconds\n', t.total - t.solve - t.int);
    fprintf('Total:              %6.2f seconds\n', t.total);
end

fprintf('Finished [NX,NZ]=[%3dx%3d], A=%g, c=%g m/s, wave amplitude=%g m\n',NX,NZ,A,c,wave_ampl);

% Cleanup unneeded variables (comment these lines for debugging)
clear uflag S uhat R rhs temp nu apedens F S1 S2 flag idx uhatz eta0x eta0z t0 eta0 lambda0 c0
