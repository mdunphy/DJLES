clear
%%% Specify the parameters of the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Approximate reproduction of the smoothed curve from Figure 9a from
% Pineda et al (2015). Courtesy Jorge Magalhães.
data = load('case_pineda_cast1.txt'); zdata = data(:,1); rhodata = data(:,2);

rho0 = max(rhodata);  % Used in djles_common to compute N2(z)

% Use MATLAB's gradient function to find the derivative drho/dz
rhozdata = gradient(rhodata,zdata);

% Now build piecewise interpolating polynomials from the data,
% and convert them to function handles for the solver
method='pchip';  % pchip respects monotonicity
warning('off','MATLAB:interp1:ppOutput')
rho  = @(z) ppval(interp1(zdata,rhodata ,method,'pp'), z);
rhoz = @(z) ppval(interp1(zdata,rhozdata,method,'pp'), z);

% The velocity profile (zero for this case) (m/s)
Ubg=@(z) 0*z; Ubgz=@(z) 0*z; Ubgzz=@(z) 0*z;

A  =  5e3;  % APE for wave (kg m/s^2)
L  = 1200;  % domain width (m)
H  =   57;  % domain depth (m), estimated from Pineda et al's Figure 9a
NX =   64;  % grid
NZ =   32;  % grid

epsilon=1e-3;   % use a larger epsilon for intermediate waves

%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = clock;

djles_refine_solution
djles_diagnostics; djles_plot

% Raise amplitude in a few increments
for A=linspace(1e4, 3.62e5, 6)
    % Find the solution of the DJL equation
    djles_refine_solution
    djles_diagnostics; djles_plot
end

% Increase resolution, reduce epsilon, iterate to convergence
NX=64; NZ=64; epsilon=1e-5;
djles_refine_solution
djles_diagnostics; djles_plot

NX=128; NZ=128; epsilon=1e-6;
djles_refine_solution
djles_diagnostics; djles_plot

NX=256; NZ=256; epsilon=1e-7;
djles_refine_solution
djles_diagnostics; djles_plot

NX=512; NZ=512; epsilon=1e-8;
djles_refine_solution
djles_diagnostics; djles_plot

end_time=clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));
