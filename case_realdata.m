clear all

%%% Specify the parameters of the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load sample data file: two columns, depth and temperature
data = load('case_realdata.txt'); zdata = data(:,1); T = data(:,2);

% Find a value for z=0 by linear extrapolation
T = [T(:); interp1(zdata,T,0,'linear','extrap')];
zdata = [zdata(:); 0];

% Convert to density with linear EOS, normalized by rho0
rho0=1000; T0=10; alpha=1.7e-4;
rhodata = rho0*(1 -alpha*(T-T0))/rho0;

A  = 5;    % APE for wave (m^4/s^2)
L  = 500;  % domain width (m)
H  = 16.5; % domain depth (m)
NX = 64;   % grid
NZ = 32;   % grid

relax=0.1; % use strong underrelaxation
verbose=1;

% Use MATLAB's gradient function to find the derivative drho/dz
rhozdata = gradient(rhodata,zdata);

% Now build piecewise interpolating polynomials from the data,
% and convert them to function handles for the solver
%method='linear'; % linear works well
method='pchip';  % pchip also respects monotonicity
rho  = @(z) ppval(interp1(zdata,rhodata ,method,'pp'), z);
rhoz = @(z) ppval(interp1(zdata,rhozdata,method,'pp'), z);

% The velocity profile (zero for this case) (m/s)
Ubg=@(z) 0*z; Ubgz=@(z) 0*z; Ubgzz=@(z) 0*z;

%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = clock;
djles_refine_solution;

% Increase resolution, iterate to convergence
NX=128; NZ=128;
djles_refine_solution

% Increase to the final resolution, iterate to convergence
NX=512; NZ=512;
djles_refine_solution

end_time=clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

% Compute and plot the diagnostics
djles_diagnostics
djles_plot
