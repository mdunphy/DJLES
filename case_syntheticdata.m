clear all

%%% Specify the parameters of the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A  = 1e-4; % APE for wave (m^4/s^2)
L  = 8.0;  % domain width (m)
H  = 0.2;  % domain depth (m)
NX = 32;   % grid
NZ = 32;   % grid

% Specify analytic density and velocity profiles
a_d=0.02; z0_d=0.05; d_d=0.01;
frho=@(z) 1-a_d*tanh((z+z0_d)/d_d);

U0=0.10; zj=0.5*H; dj=0.4*H;
fUbg=@(z) U0*tanh((z+zj)/dj);

% Now sample the density & velocity profiles to get synthetic mooring data
NDATA=26;
zdata=linspace(-H, 0, NDATA).';
rhodata=frho(zdata);
ubgdata=fUbg(zdata);

% Use numerical derivatives for the gradients of the background fields
Dz  = djles_diffmatrix(H/(NDATA-1), NDATA, 1, 'not periodic');
Dzz = djles_diffmatrix(H/(NDATA-1), NDATA, 2, 'not periodic');
rhozdata =Dz *rhodata;
ubgzdata =Dz *ubgdata;
ubgzzdata=Dzz*ubgdata;

% Now we build piecewise interpolating polynomials from the data,
% and convert them to function handles for the solver
method='linear'; % linear works well
% method='pchip';  % pchip also respects monotonicity
rho  = @(z) ppval(interp1(zdata,rhodata ,method,'pp'), z);
rhoz = @(z) ppval(interp1(zdata,rhozdata,method,'pp'), z);

Utarget   = @(z) ppval(interp1(zdata,ubgdata  ,method,'pp'), z);
Utargetz  = @(z) ppval(interp1(zdata,ubgzdata ,method,'pp'), z);
Utargetzz = @(z) ppval(interp1(zdata,ubgzzdata,method,'pp'), z);

%%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = clock;

% Now solve DJL, bringing in the background velocity incrementally
for alpha=linspace(0,1,4)
    % Velocity profiles
    Ubg   = @(z) alpha*Utarget(z);
    Ubgz  = @(z) alpha*Utargetz(z);
    Ubgzz = @(z) alpha*Utargetzz(z);
    
    % Use a reduced epsilon for these intermediate waves
    epsilon=1e-3;
    
    % Iterate the DJL solution
    djles_refine_solution
    djles_diagnostics; djles_plot; % uncomment to view progress at each step
end

% Increase resolution, restore default epsilon, iterate to convergence
NX=64; NZ=64; clear epsilon
djles_refine_solution
djles_diagnostics; djles_plot; % uncomment to view progress here

% Increase resolution, iterate to convergence
NX=128; NZ=128;
djles_refine_solution
djles_diagnostics; djles_plot; % uncomment to view progress here

% Increase to the final resolution, iterate to convergence
NX=512; NZ=256;
djles_refine_solution

end_time=clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

% Compute and plot the diagnostics
djles_diagnostics
djles_plot
