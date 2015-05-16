clear all

%%% Specify the parameters of the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A  = 1e-4; % APE for wave (m^4/s^2)
L  = 8.0;  % domain width (m)
H  = 0.2;  % domain depth (m)
NX = 32;   % grid
NZ = 32;   % grid

% The unitless density profile (normalized by a reference density rho0)
a_d=0.02; z0_d=0.05; d_d=0.01;
rho =@(z) 1-a_d*tanh((z+z0_d)/d_d);
rhoz=@(z) -(a_d/d_d)*sech((z+z0_d)/d_d).^2;

% Specify general velocity profile that takes U0 as a second parameter (m/s)
zj=0.5*H; dj=0.4*H;
fUbg  =@(z,U0) U0*tanh((z+zj)/dj);
fUbgz =@(z,U0) (U0/dj)*sech((z+zj)/dj).^2;
fUbgzz=@(z,U0) (-2*U0/(dj*dj)).*(sech((z+zj)/dj).^2).*tanh((z+zj)/dj);

%%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = clock;

% Start with U0=0, raise it to U0=0.1 over 6 increments
for U0=linspace(0, 0.1, 6)
    % Velocity profile for this wave
    Ubg  =@(z) fUbg(z,U0);
    Ubgz =@(z) fUbgz(z,U0);
    Ubgzz=@(z) fUbgzz(z,U0);
    
    % Use a reduced epsilon for these intermediate waves
    epsilon=1e-3;
    
    % Find the solution of the DJL equation
    djles_refine_solution
    %djles_diagnostics; djles_plot; % uncomment to view progress at each step
end

% Increase the resolution, reduce epsilon, iterate to convergence
NX=64; NZ=64; epsilon=1e-6;
djles_refine_solution

% Increase the resolution and iterate to convergence
NX=512; NZ=256;
djles_refine_solution

end_time=clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

% Compute and plot the diagnostics
djles_diagnostics
djles_plot
