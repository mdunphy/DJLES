clear all

%%% Specify the parameters of the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Atarget = 2e-3; % APE for wave (m^4/s^2)
L  = 8.0;  % domain width (m)
H  = 0.2;  % domain depth (m)
NX = 32;   % grid
NZ = 32;   % grid

% The unitless density profile (normalized by a reference density rho0)
a_d=0.02; z0_d=0.05; d_d=0.01;
rho =@(z) 1-a_d*tanh((z+z0_d)/d_d);
rhoz=@(z) -(a_d/d_d)*sech((z+z0_d)/d_d).^2;

% The velocity profile (zero for this case) (m/s)
Ubg=@(z) 0*z; Ubgz=@(z) 0*z; Ubgzz=@(z) 0*z;

%%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = clock;

% Start with 1% of target APE, raise to target APE in 5 steps
for A=linspace(Atarget/100, Atarget, 5)
    % Find the solution of the DJL equation
    djles_refine_solution
end

% Increase the resolution, reduce epsilon, and iterate to convergence
NX=64; NZ=64; epsilon=1e-5;
djles_refine_solution

% Increase the resolution, and iterate to convergence
NX=512; NZ=256;
djles_refine_solution

end_time=clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

% Compute and plot the diagnostics
djles_diagnostics
djles_plot
