clear all

%%% Specify the parameters of the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A  = 5e-5; % APE for wave (m^4/s^2)
L  = 3.0;  % domain width (m)
H  = 0.1;  % domain depth (m)

relax=0.15; % use strong underrelaxation

% Specify the general density profile which takes d_d as a second parameter
a_d=0.02; z0_d=0.025;
frho=@(z,d_d) 1-a_d*tanh((z+z0_d)/d_d);
frhoz=@(z,d_d) -(a_d/d_d)*sech((z+z0_d)/d_d).^2;

% The velocity profile (zero for this case) (m/s)
Ubg=@(z) 0*z; Ubgz=@(z) 0*z; Ubgzz=@(z) 0*z;

%%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

start_time = clock;

% Specify resolution and pycnocline parameter d_d according to this schedule
NXlist=[  64   128    256   256     512];
NZlist=[  32    64    128   128     256];
ddlist=[0.01 0.005 0.0025 0.001 0.00075];
for ddindex=1:length(ddlist)
    % Resolution for this wave
    NX = NXlist(ddindex);
    NZ = NZlist(ddindex);
    
    % Density profile for this wave with specified d_d
    d_d  = ddlist(ddindex);
    rho  = @(z) frho(z, d_d);
    rhoz = @(z) frhoz(z, d_d);
    
    % Iterate the DJL solution
    djles_refine_solution
    djles_diagnostics; djles_plot; % uncomment to view progress at each step
end

% Reduce epsilon, iterate to convergence
epsilon=1e-5;
djles_refine_solution

% Raise resolution, iterate to convergence
NX=2048; NZ=1024;
djles_refine_solution

end_time=clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

% Compute and plot the diagnostics
djles_diagnostics
djles_plot
