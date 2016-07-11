clear all

%%% Load the Amazon data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data      = load('case_amazon.txt'); % columns are z,rho,u
zdata     = data(:,1);
rhodata0  = data(:,2);
udata     = data(:,3);

% Compute derivatives
rhozdata0 = gradient(rhodata0,zdata);
uzdata    = gradient(udata,zdata);
uzzdata   = gradient(uzdata,zdata);

% Density adjustments (to ensure that drho/dz is negative)
rhodata=rhodata0;                   % copy original density
rhodata(end)=rhodata(end-1)*0.9999; % replace surface value to ensure positive N^2

rb=0.5*mean(rhozdata0(zdata>-55 & zdata<-30));  % define a mean density gradient from z=(-55,-30) m
ri=find(zdata<=-55);                            % find indexes for z <= -55 m
rhodata(ri)=rb*(zdata(ri)- zdata(ri(end)+1)) + rhodata(ri(end)+1); % replace data below z=-55 m with linear density
rhozdata = gradient(rhodata,zdata);             % adjusted density derivative

% Show original (blue) and adjusted (red) density and first derivative
figure(2)
subplot(1,2,1); plot(rhodata0 ,zdata,rhodata ,zdata,'r','linewidth',2); grid on
subplot(1,2,2); plot(rhozdata0,zdata,rhozdata,zdata,'r','linewidth',2); grid on

% Build piecewise interpolating polynomials from the data
warning('off','MATLAB:interp1:ppGriddedInterpolant'); % Silences R2014a warnings
method='linear';
rho    = @(z) ppval(interp1(zdata,rhodata ,method,'pp'), z);
rhoz   = @(z) ppval(interp1(zdata,rhozdata,method,'pp'), z);
fUbg   = @(z) ppval(interp1(zdata,udata   ,method,'pp'), z);
fUbgz  = @(z) ppval(interp1(zdata,uzdata  ,method,'pp'), z);
fUbgzz = @(z) ppval(interp1(zdata,uzzdata ,method,'pp'), z);

%%%%%%% Now we have data, prepare for DJLES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L       = 600;  % domain width (m)
H       = 80;   % domain depth (m)
rho0    = 1000; % reference density (kg/m^3)
verbose = 1;

start_time = clock;
% Start at low resolution/epsilon, no background velocity, and raise amplitude
Ubg=@(z) 0*z; Ubgz=@(z) 0*z; Ubgzz=@(z) 0*z; 
NX = 32; NZ = 32; epsilon=1e-3;
A1=1e5; A2=1e6;
for A=linspace(A1, A2, 3)
    djles_refine_solution
end

% Improve resolution, raise amplitude
NX=64; NZ=64;
A1=1e6; A2=5e6;
for A=linspace(A1, A2, 3)
    djles_refine_solution
end

% Now bring in the background velocity incrementally
for alpha=linspace(0.25,1,4)
    % Velocity profile for this wave
    Ubg  =@(z) alpha*fUbg(z);
    Ubgz =@(z) alpha*fUbgz(z);
    Ubgzz=@(z) alpha*fUbgzz(z);

    djles_refine_solution
end

% Improve resolution and epsilon, raise amplitude
NX=128; NZ=128; epsilon=1e-4; relax=0.4;
A1=5e6; A2=7.5e6;
for A=linspace(A1, A2, 4)
    djles_refine_solution
end

% Final wave: increase resolution, iterate to convergence
NX=256; NZ=256; epsilon=1e-6;
djles_refine_solution

% Compute diagnostics, plot wave
djles_diagnostics
djles_plot

% The wave is now near its amplitude limit and is beginning to transition
% to a broad flat crested wave. Use a larger domain for larger vales of APE.

end_time=clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));
