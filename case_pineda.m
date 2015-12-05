clear
% Load density data. The text file contains an approximate reproduction of
% the smoothed curve in Figure 9 (a) from Pineda et al (2015).
% Courtesy of Jorge Magalhães and José da Silva.
data = load('case_pineda_cast1.txt'); zdata = data(:,1); rhodata = data(:,2);

rho0 = max(rhodata);  % Used in djles_common to compute N2(z)

% Use MATLAB's gradient function to find the derivative drho/dz
rhozdata = gradient(rhodata,zdata);

% Now build piecewise interpolating polynomials from the data,
% and convert them to function handles for the solver
method='pchip';  % pchip respects monotonicity
warning('off','MATLAB:interp1:ppGriddedInterpolant'); % Silences R2014a warnings
rho  = @(z) ppval(interp1(zdata,rhodata ,method,'pp'), z);
rhoz = @(z) ppval(interp1(zdata,rhozdata,method,'pp'), z);

% The velocity profile (zero for this case) (m/s)
Ubg=@(z) 0*z; Ubgz=@(z) 0*z; Ubgzz=@(z) 0*z;

L  =   1200;  % domain width (m)
H  =     57;  % domain depth (m), estimated from Pineda et al's Figure 9 (a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find the wave showcased in Pineda et al. (2015) Figure 11 %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = clock;

% Set initial resolution and large epsilon for intermediate waves
NX=32; NZ=32; epsilon=1e-3;

% Raise amplitude in a few increments
for A=linspace(1e4, 3.62e5, 6)  % APE (kg m/s^2)
    djles_refine_solution
end

% Increase resolution, reduce epsilon, iterate to convergence
NX=64; NZ=64; epsilon=1e-4;
djles_refine_solution

NX=128; NZ=128; epsilon=1e-5;
djles_refine_solution

NX=256; NZ=256; epsilon=1e-6;
djles_refine_solution

NX=512; NZ=512; epsilon=1e-7;
djles_refine_solution

end_time=clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

% Compute diagnostics, plot wave
djles_diagnostics
djles_plot
djles_pressure

% Construct Pineda et al. (2015) Figure 11
figure(11)
subplot(4,1,1:2)
contour(XC-L/2,ZC,density,[1022.3:0.3:1024.7],'k')
xlim([-1 1]*300); title('Density contours')

subplot(4,1,3)
plot(xc-L/2,interp2(XC,ZC,u,xc,-H+1),'k')
xlim([-1 1]*300); grid on; title('U at 1 mab');

subplot(4,1,4)
plot(xc-L/2,interp2(XC,ZC,p,xc,-H+1),'k')
xlim([-1 1]*300); grid on; title('Pressure minus background pressure')

clear eta c

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Find the solid curves shown in Pineda et al. (2015) Figure 10 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = clock;

verbose=0;

Alist=logspace(log10(1e3), log10(1e6), 11);
WArec=zeros(size(Alist));   % wave amplitude
Crec=zeros(size(Alist));    % wave phase speed
Urec=zeros(size(Alist));    % velocity at 1 metre above bottom
Prec=zeros(size(Alist));    % pressure at 1 metre above bottom

% Solve the wave at each amplitude in Alist
% Each wave is re-used as initial guess for subsequent waves
for ai=1:length(Alist)
    A=Alist(ai);

    % Change resolution, reduce epsilon, iterate to convergence
    NX=32; NZ=32; epsilon=1e-3;
    djles_refine_solution

    NX=64; NZ=64; epsilon=1e-4;
    djles_refine_solution

    NX=128; NZ=128; epsilon=1e-5;
    djles_refine_solution

    % Compute and record quantities
    djles_diagnostics
    djles_pressure

    WArec(ai) = wave_ampl;
    Crec(ai) = c;

    u1mab = interp2(XC,ZC,u,xc,-H+1);
    [~,idx] = max(abs(u1mab));
    Urec(ai) = u1mab(idx);

    p1mab = interp2(XC,ZC,p,xc,-H+1);
    [~,idx] = max(abs(p1mab));
    Prec(ai) = p1mab(idx);
end
end_time=clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

% Construct Pineda et al. (2015) Figure 10
figure(10)
subplot(3,1,1)
plot(-WArec, Crec,'k'); xlim([0 23]); grid on; title('c (m/s)');

subplot(3,1,2)
plot(-WArec, Urec,'k'); xlim([0 23]); grid on; title('U at 1 mab (m/s)');

subplot(3,1,3)
plot(-WArec, Prec,'k'); xlim([0 23]); grid on; title('P at 1 mab (Pa)');
