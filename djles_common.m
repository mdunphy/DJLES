% djles_common.m
% This file sets some default parameters and generates the grid/wavenumbers.

%%% Specify the default  parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting any of these options in the case file overrides the value here.

% Default min and max number of iterations for the iterative procedure.
if ~exist('min_iteration','var'), min_iteration=10; end
if ~exist('max_iteration','var'), max_iteration=2000; end

% Default number of Legendre points for Gauss quadrature
if ~exist('NL','var'), NL=20; end

% Underrelaxation factor, 0 < relax <= 1
% Setting this is often helps in finding a solution and/or speeding
% convergence. This value is interpreted as the fraction of the new
% field to keep after an iteration completes, that is,
%   val = (1-relax)*oldval + (relax)*newval
% Setting relax=0 prevents the iterations from proceeding, and setting
% relax=1 simply disables underrelaxation.
% Values between 0.3 and 0.8 seem to work fairly well.
if ~exist('relax','var'), relax=0.5; end

% gravitational constant (m/s^2)
if ~exist('g','var'), g=9.81; end

% Verbose flag, set to 1 to see the progress during solving
if ~exist('verbose','var'), verbose=1; end

% Convergence criteria: Stop iterating when the relative difference between
% successive iterations differs by less than epsilon
if ~exist('epsilon','var'), epsilon=1e-4; end

%%% Generate the grid and wavenumbers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare the grids: (xc,zc) are cell centred, (xe,ze) are cell edges
dx=L/NX; dz=H/NZ;
xc = (0.5:NX)*dx; zc = (0.5:NZ)*dz -H; [XC, ZC] = meshgrid(xc,zc);
xe = (0:NX)*dx;   ze = (0:NZ)*dz -H;   [XE, ZE] = meshgrid(xe,ze);

% Get Legendre points and weights for Gauss quadrature
[zl,wl]=gauss(NL); zl=(zl+1)/2; wl=wl/2;

% Get 2D weights for sine rectangle quadrature
wsine = djles_sinequadrature(NX,NZ,H,L);

% Create the wavenumber vectors and inverse Laplacian operator
ks = (pi/L) * [0:(NX-1) -NX:-1];
ms = (pi/H) * [0:(NZ-1) -NZ:-1]';
INVLAP = -1./bsxfun(@plus,ms.^2,ks.^2);
INVLAP(isinf(INVLAP))=0;
