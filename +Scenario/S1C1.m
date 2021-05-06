%% Input parameters for Biomass-nutrient-flow solver
%% flags
flag.restart = 0;   % restart trigger
flag.adv = 0;       % advection trigger
flag.DR = 1;        % diffusion-reaction trigger

%% constraints
num.T = 1e2;        % max time in [h]
num.dtmax = 1e-3;   % max time step size [h]
num.tnmax = 1e3;    % max number of time steps
num.flow = 1;       % how often we solve for flow
num.tau = 1;        % number of minor time steps within each major time steps
num.rec = 50;       % how often we record data
num.plot = 50;      % how often we generate figure

%% unit conversion factors from [m,s,kg] to ...
unit.str_len = 'mm'; unit.length = 1e3; 
unit.str_time = 'h'; unit.time = 1./(3600);
unit.str_mass = 'g'; unit.mass = 1e3;

%% space in [m]
x.lo = 0;       x.hi = 100e-6;     x.n = 100;
y.lo = 0;       y.hi = 100e-6;     y.n = 100;

%% geometry
pm = 3;
switch pm
    case 0; porous_medium = 'channel';
    case 1; porous_medium = 'rock1';
    case 2; porous_medium = 'double_channel';
    case 3; porous_medium = 'corner_nonsymmetric';
    case 4; porous_medium = 'Symmetric';
    case 5; porous_medium = 'Nonsymmetric';
    case 6; porous_medium = 'porous1';
end

obst = 3;
switch obst
    case 0; obstacle = 'none';
    case 1; obstacle = 'obstacle1';
    case 2; obstacle = 'channel_biofilm2';
    case 3; obstacle = 'corner_nonsymmetric_biofilm';
    case 4; obstacle = 'Symmetric_biofilm';
    case 5; obstacle = 'Nonsymmetric_biofilm';
    case 6; obstacle = 'biofilm2';
end 

kb = 1e-9; % permeability of obstacle

%% parameters for flow and advection
% Dirichlet inlet, natural outflow, no-slip walls
flow.um = 1e-9;   % mean flow velocity
flow.uf = @(LB,UB,val,Um) -6.*Um.*(LB-val).*(UB-val)./(LB-UB).^2; % fully developed parabolic velocity
flow.dir = 3; % flow direction; 1: B->T, 2: T->B, 3: L->R, 4: R->L
flow.cfl = 0.95;
flow.mu = 8.9e-4; % dynamic viscosity of water [Pa-s]

%% parameters for diffusion-reaction
% biomass
b.kappa = .5;
b.star = 1;
b.nu = 0.9;
b.init = 0.8.*b.star;
b.eberla = 2;
b.dB0 = 1e-4.*unit.time./unit.length.^2;
b.eberlstar = 1.01;
% nutrient
n.kappa = 2.*unit.time;
n.inlet = 1;
if flag.adv == 1; n.init = 0; else; n.init = n.inlet; end
n.N0 = 1.18e-3;
n.df = 6.*unit.time./unit.length.^2;
n.db = 0.1.*n.df;