%% load values from input card
inputcard = sprintf('S%dC%d.m',ns,nc);
run(sprintf('Scenario.%s',inputcard));
%% folder set up ==========================================================
if flag.restart == 1
    dir.scenario = sprintf('Simulation/Scenario%d',ns);
    dir.case = sprintf('%s/Case%d',dir.scenario,nc);
    ver = 1;
    while isfolder(sprintf('%s/V%d',dir.case,ver)); ver = ver + 1; end
    dir.version = sprintf('%s/V%d',dir.case,ver-1);
else
    if not(isfolder('Simulation')); mkdir('Simulation'); end
    dir.scenario = sprintf('Simulation/Scenario%d',ns);
    if not(isfolder(dir.scenario)); mkdir(dir.scenario); end
    dir.case = sprintf('%s/Case%d',dir.scenario,nc);
    if not(isfolder(dir.case)); mkdir(dir.case); end
    ver = 1;
    while isfolder(sprintf('%s/V%d',dir.case,ver)); ver = ver + 1; end
    dir.version = sprintf('%s/V%d',dir.case,ver);
    mkdir(dir.version);
end
% save a copy of inputcard into recording folder
copyfile(sprintf('+Scenario/%s',inputcard),sprintf('%s/',dir.version));
%% unit conversion ========================================================
x.lo = x.lo.*unit.length; x.hi = x.hi.*unit.length;
y.lo = y.lo.*unit.length; y.hi = y.hi.*unit.length;
%
G.kb = kb.*unit.length.^2;
%
flow.um = flow.um.*unit.length./unit.time;
flow.mu = flow.mu.*unit.mass./(unit.length.*unit.time);
% 
% b.star = b.star.*unit.mass./unit.length.^3;
b.star2 = b.nu.*b.star;
% b.init = b.init.*unit.mass./unit.length.^3;
b.dB0 = b.dB0.*unit.length.^2./unit.time;
%
n.kappa = n.kappa./unit.time;
% n.inlet = n.inlet.*unit.mass./unit.length.^3;
% n.init = n.init.*unit.mass./unit.length.^3;
% n.N0 = n.N0.*unit.mass./unit.length.^3;
n.df = n.df.*unit.length.^2./unit.time;
n.db = n.db.*unit.length.^2./unit.time;
%% generate mesh ==========================================================
% extended domain
x.nw = 0; y.nw = 0;
% if flag.adv == 1 && (pm > 1 || (obst > 0 && G.kb == 0))
%     if flow.dir <= 2; y.nw = floor(y.n./2); else; x.nw = floor(x.n./2); end
% end
% mesh size
x.h = (x.hi-x.lo)./x.n;     y.h = (y.hi-y.lo)./y.n;
% indices for domain of interest
x.ind = [x.nw+1:x.nw+x.n]'; y.ind = [y.nw+1:y.nw+y.n]';
% width of extr flow region
x.w = x.h.*x.nw; y.w = y.h.*y.nw;
% dimension of full medium
x.a = x.lo - x.w; x.b = x.hi + x.w; x.L = x.b - x.a;
y.a = y.lo - y.w; y.b = y.hi + y.w; y.L = y.b - y.a;
% total number of mesh
x.N = x.n + 2.*x.nw; y.N = y.n + 2.*y.nw;
% discretize x- and y- intevals
x.e = [x.a:x.h:x.b]'; x.ce = 0.5.*(x.e(1:end-1)+x.e(2:end)); x.c = x.ce(x.ind);
y.e = [y.a:y.h:y.b]'; y.ce = 0.5.*(y.e(1:end-1)+y.e(2:end)); y.c = y.ce(y.ind);
% generate 2D mesh at cell centers
[x.p,y.p] = meshgrid(x.c,y.c);
%% flow boundary conditions ===============================================
% [north, west, south, east]
if flag.adv == 1
    p.N = x.N.*y.N;
    p.nind = [p.N-x.N+1:p.N];   p.sind = [1:x.N]';
    p.wind = [1:x.N:p.N]';      p.eind = [x.N:x.N:p.N]';
    u.N = (x.N+1).*y.N;         v.N = x.N.*(y.N+1);
    u.nind = [u.N-x.N:u.N]';    v.nind = [v.N-x.N+1:v.N]';
    u.sind = [1:x.N+1]';        v.sind = [1:x.N]';
    u.wind = [1:x.N+1:u.N]';    v.wind = [1:x.N:v.N]';
    u.eind = [x.N+1:x.N+1:u.N]';v.eind = [x.N:x.N:v.N]';
    u.NBC = sparse(x.N+1,1);    v.NBC = sparse(x.N,1);
    u.WBC = sparse(y.N,1);      v.WBC = sparse(y.N+1,1);
    u.SBC = sparse(x.N+1,1);    v.SBC = sparse(x.N,1);
    u.EBC = sparse(y.N,1);      v.EBC = sparse(y.N+1,1);
    switch flow.dir
        case 1 % bottom to top
            flow.bdry_id =[-1, 1, 1, 1];
            flow.inlet_id = 3;
            v.SBC = flow.uf(x.a,x.b,x.ce,flow.um);
        case 2 % top to bottom
            flow.bdry_id =[1, 1, -1, 1];
            flow.inlet_id = 1;
            v.NBC = -flow.uf(x.a,x.b,x.ce,flow.um);
        case 3 % left to right
            flow.bdry_id =[1, 1, 1, -1];
            flow.inlet_id = 2;
            u.WBC = flow.uf(y.a,y.b,y.ce,flow.um);
        case 4 % right to left
            flow.bdry_id =[1, -1, 1, 1];
            flow.inlet_id = 4;
            u.EBC = -flow.uf(y.a,y.b,y.ce,flow.um);
    end
end
%% interpret geometry =====================================================
% read image
% G.rock = sparse(x.n,y.n)+1;
G.rock = flipud(round(im2double(rgb2gray(imresize(imread(...
            sprintf('Porous_Medium/%s.jpg',porous_medium)),[y.n,x.n])))))';
OB = flipud(round(im2double(rgb2gray(imresize(imread(...
            sprintf('Obstacle/%s.jpg',obstacle)),[y.n,x.n])))))';    
% indices for rock (0), obstacle(2), flow(1)        
G.omega = sparse(x.n,y.n)+1;
G.omega(OB == 0) = 2;
G.omega(G.rock == 0) = 0;
G.rock_id = find(G.omega == 0);
G.flow_id = find(G.omega > 0);
G.u = [G.rock(1,:); 2./(1./G.rock(1:end-1,:)+1./G.rock(2:end,:));G.rock(end,:)];
G.v = [G.rock(:,1), 2./(1./G.rock(:,1:end-1)+1./G.rock(:,2:end)),G.rock(:,end)];
%% Biomass and nutrient ===================================================
% initial conditions
b.IC = 0*G.omega; b.IC(G.omega == 2) = b.init; b.IC = b.IC';
n.IC = sparse(y.n,x.n) + n.init;
% boundary conditions
b.NBC = 0*x.c; b.SBC = 0*x.c; b.WBC = 0*y.c; b.EBC = 0*y.c;
n.NBC = 0*x.c; n.SBC = 0*x.c; n.WBC = 0*y.c; n.EBC = 0*y.c;
b.NBC = 0*x.c; b.SBC = 0*x.c; b.WBC = 0*y.c; b.EBC = 0*y.c;
n.NBC = 0*x.c; n.SBC = 0*x.c; n.WBC = 0*y.c; n.EBC = 0*y.c;
b.ind = -ones(1,4); n.ind = -ones(1,4); % BC indicator
if flag.adv == 1
    switch flow.dir
        case 1; b.ind(3) = 1; n.ind(3) = 1; n.SBC = n.SBC + n.inlet;        
        case 2; b.ind(1) = 1; n.ind(1) = 1; n.NBC = n.NBC + n.inlet;
        case 3; b.ind(2) = 1; n.ind(2) = 1; n.WBC = n.WBC + n.inlet;
        case 4; b.ind(4) = 1; n.ind(4) = 1; n.EBC = n.EBC + n.inlet;
    end
else
    n.ind = -n.ind;
    n.NBC = n.NBC + n.inlet; n.SBC = n.SBC + n.inlet;
    n.WBC = n.WBC + n.inlet; n.EBC = n.EBC + n.inlet;
end
%% color code
cc.purple = [128 0 128]./255;
cc.purple2 = [204 0 204]./255;
cc.gray = [128 128 128]./255;
%% Record data ============================================================
diary(sprintf('%s/record.txt',dir.version));
diary on
fprintf('Scenario %d case %d: ',ns,nc);
if flag.adv == 1; fprintf('Adv-Diff-Reac\n'); else; fprintf('Diff-Reac\n'); end
fprintf('Units length [%s], time [%s], mass [%s]\n',unit.str_len,unit.str_time,unit.str_mass);
fprintf('Geometry --------------------------------\n');
fprintf('Omega = (%g,%g) x (%g,%g)\n', x.lo,x.hi,y.lo,y.hi);
fprintf('[Nx,Ny] = [%d,%d]\n',x.n,y.n);
fprintf('Porous medium = %s\n', porous_medium);
fprintf('Obstacle = %s\n',obstacle);
fprintf('kb = %.4e [%s^2]\n',G.kb,unit.str_len);
if flag.adv == 1
    fprintf('Flow ------------------------------------\n');
    fprintf('mu = %.4e\n',flow.mu);
    fprintf('flow direction = %d\n',flow.dir);
    fprintf('mean velocity = %.2e\n',flow.um);
end
% fprintf('Biomass ---------------------------------\n');
% fprintf('B* = %g\n', b.star);
% fprintf('B_* = %g\n',b.star2);
% fprintf('B_o = %g\n',b.init);
% fprintf('d_B = %.4e\n',b.dB0);
% fprintf('alpha = %g\n',b.eberla);
% fprintf('Utilitzation parameter = %.2e\n',b.kappa);
% fprintf('Nutrient --------------------------------\n');
% fprintf('N_inlet = %g\n', n.inlet);
% fprintf('N_o = %g\n',n.init);
% fprintf('Dimensionless factor in g(N) = %.2e\n',n.kappa);
% fprintf('Monod constant = %g\n',n.N0);
% fprintf('recording -------------------------------\n');
% fprintf('num.flow = %d\n',num.flow);
% fprintf('num.tau = %d\n',num.tau);
% fprintf('num.plot = %d\n',num.plot);
% fprintf('max_tn = %d\n',num.tnmax);
% fprintf('max_dt = %d\n',num.dtmax);
% fprintf('restart = %d\n',flag.restart);
% fprintf('advection = %d\n',flag.adv);
% fprintf('diffusion-reaction = %d\n',flag.DR);
fprintf('=========================================\n');
diary off