%% Solver
close('all','hidden');
clear all;

%% choose your senario (ns) and case (nc)
ns = 0;
nc = 0;

%% import data
supp.input_data;

start_clock = tic;    
%% initialize
T = num.T;
if flag.restart == 1
    tn = 0;
    while isfolder(sprintf('%s/plt%d',dir.version,tn)); tn = tn + 1; end
    tn = tn - 1;
    dir.tn = sprintf('%s/plt%d',dir.version,tn);
    B0 = load(sprintf('%s/B.txt',dir.tn));
    N0 = load(sprintf('%s/N.txt',dir.tn));
    L0 = load(sprintf('%s/L.txt',dir.tn));
    t = load(sprintf('%s/T.txt',dir.tn));
else
    B0 = b.IC;
    N0 = n.IC;
    L0 = 0*B0;
    tn = 0;
    t = 0;
end

G.b0 = sparse(y.n,x.n);
G.b0(B0 >= b.star2 & B0 <= b.star) = 1;
G.bn = G.b0;
flag.Gb = [];

fprintf('tn = %d, t = %g [%s]\n', tn, t, unit.str_time);

%% iteration
if G.kb > 0
    supp.BN_Brinkman;
else
    supp.BN_Stokes;
end

%% record at final time T

end_clock = toc(start_clock);
close('all','hidden');

diary(sprintf('%s/recoard.txt',dir.version));
diary on
fprintf('Total clock time = %.4f minutes\n',end_clock./60);
fprintf('========================================\n\n');
diary off

supp.record;
supp.figGen;
