function [K] = Flow(ns,nc,nd,ne)
%=========================================================================%
% Solves 2d flow models
%-------------------------------------------------------------------------%
switch nargin
    case 2; nd = 3; ne = 1; 
    case 3; ne = 1;
end
nf = 1;
%% import data
supp.input_data;
start_clock = tic;

B0 = b.IC;
N0 = n.IC;
L0 = 0*B0;
tn = 0; 
t = 0;
    
G.b0 = sparse(y.n,x.n);
G.b0(B0 >= b.star2 & B0 <= b.star) = 1;

supp.UpscaledK;

end_clock = toc(start_clock);
close('all','hidden');

diary(sprintf('%s/record.txt',dir.version));
diary on
if G.kb > 0; fprintf('Heterogeneous Brinkman Flow\n');
else; fprintf('Stokes flow\n'); end
if ne == 1
    fprintf(sprintf('K = %.4e %s^2\n',full(K),unit.str_len));
else
    fprintf(sprintf('K = [%.4e,%.4e;%.4e,%.4e] %s^2\n',...
            K(1,1),K(1,2),K(2,1),K(2,2),unit.str_len));
end
fprintf('Total clock time = %.4f minutes\n',end_clock./60);
fprintf('========================================\n\n');
diary off
end