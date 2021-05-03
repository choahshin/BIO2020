%% diffusion-reaction with constraints
function [Bnew,Nnew,Lnew] = DR2D(x,y,dt,Bold,Nold,Lold,G,b,n,flag,flow)
% B = biomass
% N = nutrient
% L = Lagrange multiplier

%% reorient
B0 = Bold';
N0 = Nold';
L0 = Lold';
%% diffusion coefficients
% biomass (Eberl)
Bbar = B0./b.star;
G.dB = b.dB0.*(1+Bbar./(b.eberlstar -Bbar)).^b.eberla;
% arithmetic average
G.dBu = [G.dB(1,:);(G.dB(1:end-1,:) + G.dB(2:end,:))./2;G.dB(end,:)];
G.dBv = [G.dB(:,1),(G.dB(:,1:end-1) + G.dB(:,2:end))./2,G.dB(:,end)];
% assign 0 for no-flow cell edges
G.dBu(G.u == 0) = 0;
G.dBv(G.v == 0) = 0;
if flag.adv == 1
    if flow.dir <= 2; G.dBu([1,end],:) = 0; else; G.dBv(:,[1,end]) = 0; end
else
    G.dBu([1,end],:) = 0; G.dBv(:,[1:end]) = 0;
end
% nutrient
G.dN = n.df.*(b.star-B0)+ n.db.*B0;
% G.dN = n.df.*G.rock;
% G.dN(G.b0' == 1) = n.db;
G.dN(G.rock_id) = 0;
% harmonic average
G.dNu = [G.dN(1,:);2./(1./G.dN(1:end-1,:)+1./G.dN(2:end,:));G.dN(end,:)];
G.dNv = [G.dN(:,1),2./(1./G.dN(:,1:end-1)+1./G.dN(:,2:end)),G.dN(:,end)];
if flag.adv == 1
    if flow.dir <= 2; G.dNu([1,end],:) = 0; else; G.dNv(:,[1,end]) = 0; end
end
%% stiffness matrices
Np = x.n.*y.n; hx2 = 1./x.h.^2; hy2 = 1./y.h.^2;
% index
nind = [Np - x.n+1:Np]'; sind = [1:x.n]';
wind = [1:x.n:Np]';      eind = [x.n:x.n:Np]';
% -div(d_B grad(B)) = (b.A)B
% -dif(d_N grad(N)) = (n.A)N
b.Au = []; n.Au = [];
for j = 1:y.n
    dbtu = G.dBu(:,j); dbtu = dbtu(:);
    bAut = spdiags([-dbtu(2:end),dbtu(1:end-1)+dbtu(2:end),-dbtu(1:end-1)],...
                   [-1:1],x.n,x.n);
    b.Au = blkdiag(b.Au,bAut);
    dntu = G.dNu(:,j); dntu = dntu(:);
    nAut = spdiags([-dntu(2:end),dntu(1:end-1)+dntu(2:end),-dntu(1:end-1)],...
                   [-1:1],x.n,x.n);
    n.Au = blkdiag(n.Au,nAut);
end

b.F = sparse(Np,1);
n.F = sparse(Np,1);

dbtv = G.dBv(:);
b.Av = spdiags([-dbtv(x.n+1:end),dbtv(1:end-x.n)+dbtv(x.n+1:end),-dbtv(1:end-x.n)],...
               [-x.n,0,x.n],Np,Np);
dntv = G.dNv(:);
n.Av = spdiags([-dntv(x.n+1:end),dntv(1:end-x.n)+dntv(x.n+1:end),-dntv(1:end-x.n)],...
               [-x.n,0,x.n],Np,Np);
%% boundary conditions
% north
if b.ind(1) == -1
    b.Av(nind,nind) = b.Av(nind,nind) - G.dBv(:,end).*speye(x.n);
else
    b.F(nind) = b.F(nind) + G.dBv(:,end).*hy2.*b.NBC(:);
end
if n.ind(1) == -1
    n.Av(nind,nind) = n.Av(nind,nind) - G.dNv(:,end).*speye(x.n);
else
    n.F(nind) = n.F(nind) + G.dNv(:,end).*hy2.*n.NBC(:);
end
% west
if b.ind(2) == -1
    b.Av(wind,wind) = b.Au(wind,wind) - G.dBu(1,:)'.*speye(y.n);
else
    b.F(wind) = b.F(wind) + G.dBu(1,:)'.*hx2.*b.WBC(:);
end
if n.ind(2) == -1
    n.Av(wind,wind) = n.Au(wind,wind) - G.dNu(1,:)'.*speye(y.n);
else
    n.F(wind) = n.F(wind) + G.dNu(1,:)'.*hx2.*n.WBC(:);
end
% south
if b.ind(3) == -1
    b.Av(sind,sind) = b.Av(sind,sind) - G.dBv(:,1).*speye(x.n);
else
    b.F(sind) = b.F(sind) + G.dBv(:,1).*hy2.*b.SBC(:);
end
if n.ind(3) == -1
    n.Av(sind,sind) = n.Av(sind,sind) - G.dNv(:,1).*speye(x.n);
else
    n.F(sind) = n.F(sind) + G.dNv(:,1).*hy2.*n.SBC(:);
end
% east
if b.ind(4) == -1
    b.Av(eind,eind) = b.Au(eind,eind) - G.dBu(end,:)'.*speye(y.n);
else
    b.F(eind) = b.F(eind) + G.dBu(end,:)'.*hx2.*b.EBC(:);
end
if n.ind(4) == -1
    n.Av(eind,eind) = n.Au(eind,eind) - G.dNu(end,:)'.*speye(y.n);
else
    n.F(eind) = n.F(eind) + G.dNu(end,:)'.*hx2.*n.EBC(:);
end
% combine
b.A = b.Au.*hx2 + b.Av.*hy2;
n.A = n.Au.*hx2 + n.Av.*hy2;
b.F = dt.*b.F + B0(:);
n.F = dt.*n.F + N0(:);

%% nonlinear equation F([B;N;L]) = 0
% monod function
g = @(nval) n.kappa.*nval./(nval+n.N0);
dg = @(nval) n.kappa.*n.N0./(nval+n.N0).^2;
% (I+dt A)
b.AA = speye(Np) + dt.*b.A;
n.AA = speye(Np) + dt.*n.A;

% select unknowns
np = length(G.flow_id);
b.AA = b.AA(G.flow_id,G.flow_id); n.AA = n.AA(G.flow_id,G.flow_id);
b.F = b.F(G.flow_id);             n.F = n.F(G.flow_id);
B1 = B0(G.flow_id);     N1 = N0(G.flow_id);     L1 = L0(G.flow_id);
% nonlinear function F
f1 = @(bnl) b.AA*bnl(1:np) - dt.*b.kappa.*bnl(1:np).*...
            g(bnl(np+1:2.*np)) - dt.*bnl(end-np+1:end) - b.F;
f2 = @(bnl) n.AA*bnl(np+1:2.*np) +dt.*bnl(1:np).*...
            g(bnl(np+1:2.*np)) - n.F;
f3 = @(bnl) max(bnl(1:np)-b.star,bnl(end-np+1:end));
f = @(bnl) x.h.*y.h.*[f1(bnl);f2(bnl);f3(bnl)];

% Jacobian
j11 = @(bnl) b.AA - dt.*b.kappa.*g(bnl(np+1:2.*np)).*speye(np);
j12 = @(bnl) -dt.*b.kappa.*bnl(1:np).*dg(bnl(np+1:2.*np)).*speye(np);
j13 = @(bnl) -dt.*speye(np);
j21 = @(bnl) dt.*g(bnl(np+1:2.*np)).*speye(np);
j22 = @(bnl) n.AA + dt.*bnl(1:np).*dg(bnl(np+1:2.*np)).*speye(np);
j23 = @(bnl) sparse(np,np);
j31 = @(bnl) spdiags([(bnl(1:np)-b.star >= bnl(end-np+1:end))],0,np,np);
j32 = @(bnl) sparse(np,np);
j33 = @(bnl) spdiags([(bnl(1:np)-b.star < bnl(end-np+1:end))],0,np,np);
j = @(bnl) x.h.*y.h.*[j11(bnl), j12(bnl), j13(bnl);...
            j21(bnl), j22(bnl), j23(bnl);...
            j31(bnl), j32(bnl), j33(bnl)];
% semi-smooth Newton
iter = 0;
max_iter = 10;
bnl0 = [B1(:);N1(:);L1(:)];
abs_res = norm(f(bnl0),inf);
abs_tol = min(1e-4.*abs_res,1e-8);
rel_res = inf;
rel_tol = 1e-4;
while iter < max_iter && (abs_res > abs_tol && rel_res > rel_tol)
    iter = iter + 1;
    bnl = bnl0 - j(bnl0)\f(bnl0);
    abs_res = norm(f(bnl),inf);
    rel_res = norm(bnl-bnl0,inf);
    bnl0 = bnl;
    fprintf('\t Semismooth Newton: iter = %d, abs_res = %g, rel_res = %g \n',iter,abs_res,rel_res);
end  

bnlsol = mat2cell(bnl0,[np,np,np]);
Bnew = sparse(Np,1); Bnew(G.flow_id) = bnlsol{1}; Bnew = reshape(Bnew,x.n,y.n)';
Nnew = sparse(Np,1); Nnew(G.flow_id) = bnlsol{2}; Nnew = reshape(Nnew,x.n,y.n)';
Lnew = sparse(Np,1); Lnew(G.flow_id) = bnlsol{3}; Lnew = reshape(Lnew,x.n,y.n)';
end