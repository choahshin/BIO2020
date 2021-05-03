function [U,V,P,nu,K] = Brinkman2D(x,y,G,flow,u,v,p)
% define parameters
hx = 1./x.h; hx2 = hx.^2;
hy = 1./y.h; hy2 = hy.^2;

% forcing term
fu = sparse(u.N,1); fv = sparse(v.N,1);

%% generate matrix
% u
Auux = kron(speye(y.N),spdiags([-1,2,-1].*ones(x.N+1,1),[-1:1],x.N+1,x.N+1));
Auuy = spdiags([-1,2,-1].*ones(u.N,1),[-x.N-1,0,x.N+1],u.N,u.N);
Auu = hx2.*Auux + hy2.*Auuy;

Bu = hx.*kron(speye(y.N),spdiags([-1,1].*ones(x.N+1,1),[-1,0],x.N+1,x.N));
Fu = fu(:);

% v
Avvx = kron(speye(y.N+1),spdiags([-1,2,-1].*ones(x.N,1),[-1:1],x.N,x.N));
Avvy = spdiags([-1,2,-1].*ones(v.N,1),[-x.N,x.N],v.N,v.N);
Avv = hx2.*Avvx + hy2.*Avvy;

Bv = hy.*spdiags([-1,1].*ones(v.N,1),[-x.N,0],v.N,p.N);
Fv = fv(:);

% p
C = sparse(p.N,p.N);
Fp = sparse(p.N,1);

%% boundary conditions
%%% u
% north
Auu(u.nind,u.nind) = Auu(u.nind,u.nind) + flow.bdry_id(1).*hy2.*speye(x.N+1);
% south
Auu(u.sind,u.sind) = Auu(u.sind,u.sind) + flow.bdry_id(3).*hy2.*speye(x.N+1);
% west
if flow.bdry_id(2) == -1 % outflow BC
    Auu(u.wind,u.wind) = Auu(u.wind,u.wind) - hx2.*speye(y.N);
else
    Auu(u.wind,:) = 0;
    Auu(u.wind,u.wind) = speye(y.N);
    Bu(u.wind,:) = 0;
    Fu(u.wind) = u.WBC(:);
end
% east
if flow.bdry_id(4) == -1
    Auu(u.eind,u.eind) = Auu(u.eind,u.eind) - hx2.*speye(y.N);
else
    Auu(u.eind,:) = 0;
    Auu(u.eind,u.eind) = speye(y.N);
    Bu(u.eind,:) = 0;
    Fu(u.eind) = u.EBC(:);
end
%%% v
% west
Avv(v.wind,v.wind) = Avv(v.wind,v.wind) + flow.bdry_id(2).*hx2.*speye(y.N+1);
% east
Avv(v.eind,v.eind) = Avv(v.eind,v.eind) + flow.bdry_id(4).*hx2.*speye(y.N+1);
% north
if flow.bdry_id(1) == -1
    Avv(v.nind,v.nind) = Avv(v.nind,v.nind) - hy2.*speye(x.N);
else
    Avv(v.nind,:) = 0;
    Avv(v.nind,v.nind) = speye(x.N);
    Bv(v.nind,:) = 0;
    Fv(v.nind) = v.NBC(:);
end
% south
if flow.bdry_id(3) == -1
    Avv(v.sind,v.sind) = Avv(v.sind,v.sind) - hy2.*speye(x.N);
else
    Avv(v.sind,:) = 0;
    Avv(v.sind,v.sind) = speye(x.N);
    Bv(v.sind,:) = 0;
    Fv(v.sind) = v.SBC(:);
end
%%% p
switch flow.inlet_id
    case 1; Fp(p.nind) = v.NBC(:).*hy;
    case 2; Fp(p.wind) = -u.WBC(:).*hx;
    case 3; Fp(p.sind) = -v.SBC(:).*hy;
    case 4; Fp(p.eind) = u.EBC(:).*hx;
end
Fp = Fp./flow.mu;
%% Darcy drag term
Ge.rock = sparse(x.N,y.N) + 1;
rock = G.rock;
Ge.K = sparse(x.N,y.N) + inf;
G.K = sparse(x.n,y.n) + inf;
G.K(G.b0'==1) = G.kb;
Ge.K(x.ind,y.ind) = G.K;
Ge.Ku = [Ge.K(1,:);2./(1./Ge.K(1:end-1,:)+1./Ge.K(2:end,:));Ge.K(end,:)];
Ge.Kv = [Ge.K(:,1),2./(1./Ge.K(:,1:end-1)+1./Ge.K(:,2:end)),Ge.K(:,end)];
Auu = Auu + spdiags(1./Ge.Ku(:),0,u.N,u.N);
Avv = Avv + spdiags(1./Ge.Kv(:),0,v.N,v.N);
Ge.rock(x.ind,y.ind) = rock;
Ge.u = [Ge.rock(1,:);min(Ge.rock(1:end-1,:),Ge.rock(2:end,:));Ge.rock(end,:)];
Ge.v = [Ge.rock(:,end),min(Ge.rock(:,1:end-1),Ge.rock(:,2:end)),Ge.rock(:,end)];

%% select unknowns
indu = find(Ge.u > 0);      Lu = length(indu);
indv = find(Ge.v > 0);      Lv = length(indv);
indp = find(Ge.rock > 0);   Lp = length(indp);

A = blkdiag(Auu(indu,indu),Avv(indv,indv));
B = [Bu(indu,indp);Bv(indv,indp)]./flow.mu;

S = [A,B;B',C(indp,indp)];
F = [Fu(indu);Fv(indv);Fp(indp)];

sol = S\F;
solc = mat2cell(sol,[Lu,Lv,Lp]);

Ue = sparse(u.N,1); Ue(indu) = solc{1}; Ue = reshape(Ue,x.N+1,y.N)';
Ve = sparse(v.N,1); Ve(indv) = solc{2}; Ve = reshape(Ve,x.N,y.N+1)';
Pe = sparse(p.N,1); Pe(indp) = solc{3}; Pe = reshape(Pe,x.N,y.N)';

U = Ue(y.ind,x.ind(1):x.ind(end)+1);
V = Ve(y.ind(1):y.ind(end)+1,x.ind);
P = Pe(y.ind,x.ind);

nu = max(abs(U(:)))./x.h + max(abs(V(:)))./y.h;

%% compute K
Pp = P; Pp(isnan(P)) = 0;
if flow.dir <= 2
    Vp = V; Vp(isnan(V)) = 0;
    Vel = mean(Vp(:));
    PB = Pp(1:y.n./2,:);     PB = mean(PB(:));
    PT = Pp(y.n./2+1:end,:); PT = mean(PT(:));
    grad = 2.*(PB-PT)./(y.hi-y.lo);
else
    Up = U; Up(isnan(U)) = 0;
    Vel = mean(Up(:));
    PL = Pp(:,1:x.n./2);     PL = mean(PL(:));
    PR = Pp(:,x.n./2+1:end); PR = mean(PR(:));
    grad = 2.*(PL-PR)./(x.hi-x.lo);
end
K = flow.mu.*Vel./grad;
end
    