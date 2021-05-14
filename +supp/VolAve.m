function [Velu,Velv,Px,Py] = VolAve(x,y,U,V,P)
%% Volume averages
U0 = U; U0(isnan(U)) = 0; Velu = mean(U0(:));
V0 = V; V0(isnan(V)) = 0; Velv = mean(V0(:));
P0 = P; P0(isnan(P)) = 0;
PB = P0(1:y.n./2,:);        PB = mean(PB(:));
PT = P0(y.n./2+1:end,:);    PT = mean(PT(:));
PL = P0(:,1:x.n./2);        PL = mean(PL(:));
PR = P0(:,x.n./2+1:end);    PR = mean(PR(:));
Px = 2.*(PL-PR)./(x.hi-x.lo);
Py = 2.*(PB-PT)./(y.hi-y.lo);
end
