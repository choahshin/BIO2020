%% Upwind2D
function Qnew = Upwind2D_B(x,y,dt,U,V,Qold,b)
    Qnew = sparse(y.n,x.n);
    Qstar = Qnew;
    Qbar = Qold;
    Qbar(Qold >= b.star2) = 0;
    % x-direction
    for j = 1:y.n
        Qstar(j,:) = supp.Upwind1D(x,dt,Qold(j,:),Qbar(j,:),U(j,:),[b.ind(2),b.ind(4)],...
                           [b.WBC(j),b.EBC(j)],1);
    end
    % y-direction
    for i = 1:x.n
        Qnew(:,i) = supp.Upwind1D(y,dt,Qstar(:,i),Qbar(:,i),V(:,i),[b.ind(3),b.ind(1)],...
                          [b.SBC(i),b.NBC(i)],2);
    end
end