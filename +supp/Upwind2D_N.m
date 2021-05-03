%% Upwind2D
function Qnew = Upwind2D_N(x,y,dt,U,V,Qold,n)
    Qnew = sparse(y.n,x.n);
    Qstar = Qnew;
    Qbar = Qold;
    % x-direction
    for j = 1:y.n
        Qstar(j,:) = supp.Upwind1D(x,dt,Qold(j,:),Qbar(j,:),U(j,:),[n.ind(2),n.ind(4)],...
                           [n.WBC(j),n.EBC(j)],1);
    end
    % y-direction
    for i = 1:x.n
        Qnew(:,i) = supp.Upwind1D(y,dt,Qstar(:,i),Qbar(:,i),V(:,i),[n.ind(3),n.ind(1)],...
                          [n.SBC(i),n.NBC(i)],2);
    end
end