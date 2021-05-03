%% function upwind1D
function Qnew = Upwind1D(x,dt,Qold,Qbar,velocity,BC_ind,BC_val,dir)
    ds = x.h;
    Q0 = Qold(:);
    Q1 = Qbar(:);
    vel = velocity(:);
    
    LBC = BC_ind(1).*(BC_val(1).*max(0,BC_ind(1))+Q1(1).*min(0,BC_ind(1)));
    RBC = BC_ind(2).*(BC_val(2).*max(0,BC_ind(2))+Q1(end).*min(0,BC_ind(2)));
    
    Qtemp = [LBC;Q1;RBC];
    F = 0.5.*((vel-abs(vel)).*Qtemp(2:end)+(vel+abs(vel)).*Qtemp(1:end-1));
    
    Q = Q0 - dt./ds.*(F(2:end)-F(1:end-1));
    switch dir
        case 1 % x-direction
            Qnew = Q';
        case 2
            Qnew = Q;
    end
end
