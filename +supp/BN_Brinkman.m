%% BN with Brinkman flow
%% iteration
while (t < T && tn < num.tnmax)
    if flag.adv == 1 || flag.flow == 1
        if mod(tn,num.flow) == 0 || length(flag.Gb) == 1 || flag.flow == 1
            flag.Gbn = 0;
            % solve for flow
            [U,V,P,nu,K] = supp.Brinkman2D(x,y,G,flow,u,v,p);
            dt = min(flow.cfl./nu,num.dtmax);
            % check outflow
            switch flow.dir
                case 1; outflow = norm(V(end,:),inf);
                case 2; outflow = norm(V(1,:),inf);
                case 3; outflow = norm(U(:,end),inf);
                case 4; outflow = norm(U(1,:),inf);
            end
            if flag.flow == 1
                break;
            end
        end
        if outflow == 0; break; end
        if tn == 0
            % record data
            supp.record;
            % generate initial plots
            supp.figGen;
        end
        tn = tn + 1;
        t = t + dt;
        if t > T; dt = dt - (t - T); t = T; end
        fprintf('tn = %d, t = %g [%s]\n',tn,t,unit.str_time);

        % advection
        Bn = supp.Upwind2D_B(x,y,dt,U,V,B0,b);
        Nn = supp.Upwind2D_N(x,y,dt,U,V,N0,n);
        Ln = L0;                
    else
        if tn == 0
            % record data
            supp.record;
            % generate initial plots
            supp.figGen;
        end
        outflow = 1;
        dt = num.dtmax;
        tn = tn + 1;
        t = t + dt;
        if t > T; dt = dt - (t-T); t = T; end
        fprintf('tn = %d, t = %g [%s]\n',tn,t,unit.str_time);
        Bn = B0; Nn = N0; Ln = L0;
    end
    if outflow == 0; break; end
    if flag.DR == 1
    %% diffusion and reaction
    tau = dt./num.tau;
    for j = 1:num.tau
        [Bnew,Nnew,Lnew] = supp.DR2D(x,y,tau,Bn,Nn,Ln,G,b,n,flag,flow);
        Bn = Bnew;  Nn = Nnew;  Ln = Lnew;
    end
    end
    %% update
    B0 = Bn;    N0 = Nn;    L0 = Ln;
    G.bn = double((B0 >= b.star2 & B0 <= b.star));
    if isempty(flag.Gb) && norm((G.b0-G.bn)./G.b0,inf) > 0.05
        flag.Gb = 1;
    end
    G.b0 = G.bn;
    if mod(tn,num.rec) == 0
        %% record data
        supp.record;
    end
    if mod(tn,num.plot) == 0
        %% save figures
        supp.figGen;
    end
end