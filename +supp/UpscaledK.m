%% Script to get upscaled permeability
if G.kb > 0
    [U,V,P,nu,K] = supp.Brinkman2D(x,y,G,flow,u,v,p);
else
    [U,V,P,nu,K] = supp.Stokes2D(x,y,G,flow,u,v,p);
end
supp.record;
supp.figGen;
if ne == 2    
    nf = ne;
    if nd <= 2 
        [Velu2,Velv2,Px2,Py2] = supp.VolAve(x,y,U,V,P);
        nd = 3; 
    else
        [Velu1,Velv1,Px1,Py1] = supp.VolAve(x,y,U,V,P);
        nd = 1; 
    end
    supp.input_data;
    if G.kb > 0
        [U,V,P,nu,K] = supp.Brinkman2D(x,y,G,flow,u,v,p);
    else
        [U,V,P,nu,K] = supp.Stokes2D(x,y,G,flow,u,v,p);
    end
    supp.record;
    supp.figGen;
    if nd <= 2
        [Velu2,Velv2,Px2,Py2] = supp.VolAve(x,y,U,V,P);
        nd = 3; 
    else
        [Velu1,Velv1,Px1,Py1] = supp.VolAve(x,y,U,V,P);
        nd = 1; 
    end
    Grad = [Px1,Py1,0,0;
            0,0,Px1,Py1;
            Px2,Py2,0,0;
            0,0,Px2,Py2];
    Vm = [Velu1;Velv1;Velu2;Velv2];
    
    K = Grad\Vm.*flow.mu;
    K = full(reshape(K',2,2));
    save(sprintf('%s/K_tensor.dat',dir.tn),'K','-ascii');
end
    