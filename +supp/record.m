%% record data
dir.tn = sprintf('%s/plt%d',dir.version,tn);
if not(isfolder(dir.tn)); mkdir(dir.tn); end

X = full(x.p);
Y = full(y.p);
G0 = full(G.omega');
Bt = full(B0);
Nt = full(N0);
Lt = full(L0);
save(sprintf('%s/T.dat',dir.tn),'t','-ascii');
save(sprintf('%s/X.dat',dir.tn),'X','-ascii')
save(sprintf('%s/Y.dat',dir.tn),'Y','-ascii')
save(sprintf('%s/G.dat',dir.tn),'G0','-ascii')    
save(sprintf('%s/B.dat',dir.tn),'Bt','-ascii')
save(sprintf('%s/N.dat',dir.tn),'Nt','-ascii')
save(sprintf('%s/L.dat',dir.tn),'Lt','-ascii') 

if flag.adv == 1 || flag.flow == 1
    Ut = full(U); 
    Vt = full(V); 
    Pt = full(P);
    Kt = full(K);
    save(sprintf('%s/U.dat',dir.tn),'Ut','-ascii')
    save(sprintf('%s/V.dat',dir.tn),'Vt','-ascii')
    save(sprintf('%s/P.dat',dir.tn),'Pt','-ascii')
    save(sprintf('%s/K.dat',dir.tn),'Kt','-ascii')  
end