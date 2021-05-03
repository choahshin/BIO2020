%% record data
dir.tn = sprintf('%s/plt%d',dir.version,tn);
if not(isfolder(dir.tn)); mkdir(dir.tn); end

X = full(x.p);
Y = full(y.p);
G0 = full(G.omega');
Bt = full(B0);
Nt = full(N0);
Lt = full(L0);
save(sprintf('%s/T.txt',dir.tn),'t','-ascii');
save(sprintf('%s/X.txt',dir.tn),'X','-ascii')
save(sprintf('%s/Y.txt',dir.tn),'Y','-ascii')
save(sprintf('%s/G.txt',dir.tn),'G0','-ascii')    
save(sprintf('%s/B.txt',dir.tn),'Bt','-ascii')
save(sprintf('%s/N.txt',dir.tn),'Nt','-ascii')
save(sprintf('%s/L.txt',dir.tn),'Lt','-ascii') 

if flag.adv == 1
    Ut = full(U); 
    Vt = full(V); 
    Pt = full(P);
    Kt = full(K);
    save(sprintf('%s/U.txt',dir.tn),'Ut','-ascii')
    save(sprintf('%s/V.txt',dir.tn),'Vt','-ascii')
    save(sprintf('%s/P.txt',dir.tn),'Pt','-ascii')
    save(sprintf('%s/K.txt',dir.tn),'Kt','-ascii')  
end