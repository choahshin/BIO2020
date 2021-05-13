% figure generation
wx = 8; wy = 7;
if t == 0
    %% generate initial domain image
    plt1 = figure('visible','off');
    h1 = pcolor(x.p,y.p,G.omega');
    h1.EdgeColor = 'none';
    hold on
    colormap(plt1,[0 0 0; 255,255,255;150 150 150]./255);
    caxis([0,2])
    hc1=colorbar;
    hc1.XTickLabel = {'\Omega_r','\Omega_{w}','\Omega_b'};
    hc1.XTick = [1/3,1,5/3];
    hc1.FontSize = 14;
    hold off
    set(gca,'fontsize',14);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 wx wy]);
    print(sprintf('%s/InitialOmega',dir.tn),'-dpng','-r100');
    cla(plt1)
end
if flag.adv == 1 || flag.flow == 1
    Uave = (U(:,1:end-1)+U(:,2:end))./2;
    Vave = (V(1:end-1,:)+V(2:end,:))./2;
    velmag = sqrt(Uave.^2 + Vave.^2);
    velmax = full(max(max(velmag(G.flow_id))));
    velbar = velmag'./velmax;
    velbar(G.rock_id) = nan;
    plt2 = figure('visible','off');
    pcolor(x.p,y.p,velbar');
    shading interp
    colormap(plt2,jet)
    colorbar
    caxis([0,1]);
    hold on
    stl = streamslice(x.p,y.p,Uave,Vave,.4);
    set(stl,'Color',cc.gray,'linewidth',1.5);
    title(sprintf('$|u|$ at tn = %d, t = %.4e [%s]\n $|u|_{\\infty}$ = %.4e %s/%s\n',...
          tn,t,unit.str_time,velmax,unit.str_len,unit.str_time),'interpreter','latex');
    hold off
    set(gca,'fontsize',14);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 wx wy]);
    print(sprintf('%s/u_tn%d',dir.tn,tn),'-dpng','-r100');
    clear Uave Vave velmag velmax velbar;
    cla(plt2);
    %%      
    plt3 = figure('visible','off');
    Ptemp = P'; Ptemp(G.rock_id) = nan;
    P_max0 = full(max(max(Ptemp(G.flow_id))));
    Ptemp = Ptemp./P_max0;
    surf(x.p, y.p, Ptemp');
    colormap(plt3, jet)
    shading interp
    colorbar
    xlabel('x'); ylabel('y'); zlabel('p');
    if flow.dir <= 2; view(120,50); else; view(15,55); end
    xlim([x.c(1),x.c(end)]);
    ylim([y.c(1),y.c(end)]);
    title(sprintf('P at tn = %d, t = %.4e [%s]\n $|p|_{\\infty}$ = %.4e [Pa]\n',...
          tn,t,unit.str_time,P_max0.*unit.length.*unit.time./unit.mass),'interpreter','latex');
    set(gca,'fontsize',14);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 wx wy]);
    print(sprintf('%s/p_tn%d',dir.tn,tn),'-dpng','-r100');
    clear Ptemp P_max0;
    cla(plt3)
end
%% biomass
plt4 = figure('visible','off');
Bbar = B0'; 
Bbar(G.rock_id) = nan;
pcolor(x.p,y.p,Bbar');
shading interp
colormap(plt4,jet)
colorbar
caxis([0,1]);
title(sprintf(['Biomass, tn = %d, t = %.4e [%s]\n',...
               '$B^*$ = %g, $B_*$ = %g, $B_{max}$ = %g, $B_{min}$ = %g\n'],...
              tn, t,unit.str_time, b.star,b.star2, ...
              full(max(max(Bbar(G.flow_id)))),...
              full(min(min(Bbar(G.flow_id))))),...
      'fontsize',12,'interpreter','latex');
set(gca,'fontsize',14);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 wx wy]);
print(sprintf('%s/B_tn%d',dir.tn,tn),'-dpng','-r100');   
clear Bbar;
cla(plt4);
%% nutrient
plt5 = figure('visible','off');
Nbar = N0'; 
Nbar(G.rock_id) = nan;
pcolor(x.p,y.p,Nbar');
shading interp
colormap(plt5,jet)
colorbar
caxis([0,1]);
title(sprintf(['Nutrient, tn = %d, t = %.4e [%s]\n',...
               '$N_{in}$ = %g, $N_o$, = %g, $N_{max}$ = %g, $N_{min}$ = %g\n'],...
              tn, t, unit.str_time, n.inlet./n.inlet,n.init./n.inlet, ...
              full(max(max(Nbar(G.flow_id)))),...
              full(min(min(Nbar(G.flow_id))))),...
      'fontsize',12,'interpreter','latex');
set(gca,'fontsize',14);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 wx wy]);
print(sprintf('%s/N_tn%d',dir.tn,tn),'-dpng','-r100');          
clear Nbar;
cla(plt5);