% This script calculates various specific heat related properties
% in ACCESS-CM2 PI control or historical simulations

% addpath(genpath('~/software/matlab-utilities/'));
% startup;

plot_only = 1;
PI_or_his = 1; % 1 = PI-control, 0 = historical simualtion
mname = 'ACCESS_SpecificHeat_PIcontrol_SWP.mat';
mname = 'ACCESS_SpecificHeat_PIcontrol_SWP_2Cwarm.mat';

load(mname);

Qf_swp(Qf_swp==0) = NaN;
Qf_swp_warm(Qf_swp_warm==0) = NaN;
Qf = Qf+Qf_swp; % Correct for SWP.
Qf_warm = Qf_warm+Qf_swp_warm; % Correct for SWP.

dQ = Qf-Q;
dQ_warm = Qf_warm-Q;
Qs = dCp_dS_mean*Qs/Cp0;

figure;
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
set(gca,'Position',[0.4108 0.11 0.2134 0.3412]);

rg = 0.02;

subplot(2,3,1);
contourf(lon,lat,dQ,[-10 -0.4:0.01:0.4 10],'linestyle','none');
cb = colorbar;
hold on;
[c,h] = contour(lon,lat,dQ,[0.1 0.2 0.3],'-k');
clabel(c,h);
[c,h] = contour(lon,lat,dQ,[-0.1 -0.2 -0.3],'--k');
clabel(c,h);
ylabel(cb,'Wm$^{-2}$','Interpreter','latex');
ylim([-80 75]);
caxis([-0.4 0.4]);
set(gca,'color',0.5*[1 1 1]);
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
text(-275,67,'(a) $\Delta Q$','Backgroundcolor','w');

subplot(2,3,2);
contourf(lon,lat,dQ_warm,[-10 -0.4:0.01:0.4 10],'linestyle','none');
cb = colorbar;
hold on;
[c,h] = contour(lon,lat,dQ_warm,[0.1 0.2 0.3],'-k');
clabel(c,h);
[c,h] = contour(lon,lat,dQ_warm,[-0.1 -0.2 -0.3],'--k');
clabel(c,h);
ylabel(cb,'Wm$^{-2}$','Interpreter','latex');
ylim([-80 75]);
caxis([-0.4 0.4]);
set(gca,'color',0.5*[1 1 1]);
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
text(-275,67,'(b) $\Delta Q$ (with $2^\circ$ warming)','Backgroundcolor','w');

subplot(2,3,3);
contourf(lon,lat,dQ_warm-dQ,[-10 -0.05:0.001:0.05 10],'linestyle','none');
cb = colorbar;
hold on;
[c,h] = contour(lon,lat,dQ_warm-dQ,[0.025 0.05 0.075],'-k');
clabel(c,h);
[c,h] = contour(lon,lat,dQ_warm-dQ,[-0.025 -0.05 -0.075],'--k');
clabel(c,h);
hold on;
ylabel(cb,'Wm$^{-2}$','Interpreter','latex');
ylim([-80 75]);
caxis([-0.05 0.05]);
set(gca,'color',0.5*[1 1 1]);
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
text(-275,67,'(c) Difference','Backgroundcolor','w');

colormap(redblue);

% Histograms:
% Area-weighting:
ntot = length(find(mask));
Aavg = nansum(nansum(area.*mask))/ntot;

dQ_A = dQ.*area;
dQ_A = dQ_A(mask)/Aavg;
dQ_warm_A = dQ_warm.*area;
dQ_warm_A = dQ_warm_A(mask)/Aavg;

ymax = 0.08;
Binedges = [-100 -10:0.005:10 100];
% $$$ Binedges = [-100 -10:0.001:10 100];

subplot(2,3,4);
var = dQ_A;
histogram(var,'Normalization','probability','BinEdges',Binedges);%Limits',[-0.4 0.4]);
hold on;
plot(mean(var)*[1 1],[ymax*0.43 ymax*0.57],'-r','linewidth',3);
text(mean(var),ymax*0.6,sprintf('%3.3f',mean(var)),'HorizontalAlignment','left','color','r');%,'BackgroundColor','w','Margin',0.04);
top = prctile(var,95);
bot = prctile(var,5);
plot(top*[1 1],[ymax*0.43 ymax*0.57],'-r','linewidth',3);
text(top,ymax*0.4,sprintf('%3.3f',top),'HorizontalAlignment','left','color','r');%,'BackgroundColor','w','Margin',0.04);
plot(bot*[1 1],[ymax*0.43 ymax*0.57],'-r','linewidth',3);
text(bot,ymax*0.4,sprintf('%3.3f',bot),'HorizontalAlignment','right','color','r');%,'BackgroundColor','w','Margin',0.04);
plot([bot top],[ymax*0.5 ymax*0.5],'-r','linewidth',3);
% $$$ text(0.21,ymax*0.92,['$\mu =$ ' sprintf('%3.3f',nanmean(var)) 'Wm$^{-2}$']);
% $$$ text(0.21,ymax*0.82,['$\sigma = $ ' sprintf('%3.3f',nanstd(var)) 'Wm$^{-2}$']);
text(-0.39,ymax*0.95,'(d) $\Delta Q$ histogram');
xlabel('Cell area-weighted $\Delta Q$ (Wm$^{-2}$)');
xlim([-0.4 0.4]);
ylim([0 ymax]);
grid on;
ylabel('Normalized count');

subplot(2,3,5);
var = dQ_warm_A;
histogram(var,'Normalization','probability','BinEdges',Binedges);%Limits',[-0.4 0.4]);
hold on;
plot(mean(var)*[1 1],[ymax*0.43 ymax*0.57],'-r','linewidth',3);
text(mean(var),ymax*0.6,sprintf('%3.3f',mean(var)),'HorizontalAlignment','left','color','r');%,'BackgroundColor','w','Margin',0.04);
top = prctile(var,95);
bot = prctile(var,5);
plot(top*[1 1],[ymax*0.43 ymax*0.57],'-r','linewidth',3);
text(top,ymax*0.4,sprintf('%3.3f',top),'HorizontalAlignment','left','color','r');%,'BackgroundColor','w','Margin',0.04);
plot(bot*[1 1],[ymax*0.43 ymax*0.57],'-r','linewidth',3);
text(bot,ymax*0.4,sprintf('%3.3f',bot),'HorizontalAlignment','right','color','r');%,'BackgroundColor','w','Margin',0.04);
plot([bot top],[ymax*0.5 ymax*0.5],'-r','linewidth',3);
% $$$ text(0.21,ymax*0.92,['$\mu =$ ' sprintf('%3.3f',nanmean(var)) 'Wm$^{-2}$']);
% $$$ text(0.21,ymax*0.82,['$\sigma = $ ' sprintf('%3.3f',nanstd(var)) 'Wm$^{-2}$']);
text(-0.39,ymax*0.95,'(e) $\Delta Q$ histogram, with $2^\circ$ warming');
xlabel('Cell area-weighted $\Delta Q$ (Wm$^{-2}$)');
xlim([-0.4 0.4]);
ylim([0 ymax]);
grid on;
ylabel('Normalized count');

ymax = 0.08;
Binedges = [-100 -5:0.001:5 100];

subplot(2,3,6);
var = dQ_warm_A - dQ_A;
histogram(var,'Normalization','probability','BinEdges',Binedges);%Limits',[-0.4 0.4]);
hold on;
plot(mean(var)*[1 1],[ymax*0.43 ymax*0.57],'-r','linewidth',3);
text(mean(var),ymax*0.6,sprintf('%3.3f',mean(var)),'HorizontalAlignment','left','color','r');%,'BackgroundColor','w','Margin',0.04);
top = prctile(var,95);
bot = prctile(var,5);
plot(top*[1 1],[ymax*0.43 ymax*0.57],'-r','linewidth',3);
text(top,ymax*0.4,sprintf('%3.3f',top),'HorizontalAlignment','left','color','r');%,'BackgroundColor','w','Margin',0.04);
plot(bot*[1 1],[ymax*0.43 ymax*0.57],'-r','linewidth',3);
text(bot,ymax*0.4,sprintf('%3.3f',bot),'HorizontalAlignment','right','color','r');%,'BackgroundColor','w','Margin',0.04);
plot([bot top],[ymax*0.5 ymax*0.5],'-r','linewidth',3);
% $$$ text(0.21,ymax*0.92,['$\mu =$ ' sprintf('%3.3f',nanmean(var)) 'Wm$^{-2}$']);
% $$$ text(0.21,ymax*0.82,['$\sigma = $ ' sprintf('%3.3f',nanstd(var)) 'Wm$^{-2}$']);
text(-0.068,ymax*0.95,'(f) Difference histogram');
xlabel('Cell area-weighted $\Delta Q$ (Wm$^{-2}$)');
xlim([-0.07 0.07]);
ylim([0 ymax]);
grid on;
ylabel('Normalized count');
