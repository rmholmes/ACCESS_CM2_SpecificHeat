% This script calculates various specific heat related properties
% in ACCESS-CM2 PI control or historical simulations

plot_only = 1;
PI_or_his = 1; % 1 = PI-control, 0 = historical simualtion
mname = 'ACCESS_SpecificHeat_PIcontrol_SWP.mat';
% $$$ PI_or_his = 0; % 1 = PI-control, 0 = historical simualtion
% $$$ mname = 'ACCESS_SpecificHeat_historical.mat';

if (~plot_only)

    if (PI_or_his)
        base = '/g/data/p66/cm2704/archive/bi889/history/ocn/';
        name = 'PIcontrol';
        fname = [base 'ocean_month.nc-08500630'];
    else
        base = '/g/data/p66/cm2704/archive/bj594/history/ocn/';
        name = 'historical';
        fname = [base 'ocean_month.nc-18500630'];
    end

area = ncread(fname,'area_t');
lon = ncread(fname,'geolon_t');
lat = ncread(fname,'geolat_t');
yt_ocean = ncread(fname,'yt_ocean');
[xL,yL] = size(area);
zL = 50;
time = ncread(fname,'time');
tL = length(time);

PT = squeeze(ncread(fname,'pot_temp',[1 1 1 1],[xL yL 1 1]));
mask = ~isnan(PT);

% Some constants
rho0 = 1035;
Cp0 = 3992.10322329649; % J kg-1 degC-1 -> this is the ACCESS-OM2
                       % value. I think the same is used in CM2 but
                       % I'm not 100% sure.
                       % Note: This is for using the pre_teos10 eos
                       % but with conservative temperature as the
                       % temperature variable.
Cp0_teos10 = 3991.86795711963; % The TEOS-10 value (small correction)
Cp0_cor = Cp0/Cp0_teos10; % correction factor

PS_to_SA = 35.16504/35; % conversion factor Practical Salinity ->
                        % Absolute Salinity

% Define salinity bins:
dS = 1;
Smax = 40;
S = 0:dS:Smax;
Sa = dS/2:dS:(Smax-dS/2);
sL = length(Sa);

% Define mean values, mean Cpr and mean dCp/dS:
PT_mean = 18.1405; % mean surface PT over 599 years of PI-control
PS_mean = 34.4005; % mean surface PS over 599 years of PI-control
[~,CpR_on_Cp0] = gsw_CT_first_derivatives(PS_mean*PS_to_SA,PT_mean);
CpR_mean = CpR_on_Cp0*Cp0_teos10; % Cp at mean PT and PS
[~,CT_SA_pt_t,~] = gsw_CT_second_derivatives(PS_mean*PS_to_SA,PT_mean);
dCp_dS_mean = CT_SA_pt_t*Cp0_teos10; % dCp/dS at mean PT and PS

% Initialize variables to time average:
% $$$ Q_S = zeros(sL,1); % Q binned into salinity bins
% $$$ Qf_S = zeros(sL,1); % Qf binned into salinity bins
CpR_on_Cp0 = zeros(xL,yL); % Cp(SA,PT,0dbar)/Cp0 -> ratio of Cp for
                           % potential temperature to Cp0.
Q  = zeros(xL,yL); % Q spatial map
Qf = zeros(xL,yL); % Qf = CpR_on_Cp0*Q = the heat flux the ocean has
                   % recieved if it has interpreted its surface
                   % temperature as potential temperature but still
                   % used Cp0 for its specific heat.
Qf_swp = zeros(xL,yL); % The correction to Qf if short-wave
                       % penetration is correctly taken into
                       % account.
PT_ts = []; % Global surface mean PT
PS_ts = []; % Global surface mean PS
Q_ts = [];  % Globally summed Q time-series
Qf_ts = []; % Globally summed Qf time-series
Qf_swp_ts = []; % Globally summed Qf time-series

Qs = zeros(xL,yL); % Qs = Q*(PS-PS_mean) spatial map;
Qs_ts = []; % Qs_ts = time series of Qs globally integrated.

time = [];
OHC = [];
DT_A = [];

files = dir(base);

for fi = 1:length(files)
    if (strfind(files(fi).name,'month'))

        fname = [base files(fi).name];
        sprintf('Doing %03d of %03d',fi,length(files))
        time_t = ncread(fname,'time');
        DT_A_t = ncread(fname,'average_DT')*86400;
        
        time = cat(1,time,time_t);
        DT_A = cat(1,DT_A,DT_A_t);

        tL = length(time_t);

        % Load surface variables:
% $$$         CT = squeeze(ncread(fname,'temp',[1 1 1 1],[xL yL 1 tL]));
        PT = squeeze(ncread(fname,'pot_temp',[1 1 1 1],[xL yL 1 tL]));
        PS = squeeze(ncread(fname,'salt',[1 1 1 1],[xL yL 1 tL]));

        PT_a = squeeze(nansum(nansum(PT.*repmat(area,[1 1 tL]),1),2)./nansum(nansum(area(mask))));
        PS_a = squeeze(nansum(nansum(PS.*repmat(area,[1 1 tL]),1),2)./nansum(nansum(area(mask))));

        % Load CpR_on_Cp0:
        [~,CpR_on_Cp0_t] = gsw_CT_first_derivatives(PS*PS_to_SA,PT);
        CpR_on_Cp0_t = CpR_on_Cp0_t/Cp0_cor; % Correct for
                                           % differing Cp0s.
        
        Q_t = ncread(fname,'sfc_hflux_from_runoff')+ ...
                 ncread(fname,'sfc_hflux_coupler')+ ...
                 ncread(fname,'sfc_hflux_pme')+ ...
                 squeeze(nansum(ncread(fname,'frazil_3d'),3));
        Qf_t = CpR_on_Cp0_t.*Q_t;
        
        % Shortwave penetration correction:
        PT_3d = ncread(fname,'pot_temp',[1 1 1 1],[xL yL zL tL]);
        PS_3d = ncread(fname,'salt',[1 1 1 1],[xL yL zL tL]);
        [~,CpR_on_Cp0_3d] = gsw_CT_first_derivatives(PS_3d*PS_to_SA,PT_3d);
        CpR_on_Cp0_3d = CpR_on_Cp0_3d/Cp0_cor;
        
        sw_heat = ncread(fname,'sw_heat',[1 1 1 1],[xL yL zL tL]);
        Qf_swp_t = squeeze(nansum(CpR_on_Cp0_3d.*sw_heat,3));

% $$$         % do salinity binning:
% $$$         Q_S_t = zeros(sL,tL);
% $$$         Qf_S_t = zeros(sL,tL);
% $$$         for ti = 1:tL
% $$$             shf = Q_t(:,:,ti).*area;
% $$$             shff = Qf_t(:,:,ti).*area;
% $$$             for si = 1:sL
% $$$                 inds = PS(:,:,ti)>=S(si) & PS(:,:,ti)<S(si+1);
% $$$                 Q_S_t(si,ti) = sum(sum(shf(inds)));
% $$$                 Qf_S_t(si,ti) = sum(sum(shff(inds)));
% $$$             end
% $$$         end
% $$$         Q_S = (Q_S*sum(DT_A(1:(end-tL))) + sum(Q_S_t.*repmat(DT_A_t',[sL 1]),2))/sum(DT_A);
% $$$         Qf_S = (Qf_S*sum(DT_A(1:(end-tL))) + sum(Qf_S_t.*repmat(DT_A_t',[sL 1]),2))/sum(DT_A);

        Q = (Q*sum(DT_A(1:(end-tL))) + sum(Q_t.*repmat(permute(DT_A_t,[3 2 1]),[xL yL 1]),3))/sum(DT_A);
        CpR_on_Cp0 = (CpR_on_Cp0*sum(DT_A(1:(end-tL))) + sum(CpR_on_Cp0_t.*repmat(permute(DT_A_t,[3 2 1]),[xL yL 1]),3))/sum(DT_A);
        Qf = (Qf*sum(DT_A(1:(end-tL))) + sum(Qf_t.*repmat(permute(DT_A_t,[3 2 1]),[xL yL 1]),3))/sum(DT_A);
        Qf_swp = (Qf_swp*sum(DT_A(1:(end-tL))) + sum(Qf_swp_t.*repmat(permute(DT_A_t,[3 2 1]),[xL yL 1]),3))/sum(DT_A);

        Q_ts = cat(1,Q_ts,squeeze(nansum(nansum(Q_t.*repmat(area,[1 1 tL]),1),2)));
        Qf_ts = cat(1,Qf_ts,squeeze(nansum(nansum(Qf_t.*repmat(area,[1 1 tL]),1),2)));
        Qf_swp_ts = cat(1,Qf_swp_ts,squeeze(nansum(nansum(Qf_swp_t.*repmat(area,[1 1 tL]),1),2)));
        PT_ts = cat(1,PT_ts,PT_a);
        PS_ts = cat(1,PS_ts,PS_a);

        Qs_t = Q_t.*(PS-PS_mean);
        Qs = (Qs*sum(DT_A(1:(end-tL))) + sum(Qs_t.*repmat(permute(DT_A_t,[3 2 1]),[xL yL 1]),3))/sum(DT_A);
        Qs_ts = cat(1,Qs_ts,squeeze(nansum(nansum(Qs_t.*repmat(area,[1 1 tL]),1),2)));

        OHC_t = squeeze(Cp0*nansum(nansum(nansum(ncread(fname,'temp_rhodzt',[1 1 1 1],[xL yL zL ...
                            tL]),3).*area,1),2));
        OHC = cat(1,OHC,OHC_t);
        
        if (mod(fi,5)==0)
            save(mname,'OHC','time','DT_A','mask','Cp0','Cp0_teos10','CpR_mean', 'dCp_dS_mean','lon','lat','area','yt_ocean', ...
                 'CpR_on_Cp0','Q','Qf','Q_ts','Qs','Qf_ts','Qs_ts','S','Sa','Qf_swp','Qf_swp_ts', ...
                 'PT_ts','PS_ts');%, ...
% $$$                  ,'Q_S','Qf_S');
        end
    end
end

else
load(mname);

Qf_swp(Qf_swp==0) = NaN;
Qf = Qf+Qf_swp; % Correct for SWP.

dQ = Qf-Q;
Qs = dCp_dS_mean*Qs/Cp0;

figure;
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
set(gcf,'Position',[480.3    12.3    1537.3    1348.7]);

poss = [0.0755    0.71    0.3903    0.2738; ...
        0.5282    0.71    0.3903    0.2738; ...
        0.0755    0.41    0.3903    0.2738; ...
        0.5282    0.41    0.3903    0.2738; ...
        0.0755    0.08    0.3903    0.2738; ...
        0.5282    0.08    0.3903    0.2738];
        

% Plotting (spatial maps, 1-4):

rg = 0.02;
subplot(3,2,1);
contourf(lon,lat,CpR_on_Cp0,[0 1-rg:rg/40:1+rg 10],'linestyle','none');
cb = colorbar;
set(gca,'color',0.5*[1 1 1]);
caxis([1-rg 1+rg]);
ylim([-80 75]);
text(-275,67,'(a) $C_p(S_A,\theta,0)/C_p^0$','Backgroundcolor','w');
set(gca,'xticklabel',[]);
ylabel('Latitude ($^\circ$N)');
set(gca,'Position',poss(1,:));

subplot(3,2,2);
contourf(lon,lat,Q,[-1000 -200:5:200 1000],'linestyle','none');
cb = colorbar;
ylabel(cb,'Wm$^{-2}$','Interpreter','latex');
ylim([-80 75]);
set(gca,'color',0.5*[1 1 1]);
caxis([-200 200]);
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);
set(gca,'Position',poss(2,:));
text(-275,67,'(b) $Q$','Backgroundcolor','w');

subplot(3,2,3);
contourf(lon,lat,dQ,[-10 -0.4:0.01:0.4 10],'linestyle','none');
cb = colorbar;
hold on;
[c,h] = contour(lon,lat,dQ,[0.1 0.2 0.3],'-k');
clabel(c,h);
[c,h] = contour(lon,lat,dQ,[-0.1 -0.2 -0.3],'--k');
clabel(c,h);
% $$$ ylabel(cb,'Wm$^{-2}$');
ylim([-80 75]);
caxis([-0.4 0.4]);
set(gca,'color',0.5*[1 1 1]);
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
set(gca,'Position',poss(3,:));
text(-275,67,'(c) $\Delta Q$','Backgroundcolor','w');

Qvar = Qs;
lab = 'Q_S';
clev = [-0.4:0.01:0.4];
% $$$ Qvar = Qf_swp;
% $$$ lab = 'Q_{SWP}';
% $$$ clev = [-0.05:0.001:0.05];
subplot(3,2,4);
contourf(lon,lat,Qvar,[-10 clev 10],'linestyle','none');
cb = colorbar;
hold on;
[c,h] = contour(lon,lat,Qvar,[0.1 0.2 0.3],'-k');
clabel(c,h);
[c,h] = contour(lon,lat,Qvar,[-0.1 -0.2 -0.3],'--k');
clabel(c,h);
ylabel(cb,'Wm$^{-2}$','Interpreter','latex');
ylim([-80 75]);
caxis([clev(1) clev(end)]);
xlabel('Longitude ($^\circ$E)');
set(gca,'yticklabel',[]);
set(gca,'Position',poss(4,:));
set(gca,'color',0.5*[1 1 1]);
text(-275,67,['(d) $\Delta ' lab '$'],'Backgroundcolor','w');

colormap(cmocean('balance'));

% Histograms:
% Area-weighting:
ntot = length(find(mask));
Aavg = nansum(nansum(area.*mask))/ntot;

dQ_A = dQ.*area;
dQ_A = dQ_A(mask)/Aavg;
Qs_A = Qvar.*area;
Qs_A = Qs_A(mask)/Aavg;

ymax = 0.08;
Binedges = [-100 -10:0.005:10 100];
% $$$ Binedges = [-100 -10:0.001:10 100];

% $$$ % histogram check:
% $$$ subplot(2,2,1);
% $$$ pcolPlot(lon,lat,area/Aavg);
% $$$ title('Area-weighting factor');
% $$$ subplot(2,2,2);
% $$$ pcolPlot(lon,lat,dQ.*area/Aavg);
% $$$ title('Area-weighted dQ');
% $$$ caxis([-1 1]);
% $$$ colormap('redblue');
% $$$ set(gca,'color','k');
% $$$ subplot(2,2,3);
% $$$ var = dQ_A;
% $$$ h = histogram(var,'Normalization','probability','BinEdges',Binedges,'DisplayStyle', 'stairs','EdgeColor','b');%Limits',[-0.7 0.7]);
% $$$ hold on;
% $$$ var = dQ(~isnan(dQ));
% $$$ h = histogram(var,'Normalization','probability','BinEdges',Binedges,'DisplayStyle', 'stairs','EdgeColor','r');%Limits',[-0.7 0.7]);
% $$$ xlim([-0.8 0.8]);
% $$$ legend('Area-weighted dQ','Not Area-Weighted dQ');
% $$$ xlabel('Wm-2');
% $$$ ylabel('Normalized Count');
% $$$ 
% $$$ subplot(2,2,4);
% $$$ var = dQ_A;
% $$$ h = histogram(var,'Normalization','probability','BinEdges',Binedges,'DisplayStyle', 'stairs','EdgeColor','b');%Limits',[-0.7 0.7]);
% $$$ hold on;
% $$$ var = dQ(~isnan(dQ));
% $$$ h = histogram(var,'Normalization','probability','BinEdges',Binedges,'DisplayStyle', 'stairs','EdgeColor','r');%Limits',[-0.7 0.7]);
% $$$ xlim([-0.8 0.8]);
% $$$ ylim([0 0.002]);

subplot(3,2,5);
var = dQ_A;
histogram(var,'Normalization','probability','BinEdges',Binedges);%Limits',[-0.4 0.4]);
hold on;
plot(mean(var)*[1 1],[ymax*0.43 ymax*0.57],'-r','linewidth',3);
text(mean(var),ymax*0.6,sprintf('%3.3f',mean(var)),'HorizontalAlignment','left','color','r');%,'BackgroundColor','w','Margin',0.04);
top = prctile(var,95);
bot = prctile(var,5);
XSXSXplot(top*[1 1],[ymax*0.43 ymax*0.57],'-r','linewidth',3);
text(top,ymax*0.4,sprintf('%3.3f',top),'HorizontalAlignment','left','color','r');%,'BackgroundColor','w','Margin',0.04);
plot(bot*[1 1],[ymax*0.43 ymax*0.57],'-r','linewidth',3);
text(bot,ymax*0.4,sprintf('%3.3f',bot),'HorizontalAlignment','right','color','r');%,'BackgroundColor','w','Margin',0.04);
plot([bot top],[ymax*0.5 ymax*0.5],'-r','linewidth',3);
% $$$ text(0.21,ymax*0.92,['$\mu =$ ' sprintf('%3.3f',nanmean(var)) 'Wm$^{-2}$']);
% $$$ text(0.21,ymax*0.82,['$\sigma = $ ' sprintf('%3.3f',nanstd(var)) 'Wm$^{-2}$']);
text(-0.39,ymax*0.95,'(e) $\Delta Q$ histogram');
xlabel('Cell area-weighted $\Delta Q$ (Wm$^{-2}$)');
xlim([-0.4 0.4]);
ylim([0 ymax]);
grid on;
set(gca,'Position',poss(5,:));
ylabel('Normalized count');

% $$$ % log scale:
% $$$ set(gca,'yscale','log');
% $$$ ylim([1e-5 0.3]);

subplot(3,2,6);
var = Qs_A;
histogram(var,'Normalization','probability','BinEdges',Binedges);%Limits',[-0.4 0.4]);%,'BinLimits',[-0.4 0.4]);
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
text(-0.39,ymax*0.95,['(f) $\Delta ' lab '$ histogram']);
xlabel(['Cell area-weighted $\Delta ' lab '$ (Wm$^{-2}$)']);
xlim([clev(1) clev(end)]);
ylim([0 ymax]);
grid on;
set(gca,'Position',poss(6,:));

% Latitudinal plot of dQ:
figure;

dyt_ocean = diff(yt_ocean);
dyt_ocean = cat(1,dyt_ocean(1),(dyt_ocean(1:end-1)+dyt_ocean(2:end))/2,dyt_ocean(end));
plot(yt_ocean,filter_field(nansum(dQ.*area,1)'./dyt_ocean,3,'-t'),'-k','linewidth',2);
xlabel('Latitude ($^\circ$N)');
ylabel('$\Delta Q$ (W / $\circ$ latitude)');
xlim([-80 80]);
grid on;



% $$$ % log scale:
% $$$ set(gca,'yscale','log');
% $$$ ylim([1e-5 0.3]);


% $$$ %% Time series:
% $$$ figure;
% $$$ set(gcf,'Position',[480.3    303.0    1537.3    879.3]);
% $$$ 
% $$$ time = time/365.25;
% $$$ tL = length(time);
% $$$ OHC = OHC-OHC(1);
% $$$ 
% $$$ % time-integrate others:
% $$$ Q = cumsum(Q_ts.*DT_A);
% $$$ Qf = cumsum(Qf_ts.*DT_A);
% $$$ Qs = cumsum(Qs_ts.*DT_A);
% $$$ 
% $$$ % Annual average:
% $$$ tL_a = floor(tL/12);
% $$$ OHC_a = zeros(tL_a,1);
% $$$ Q_a = OHC_a;
% $$$ Qf_a = OHC_a;
% $$$ Qs_a = OHC_a;
% $$$ for ti=1:tL_a
% $$$     OHC_a(ti) = mean(OHC(((ti-1)*12+1):(ti*12)));
% $$$     Q_a(ti) = mean(Q(((ti-1)*12+1):(ti*12)));
% $$$     Qf_a(ti) = mean(Qf(((ti-1)*12+1):(ti*12)));
% $$$     Qs_a(ti) = mean(Qs(((ti-1)*12+1):(ti*12)));
% $$$ end
% $$$ if (PI_or_his)
% $$$     time_a = 1:tL_a;
% $$$ else
% $$$     time_a = 1850:(1850+tL_a-1);
% $$$ end
% $$$ 
% $$$ % zero initial values again:
% $$$ OHC_a = OHC_a-OHC_a(1);
% $$$ Q_a = Q_a-Q_a(1);
% $$$ Qf_a = Qf_a-Qf_a(1);
% $$$ Qs_a = Qs_a-Qs_a(1);
% $$$ 
% $$$ subplot(2,1,1);
% $$$ plot(time_a,OHC_a,'linewidth',2);
% $$$ hold on;
% $$$ plot(time_a,Q_a,'--r','linewidth',2);
% $$$ plot(time_a,Qf_a,'-b','linewidth',2);
% $$$ xlabel('Year');
% $$$ xlim([time_a(1) time_a(end)]);
% $$$ ylabel('Heat Content Anomaly [J]');
% $$$ legend('Total OHC',...
% $$$        'Area and time-integrated $Q$',...
% $$$        'Area and time-integrated $Q^f$');
% $$$ grid on;
% $$$ 
% $$$ subplot(2,1,2);
% $$$ plot(time_a,Qf_a-Q_a,'-k','linewidth',2);
% $$$ hold on;
% $$$ plot(time_a,Qs_a,'-b','linewidth',2);
% $$$ plot(time_a,Qf_a-Q_a-Qs_a,'-r','linewidth',2);
% $$$ xlabel('Year');
% $$$ xlim([time_a(1) time_a(end)]);
% $$$ ylabel('Heat Content Anomaly [J]');
% $$$ legend('Area and time-integrated $\Delta Q$',...
% $$$        'Area and time-integrated $Q_S$',...
% $$$        'Area and time-integrated $Q_T=\Delta Q-Q_S$');
% $$$ grid on;
% $$$ 
% $$$ %% Heat flux binned in salinity coordinates:
% $$$ PW = 1e15;
% $$$ 
% $$$ figure;
% $$$ subplot(2,1,1);
% $$$ plot(Sa,Q_S/PW);
% $$$ ylabel('Heat flux (PW) into each 0.5psu salinity band');
% $$$ xlabel('Practical Salinity (psu)');
% $$$ 
% $$$ subplot(2,1,2);
% $$$ plot(S,cat(1,0,cumsum(Q_S)/1e15));
% $$$ ylabel('Heat flux (PW) into sea water with salinity less than S');
% $$$ xlabel('Practical Salinity S (psu)');
% $$$ % $$$ 
% $$$ % $$$ 
% $$$ % $$$ total_nowgt = nansum(nansum(Q.*area));
% $$$ % $$$ total_wgt = nansum(nansum(Qf.*area));
% $$$ % $$$ 
% $$$ % $$$ sprintf('%3.4fPW',total_nowgt/1e15)
% $$$ % $$$ sprintf('%3.4fPW',total_wgt/1e15)
end
