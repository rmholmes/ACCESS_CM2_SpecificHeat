% This script calculates RMSE errors for calculating Conservative
% temperature from monthly-averaged potential temperature and
% salinity compared to monthly-averages of the instantaneous
% Conservative Temperature from ACCESS-CM2.

plot_only = 1;
PI_or_his = 1; % 1 = PI-control, 0 = historical simualtion
mname = 'ACCESS_SpecificHeat_PIcontrol_CTptERROR.mat';

%%% This code block generates the post-processed "mname" .mat file
%%% from ACCESS-CM2 output located on NCI (not neccessary if you
%%% are already in possesion of the processed mname file). 
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
Z = ncread(fname,'st_ocean');
[xL,yL] = size(area);
zL = 50;
time = ncread(fname,'time');
tL = length(time);

PS_to_SA = 35.16504/35; % conversion factor Practical Salinity ->
                        % Absolute Salinity

% Initialize variables to time average:

CT_from_pt_SqEr = zeros(xL,yL,zL);
time = [];
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

        % Load variables:
        CT = ncread(fname,'temp',[1 1 1 1],[xL yL zL tL])-273.15;
        PT = ncread(fname,'pot_temp',[1 1 1 1],[xL yL zL tL]);
        PS = ncread(fname,'salt',[1 1 1 1],[xL yL zL tL]);

        CT_from_PT = gsw_CT_from_pt(PS*PS_to_SA,PT);

        ERROR = (CT-CT_from_PT).^2;
        
        CT_from_pt_SqEr = (CT_from_pt_SqEr*sum(DT_A(1:(end-tL))) + sum(ERROR.*repmat(permute(DT_A_t,[4 3 2 1]),[xL yL zL 1]),4))/sum(DT_A);
        
        if (mod(fi,5)==0)
            save(mname,'time','DT_A','lon','lat','area','CT_from_pt_SqEr','Z');
        end
    end
end

%%% This code block plots the figures given the post-processed "mname" .mat file.
%%%
%%%

else
load(mname);

figure;
set(gcf,'defaulttextfontsize',15);
set(gcf,'defaultaxesfontsize',15);
set(gcf,'Position',[232.5000   34.5000  830.0000  802.0000]);

poss = [0.0755    0.71    0.3903    0.2738; ...
        0.5282    0.71    0.3903    0.2738; ...        0.0755    0.41    0.3903    0.2738; ...
        0.5282    0.41    0.3903    0.2738; ...
        0.0755    0.08    0.3903    0.2738; ...
        0.5282    0.08    0.3903    0.2738];

Er = sqrt(CT_from_pt_SqEr);

subplot(2,1,1);
contourf(lon,lat,Er(:,:,1),[0:0.1e-3:3e-3],'linestyle','none');
xlabel('Longitude ($^\circ$E)');
ylabel('Latitude ($^\circ$N)');
caxis([0 3e-3]);
set(gca,'color','k');
text(-278,69,'(a)','Backgroundcolor','w','margin',0.1);
cb = colorbar;
ylabel(cb,'K','Interpreter','latex');
ylim([-80 80]);

[X,Y] = ndgrid(lat(1,:),Z);
subplot(2,1,2);
contourf(X,Y,squeeze(nanmean(Er,1)),[0:0.02e-3:1e-3],'linestyle','none');
xlabel('Latitude ($^\circ$N)');
ylabel('Depth (m)');
caxis([0 1.5e-3]);
xlim([-80 65]);
ylim([0 3000]);
set(gca,'ydir','reverse');
set(gca,'color','k');
text(-79,200,'(b)','Backgroundcolor','w','margin',0.1);
cb = colorbar;
ylabel(cb,'K','Interpreter','latex');
colormap(flipud(pink));

end
