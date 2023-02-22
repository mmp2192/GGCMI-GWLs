%		    mmp_agmip305
%
%  This script creates maps showing the shape of the
%  production change curve for each global warming
%  level using GGCMI Phase 3 simulations. This version
%  uses pre-processed GGCMI Phase 3 outputs. This is an
%  updated version of mmp_agmip205.m.
%
%  -- Creates a line plot showing global mean temperature
%  anomalies for all ISIMIP GCMs and scenarios.
%
%  -- Creates a scatter plot showing global warming
%  level (x-axis) against CO2 level for each GCM and
%  scenario.
%
%  -- Creates a scatter plot showing global warming
%  level (x-axis) against the middle year of the crossing
%  period for that GWL for each GCM and scenario.
%
%  -- Creates a map of global mean production anomaly
%  for all crops across all ISIMIP GCMs and scenarios as
%  the median across all GGCMs.
%
%  -- Creates a categorical map showing the shape of
%  the median GGCMI production trajectory as GWLs increase,
%  for each crop and scenario. Categories can be one of the
%  following:
%      Linear Increase: production peaks at highest GWL and
%          is lowest at lowest GWL
%      Linear Decrease: production peaks at lowest GWL and
%          is highest at lowest GWL
%      Peak Increase: production peaks between lowest and
%          highest GWLs, with production values at highest
%          GWL above production values at lowest GWL
%      Peak Decrease: production peaks between lowest and
%          highest GWLs, with production values at highest
%          GWL below production values at lowest GWL
%      Dip Increase: lowest production values are between
%          lowest and highest GWLs, with production values
%          at highest GWL above production values at lowest GWL
%      Dip Decrease: lowest production values are between
%          lowest and highest GWLs, with production values
%          at highest GWL below production values at lowest GWL
%      Linear Neutral: production at highest GWL is within
%          range of production values at lowest GWL with no
%          significant changes between the two
%      Peak Neutral: production peaks between lowest and
%          highest GWLs, with production values at highest
%          GWL within range of production values at lowest GWL
%      Dip Neutral: lowest production values are between
%          lowest and highest GWLs, with production values
%          at highest GWL within range of production values
%          at lowest GWL
%  Note that "within range" means we use a window of 3% production
%  anomaly before categorizing production values as significantly
%  above or below starting/ending points.
%
%  INPUTS: none
%  OUTPUTS: none
%
%
%  This is written to be run on Discover.
%
%
%		                             author: Meridel Phillips
%                                            mmp2192@columbia.edu
%				             date:  8/17/2022
%
function mmp_agmip305();
%--------------------------------------------------
%--------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Begin Debug
%% End Debug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modeltypes = {'acea','crover','cygma1p74','dssat-pythia','epic-iiasa','isam','ldndc','lpjml','pdssat','pepic','promet','simplace-lintul5'};
modeltypesnice = {'ACEA','CROVER','CYGMA','DSSAT-Pythia','EPIC-IIASA','ISAM','LandscapeDNDC','LPJmL','pDSSAT','PEPIC','PROMET','SIMPLACE-LINTUL5'};
scentypes = {'ssp126','ssp370','ssp585'};
scentypesnice = {'SSP126','SSP370','SSP585'};
gcmtypes = {'gfdl-esm4','ipsl-cm6a-lr','mpi-esm1-2-hr','mri-esm2-0','ukesm1-0-ll'};
gcmtypesnice = {'GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL'};
croptypes = {'mai','whe','soy','ric'};
croptypesnice = {'Maize','Wheat','Soybeans','Rice'};
irrtypes = {'rf','ir'};

cd /discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/

lat = importdata('AgGRIDlat.mat');
lon = importdata('AgGRIDlon.mat');

allcroparea = ones(360,720,length(croptypes))*NaN;
area = ncread('maize.nc4','sum');
allcroparea(:,:,1) = area';
area = ncread('wheat.nc4','sum');
allcroparea(:,:,2) = area';
area = ncread('soy.nc4','sum');
allcroparea(:,:,3) = area';
area = ncread('rice.nc4','sum');
allcroparea(:,:,4) = area';

allyears = [1979:2099];

gwlevels = [0.69 1 1.5 2 2.5 3 3.5 4];

startyy = importdata('GWL_StartYear.mat');
endyy = importdata('GWL_EndYear.mat');

%% annual global mean temperature time series

allgmt = ones(length(gcmtypes),length(scentypes),3012)*NaN;
for xgcm=1:length(gcmtypes),
   thisgcm = gcmtypes{xgcm};
   thisgcmnice = gcmtypesnice{xgcm};
   for xscen=1:length(scentypes),
      thisscen = scentypes{xscen};

      globaltas = importdata(['global_tas_mon_' thisgcmnice '_' thisscen '.mat']);
      allgmt(xgcm,xscen,:) = globaltas;
      clear globaltas
   end;
end;

startyyidx = [1:12:3012];
annualgmt = ones(length(gcmtypes),length(scentypes),length(startyyidx))*NaN;
for xgcm=1:length(gcmtypes),
   for xscen=1:length(scentypes),
      for yy=1:length(startyyidx),
         annualgmt(xgcm,xscen,yy) = nanmean(squeeze(allgmt(xgcm,xscen,(startyyidx(yy):(startyyidx(yy)+11)))));
      end;
   end;
end;

annualgmt = annualgmt-273.15;

annualgmtbase = squeeze(annualgmt(:,:,1:51));
annualgmtbasemean = squeeze(nanmean(annualgmtbase,3));

annualgmtanom = annualgmt*NaN;
fullyears = [1850:2100];
for xgcm=1:length(gcmtypes),
   for xscen=1:length(scentypes),
      for yy=1:length(fullyears),
         annualgmtanom(xgcm,xscen,yy) = annualgmt(xgcm,xscen,yy)-annualgmtbasemean(xgcm,xscen);
      end;
   end;
end;

temp20 = annualgmtanom*NaN;
for xgcm=1:length(gcmtypes),
   for xscen=1:length(scentypes),
      for yy=11:(length(fullyears)-9),
         temp20(xgcm,xscen,yy) = squeeze(nanmean(squeeze(annualgmtanom(xgcm,xscen,((yy-10):(yy+9))))));
      end;
   end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  add HadCRUT4 data
%%%%  annual median over regions relative to 1961-1990 mean

crufile = csvread('HadCRUT.5.0.1.0.analysis.summary_series.global.annual.csv',1,0);
cruyears = squeeze(crufile(:,1));
crutemps = squeeze(crufile(:,2));
cruyyidx = 1:51;  %% index of 1850-1900
crudata = crutemps-(nanmean(crutemps(cruyyidx)));

cols = ones(5,3)*NaN;
cols(1,:) = [0 0 128];
cols(2,:) = [245 130 48];
cols(3,:) = [0 130 200];
cols(4,:) = [202 202 43];
cols(5,:) = [128 0 0];
cols = cols/255;

clear h
f = figure;
for xgcm=1:length(gcmtypes),
   h(xgcm) = plot(fullyears,squeeze(annualgmtanom(xgcm,1,:)),'-','color',cols(xgcm,:),'linewidth',1.2); hold on;
   plot(fullyears,squeeze(annualgmtanom(xgcm,2,:)),'--','color',cols(xgcm,:),'linewidth',1.2); hold on;
   plot(fullyears,squeeze(annualgmtanom(xgcm,3,:)),':','color',cols(xgcm,:),'linewidth',1.2); hold on;
end;
h(xgcm+1) = plot(fullyears,(squeeze(annualgmtanom(xgcm,2,:))*NaN),'-','color',[0.4 0.4 0.4],'linewidth',1.2); hold on;
h(xgcm+2) = plot(fullyears,(squeeze(annualgmtanom(xgcm,2,:))*NaN),'--','color',[0.4 0.4 0.4],'linewidth',1.2); hold on;
h(xgcm+3) = plot(fullyears,(squeeze(annualgmtanom(xgcm,2,:))*NaN),':','color',[0.4 0.4 0.4],'linewidth',1.2); hold on;
h(xgcm+4) = plot(cruyears,crudata,'-','color','k','linewidth',1.6); hold on; %[1 0.882352 0.098039]
plot([1850 2100],[0.69 0.69],'--','color',[0.6 0.6 0.6]); hold on;
plot([1850 2100],[1 1],'--','color',[0.6 0.6 0.6]); hold on;
plot([1850 2100],[1.5 1.5],'--','color',[0.6 0.6 0.6]); hold on;
plot([1850 2100],[2 2],'--','color',[0.6 0.6 0.6]); hold on;
plot([1850 2100],[2.5 2.5],'--','color',[0.6 0.6 0.6]); hold on;
plot([1850 2100],[3 3],'--','color',[0.6 0.6 0.6]); hold on;
plot([1850 2100],[3.5 3.5],'--','color',[0.6 0.6 0.6]); hold on;
plot([1850 2100],[4 4],'--','color',[0.6 0.6 0.6]); hold on;
ylabel(['Global Warming Level (' char(176) 'C)']);
%title(['Annual Global Mean Temperature Anomaly from 1850-1900 Mean']);
ah1 = gca;
ah2 = axes('position',get(gca,'position'), 'visible','off');
legendnames = {'GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL','SSP1-2.6','SSP3-7.0','SSP5-8.5','Obs (HadCRUT5)'};
legend(h(:),legendnames,'location','northwest');
%legendnames1 = {'GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL'};
%legendnames2 = {'SSP1-2.6','SSP3-7.0','SSP5-8.5','Obs (HadCRUT5)'};
%legend(h(1:5), 'location', [0.2 0.75 0.15 0.05], legendnames1); hold on;
%legend(h(6:end), 'location', [0.4 0.75 0.15 0.05], legendnames2); hold on;
print(f,'-dpng',['/home/mmphill2/figures/GGCMI/global_tas_anomaly_allgcm_1850-2100_legend.png']);
print(f,'-depsc',['/home/mmphill2/figures/GGCMI/global_tas_anomaly_allgcm_1850-2100_legend.eps']);
close all

%% CO2 response by GWL

co2levels = ones(336,2,3)*NaN;
ssp1co2 = csvread('CO2_RCP2.6.csv',1,0);
ssp3co2 = csvread('CO2_RCP7.0.csv',1,0);
ssp5co2 = csvread('CO2_RCP8.5.csv',1,0);
co2levels(:,:,1) = ssp1co2;
co2levels(:,:,2) = ssp3co2;
co2levels(:,:,3) = ssp5co2;
co2years = 1765:2100;

co2gwl = startyy*NaN;
yeargwl = startyy*NaN;
for xgcm=1:length(gcmtypes),
   for gw=1:length(gwlevels),
      for xscen=1:length(scentypes),
         yystart = startyy(xgcm,gw,xscen);
	       startidx = find(co2years == yystart);
	       yyend = endyy(xgcm,gw,xscen);
	       endidx = find(co2years == yyend);
	       co2gwl(xgcm,gw,xscen) = squeeze(nanmean(squeeze(co2levels(startidx:endidx,2,xscen))));
	       yeargwl(xgcm,gw,xscen) = nanmedian([yystart:yyend]);
      end;
   end;
end;

col = ones(5,3)*NaN;
col(1,:) = [0 0 128];
col(2,:) = [245 130 48];
col(3,:) = [0 130 200];
col(4,:) = [202 202 43];
col(5,:) = [128 0 0];
col = col/255;

f = figure;
clear h
for xgcm=1:length(gcmtypes),
   for gw=1:length(gwlevels),
      plot(gwlevels(gw),co2gwl(xgcm,gw,1),'s','markersize',9,'color',squeeze(col(xgcm,:))); hold on;
      plot(gwlevels(gw),co2gwl(xgcm,gw,2),'o','markersize',9,'color',squeeze(col(xgcm,:))); hold on;
      plot(gwlevels(gw),co2gwl(xgcm,gw,3),'^','markersize',9,'color',squeeze(col(xgcm,:))); hold on;
   end;
end;
h(1) = plot(1,NaN,'.','markersize',20,'color',squeeze(col(1,:))); hold on;
h(2) = plot(1,NaN,'.','markersize',20,'color',squeeze(col(2,:))); hold on;
h(3) = plot(1,NaN,'.','markersize',20,'color',squeeze(col(3,:))); hold on;
h(4) = plot(1,NaN,'.','markersize',20,'color',squeeze(col(4,:))); hold on;
h(5) = plot(1,NaN,'.','markersize',20,'color',squeeze(col(5,:))); hold on;
h(6) = plot(1,NaN,'s','markersize',8,'color',[0.5 0.5 0.5]); hold on;
h(7) = plot(1,NaN,'o','markersize',8,'color',[0.5 0.5 0.5]); hold on;
h(8) = plot(1,NaN,'^','markersize',8,'color',[0.5 0.5 0.5]); hold on;
xlim([0.5 4.2])
xlabel(['Global Warming Level (' char(176) 'C)']);
ylabel(['CO_2 (ppm)']);
%title(['CO2 Level by Global Warming Level for All GCMs/Scenarios']);
legendnames = {'GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL','SSP1-2.6','SSP3-7.0','SSP5-8.5'};
legend(h(:),legendnames,'location','northwest');
print(f,'-dpng',['/home/mmphill2/figures/GGCMI/globalco2_gwl.png']);
print(f,'-depsc',['/home/mmphill2/figures/GGCMI/globalco2_gwl.eps']);
close all

%%%%% crossing year of each GWL
%%%%% version with GMT data

starty = find(fullyears==1980);
subsmoothedtemp = squeeze(temp20(:,:,starty:end));

f = figure;
clear h
for xgcm=1:length(gcmtypes),
   for gw=1:length(gwlevels),
      plot(gwlevels(gw),yeargwl(xgcm,gw,1),'s','markersize',9,'color',squeeze(col(xgcm,:))); hold on;
      plot(gwlevels(gw),yeargwl(xgcm,gw,2),'o','markersize',9,'color',squeeze(col(xgcm,:))); hold on;
      plot(gwlevels(gw),yeargwl(xgcm,gw,3),'^','markersize',9,'color',squeeze(col(xgcm,:))); hold on;
   end;
end;
for xgcm=1:length(gcmtypes),
   gcmsmoothedtemp = squeeze(subsmoothedtemp(xgcm,1,:))';
   line(gcmsmoothedtemp,1980:2100,'LineStyle','-','LineWidth',1.3,'color',squeeze(col(xgcm,:))); hold on;
   gcmsmoothedtemp = squeeze(subsmoothedtemp(xgcm,2,:))';
   line(gcmsmoothedtemp,1980:2100,'LineStyle','--','LineWidth',1.3,'color',squeeze(col(xgcm,:))); hold on;
   gcmsmoothedtemp = squeeze(subsmoothedtemp(xgcm,3,:))'; %'
   line(gcmsmoothedtemp,1980:2100,'LineStyle',':','LineWidth',1.3,'color',squeeze(col(xgcm,:))); hold on;
end;
h(1) = plot(1,NaN,'.','markersize',20,'color',squeeze(col(1,:))); hold on;
h(2) = plot(1,NaN,'.','markersize',20,'color',squeeze(col(2,:))); hold on;
h(3) = plot(1,NaN,'.','markersize',20,'color',squeeze(col(3,:))); hold on;
h(4) = plot(1,NaN,'.','markersize',20,'color',squeeze(col(4,:))); hold on;
h(5) = plot(1,NaN,'.','markersize',20,'color',squeeze(col(5,:))); hold on;
h(6) = plot(1:2,[NaN NaN],'-s','LineWidth',1.3,'color','k'); hold on;
h(7) = plot(1:2,[NaN NaN],'--o','LineWidth',1.3,'color','k'); hold on;
h(8) = plot(1:2,[NaN NaN],':^','LineWidth',1.3,'color','k'); hold on;
xlim([0.5 4.2])
plot([0.69 0.69],[1980 2100],'--','color',[0.6 0.6 0.6]); hold on;
plot([1 1],[1980 2100],'--','color',[0.6 0.6 0.6]); hold on;
plot([1.5 1.5],[1980 2100],'--','color',[0.6 0.6 0.6]); hold on;
plot([2 2],[1980 2100],'--','color',[0.6 0.6 0.6]); hold on;
plot([2.5 2.5],[1980 2100],'--','color',[0.6 0.6 0.6]); hold on;
plot([3 3],[1980 2100],'--','color',[0.6 0.6 0.6]); hold on;
plot([3.5 3.5],[1980 2100],'--','color',[0.6 0.6 0.6]); hold on;
plot([4 4],[1980 2100],'--','color',[0.6 0.6 0.6]); hold on;
plot([0.5 4.2],[2091 2091],'-','linewidth',1,'color',[1 1 1]); hold on;
plot([0.5 4.2],[2100 2100],'-','linewidth',1,'color',[1 1 1]); hold on;
%r1 = [0.5 4.2];
%r2 = [r1, fliplr(r1)];
%inBetween = [[2091 2091], fliplr([2100 2100])];
%fill(r2, inBetween,[0.85 0.85 0.85]); hold on;
xlabel(['Global Warming Level (' char(176) 'C)']);
ylabel(['Year']);
%title(['Middle Year of Crossing by Global Warming Level for All GCMs/Scenarios']);
legendnames = {'GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL','SSP1-2.6','SSP3-7.0','SSP5-8.5'};
legend(h(:),legendnames,'location','southeast');
print(f,'-dpng',['/home/mmphill2/figures/GGCMI/globalyear_gwl_gmtdata.png']);
print(f,'-depsc',['/home/mmphill2/figures/GGCMI/globalyear_gwl_gmtdata.eps']);

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% global production change by year

% For each GCM/SSP line use median across GGCMs at each year
% so it looks the same as the GMT plot
% use more inclusive GGCM subset
% so SSP370 line may be different than 126 and 585

for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   thiscropnice = croptypesnice{xcrop};

   if (strcmp(thiscrop,'mai')),
      modeltypes2 = {'crover','cygma1p74','epic-iiasa','isam','ldndc','lpjml','pepic','pdssat','promet','simplace-lintul5'};
   end;
   if (strcmp(thiscrop,'ric')),
      modeltypes2 = {'cygma1p74','epic-iiasa','isam','ldndc','lpjml','pepic','pdssat','promet'};
   end;
   if (strcmp(thiscrop,'soy')),
      modeltypes2 = {'crover','cygma1p74','epic-iiasa','isam','ldndc','lpjml','pepic','pdssat','promet','simplace-lintul5'};
   end;
   if (strcmp(thiscrop,'whe')),
      modeltypes2 = {'crover','epic-iiasa','isam','ldndc','lpjml','pepic','pdssat','promet','simplace-lintul5'};
   end;
   allproductanom = ones(length(modeltypes2),length(gcmtypes),length(allyears),length(scentypes))*NaN;

   for xmodel=1:length(modeltypes2),
      thisggcm = modeltypes2{xmodel};

      scendata = ['total_global_production_' thisggcm '_' thiscrop '_processed.mat'];
      if(exist(scendata)),
        allscendata = importdata(['total_global_production_' thisggcm '_' thiscrop '_processed.mat']);
        for xgcm=1:length(gcmtypes),
           for xscen=1:length(scentypes),
              for xyear=1:length(allyears),
                 %% take 1983-2013 mean
                 base = nanmean(squeeze(allscendata(xgcm,5:35,xscen)));
                 %% calculate anomaly
                 allproductanom(xmodel,xgcm,xyear,xscen) = (((squeeze(allscendata(xgcm,xyear,xscen))/base)-1)*100);
              end;
            end;
         end;
      end;
   end;
   medianproductanom = squeeze(nanmedian(allproductanom));
   medianproductanom = permute(medianproductanom,[1 3 2]);

   cols = ones(5,3)*NaN;
   cols(1,:) = [0 0 128];
   cols(2,:) = [245 130 48];
   cols(3,:) = [0 130 200];
   cols(4,:) = [202 202 43];
   cols(5,:) = [128 0 0];
   cols = cols/255;

   clear h
   f = figure;
   for xgcm=1:length(gcmtypes),
      h(xgcm) = plot(allyears,squeeze(medianproductanom(xgcm,1,:)),'-','color',cols(xgcm,:),'linewidth',1.2); hold on;
      plot(allyears,squeeze(medianproductanom(xgcm,2,:)),'--','color',cols(xgcm,:),'linewidth',1.2); hold on;
      plot(allyears,squeeze(medianproductanom(xgcm,3,:)),':','color',cols(xgcm,:),'linewidth',1.2); hold on;
   end;
   h(xgcm+1) = plot(allyears,(squeeze(medianproductanom(xgcm,2,:))*NaN),'-','color',[0.4 0.4 0.4],'linewidth',1.2); hold on;
   h(xgcm+2) = plot(allyears,(squeeze(medianproductanom(xgcm,2,:))*NaN),'--','color',[0.4 0.4 0.4],'linewidth',1.2); hold on;
   h(xgcm+3) = plot(allyears,(squeeze(medianproductanom(xgcm,2,:))*NaN),':','color',[0.4 0.4 0.4],'linewidth',1.2); hold on;
   plot([1979 2099],[0 0],':','color','k'); hold on;
   ylabel(['Global Production Change (%)']);
   xlim([1979 2099])
   ah1 = gca;
   ah2 = axes('position',get(gca,'position'), 'visible','off');
   %legendnames = {'GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL','SSP1-2.6','SSP3-7.0','SSP5-8.5','Obs (HadCRUT5)'};
   %legend(h(:),legendnames,'location','northwest');
   %legendnames1 = {'GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL'};
   %legendnames2 = {'SSP1-2.6','SSP3-7.0','SSP5-8.5','Obs (HadCRUT5)'};
   %legend(h(1:5), 'location', [0.2 0.75 0.15 0.05], legendnames1); hold on;
   %legend(h(6:end), 'location', [0.4 0.75 0.15 0.05], legendnames2); hold on;
   print(f,'-dpng',['/home/mmphill2/figures/GGCMI/global_production_anomaly_' thiscrop '_allgcm_1850-2100.png']);
   print(f,'-depsc',['/home/mmphill2/figures/GGCMI/global_production_anomaly_' thiscrop '_allgcm_1850-2100.eps']);
   close all

end;

%%%%% GWL RESPONSE PATTERNS

for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   if (strcmp(thiscrop,'mai')),
      modeltypes2 = {'crover','cygma1p74','epic-iiasa','isam','ldndc','lpjml','pepic','pdssat','promet','simplace-lintul5'};
   end;
   if (strcmp(thiscrop,'ric')),
      modeltypes2 = {'cygma1p74','epic-iiasa','isam','ldndc','lpjml','pepic','pdssat','promet'};
   end;
   if (strcmp(thiscrop,'soy')),
      modeltypes2 = {'crover','cygma1p74','epic-iiasa','isam','ldndc','lpjml','pepic','pdssat','promet','simplace-lintul5'};
   end;
   if (strcmp(thiscrop,'whe')),
      modeltypes2 = {'crover','epic-iiasa','isam','ldndc','lpjml','pepic','pdssat','promet','simplace-lintul5'};
   end;

   allgwlmap = ones(360,720,length(gwlevels),length(modeltypes2),length(gcmtypes))*NaN;

   for xggcm=1:length(modeltypes2),
      thisggcm = modeltypes2{xggcm};
      for xgcm=1:length(gcmtypes),
         thisgcm = gcmtypes{xgcm};
         thisgcmnice = gcmtypesnice{xgcm};
         xscen = 3; %% SSP585 only
         thisscen = scentypes{xscen};
         thisscennice = scentypesnice{xscen};
         filename = [thisggcm '_' thisgcm '_total_' thiscrop '_production_' thisscen '_processed.mat'];
         if(exist(filename)),
            data = importdata([filename]);
            yystart = startyy(xgcm,1,xscen);
            startidx = find(allyears == yystart);
            yyend = endyy(xgcm,1,xscen);
            endidx = find(allyears == yyend);
            basedata = squeeze(nanmean(squeeze(data(:,:,startidx:endidx)),3));
            for gw=1:length(gwlevels),
               yystart = startyy(xgcm,gw,xscen);
	             startidx = find(allyears == yystart);
	             yyend = endyy(xgcm,gw,xscen);
	             endidx = find(allyears == yyend);
	             futdata = squeeze(nanmean(squeeze(data(:,:,startidx:endidx)),3));
	             allgwlmap(:,:,gw,xggcm,xgcm) = (((futdata./basedata)-1)*100);
	          end; % gw
	       end; % exist
	       clear data
      end; % gcm
   end; % ggcm
   clear mediangwl
   mediangwl = squeeze(nanmedian(squeeze(allgwlmap(:,:,:,:)),4));
   save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/median_total_' thiscrop '_production_' thisscen '_processed.mat'],'mediangwl');
end;  % crop


allcropmask = ones(360,720,length(croptypes))*NaN;

rfmaiarea = ncread('maize.nc4','rainfed');
irmaiarea = ncread('maize.nc4','irrigated');
for ii=1:720,
for jj=1:360,
if ((rfmaiarea(ii,jj)>10)||(irmaiarea(ii,jj)>10)),
allcropmask(jj,ii,1) = 1;
end;
end;
end;

rfwhearea = ncread('wheat.nc4','rainfed');
irwhearea = ncread('wheat.nc4','irrigated');
for ii=1:720,
for jj=1:360,
if ((rfwhearea(ii,jj)>10)||(irwhearea(ii,jj)>10)),
allcropmask(jj,ii,2) = 1;
end;
end;
end;

rfsoyarea = ncread('soy.nc4','rainfed');
irsoyarea = ncread('soy.nc4','irrigated');
for ii=1:720,
for jj=1:360,
if ((rfsoyarea(ii,jj)>10)||(irsoyarea(ii,jj)>10)),
allcropmask(jj,ii,3) = 1;
end;
end;
end;

rfricarea = ncread('rice.nc4','rainfed');
irricarea = ncread('rice.nc4','irrigated');
for ii=1:720,
for jj=1:360,
if ((rfricarea(ii,jj)>10)||(irricarea(ii,jj)>10)),
allcropmask(jj,ii,4) = 1;
end;
end;
end;

newcropmasks = importdata('GGCMI_CropMasks_processed.mat');

%%%%%%% GWL response curve types:
%% 1. overall INCREASE + LINEAR
%% 2. overall INCREASE + PEAK AND DECLINE
%% 3. overall INCREASE + DIP AND RISE
%% 4. overall NEUTRAL + LINEAR
%% 5. overall NEUTRAL + PEAK AND DECLINE
%% 6. overall NEUTRAL + DIP AND RISE
%% 7. overall DECREASE + LINEAR
%% 8. overall DECREASE + PEAK AND DECLINE
%% 9. overall DECREASE + DIP AND RISE

patterns = ones(9,3)*NaN;
patterns(1,:) = [0 0 0.8];
patterns(2,:) = [0.2 0.6 1];
patterns(3,:) = [0.2 0.8 0.6];
patterns(4,:) = [0.6 0 0.6];
patterns(5,:) = [0.8 0.6 1];
patterns(6,:) = [1 0.4 0.6];
patterns(7,:) = [0.8 0 0];
patterns(8,:) = [1 0.6 0.2];
patterns(9,:) = [0.87 0.87 0.27];

for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   thiscropnice = croptypesnice{xcrop};

   medianmap = importdata(['median_total_' thiscrop '_production_ssp585_processed.mat']);

   categories = lat*NaN;

   for ii=1:size(lon,2),
      for jj=1:size(lat,1),

         %%% maximum GWL3.5, not using GWL4

         point = squeeze(medianmap(jj,ii,1:7));
         point(isnan(point)) = [];
         point(point == Inf) = [];
         point(point == -Inf) = [];
         if (length(point)==7),
            maxidx = find(point==nanmax(point));
            minidx = find(point==nanmin(point));

            if (point(end)>=(point(1)+3)),  % overall increase

               % for repeating occurrences of ZERO
               if (length(minidx)>1),
                  if (nanmin(point)==0),
                     if (ismember(1,minidx)),
                        minidx = 1; % use start point as minimum
                     end;
                  end;
               end;

               % for repeating occurrences of MAXIMUM VALUE
               if (length(maxidx)>1),
                  if (ismember(7,maxidx)),
                      maxidx = 7; % use end point as maximum
                  else %% repeating peak in the middle of the set
                  %% for these use the first instance
                     maxidx = maxidx(1);
                  end;
               end;

               if ((minidx==1)&&(maxidx==7)),
                  categories(jj,ii) = 1; % increase + linear
               end;
               if ((maxidx~=1)&&(maxidx~=7)),
                  categories(jj,ii) = 2;  % increase + peak and decline
               end;
               if ((minidx~=1)&&(minidx~=7)),
                  categories(jj,ii) = 3;  % increase + dip and rise
               end;
               if ((maxidx~=1)&&(maxidx~=7)&&(minidx~=1)&&(minidx~=7)),
                  % increase with both dip + rise AND peak + decline
                  % so we take the one that is furthest from zero (start)
                  minval = nanmin(point);
                  maxval = nanmax(point);
                  if (abs(minval)>abs(maxval)),
                     categories(jj,ii) = 3;  % increase + dip and rise
                  end;
                  if (abs(minval)<abs(maxval)),
                     categories(jj,ii) = 2;  % increase + peak and decline
                  end;
               end;
            end;  % overall increase

            if (point(end)<=(point(1)-3)),  % overall decrease

               % for repeating occurrences of MAXIMUM
               if (length(maxidx)>1),
                  if (nanmax(point)==0),  % usually zero
                     if (ismember(1,maxidx)),
                        maxidx = 1; % use start point as maximum
                     end;
                  else   % if NOT zero
                     if (ismember(1,maxidx)),
                        maxidx = 1; % same, use start point as maximum
                     else
                        maxidx = maxidx(1); % use one of the middle points
                     end;
                  end;
               end;

               % for repeating occurrences of MINIMUM
               if (length(minidx)>1),
                  if (ismember(7,minidx)), %% often repeating end points
                     minidx = 7;
                  end;
               end;

               if ((maxidx==1)&&(minidx==7)),
                  categories(jj,ii) = 7;  % decrease + linear
               end;
               if ((maxidx~=1)&&(maxidx~=7)),
                  categories(jj,ii) = 8;  % decrease + peak and decline
               end;
               if ((minidx~=1)&&(minidx~=7)),
                  categories(jj,ii) = 9;  % decrease + dip and rise
               end;
               if ((maxidx~=1)&&(maxidx~=7)&&(minidx~=1)&&(minidx~=7)),
                  % decrease with both dip + rise AND peak + decline
                  % so we take the one that is furthest from zero (start)
                  minval = nanmin(point);
                  maxval = nanmax(point);
                  if (abs(minval)>abs(maxval)),
                     categories(jj,ii) = 9;  % decrease + dip and rise
                  end;
                  if (abs(minval)<abs(maxval)),
                     categories(jj,ii) = 8;  % decrease + peak and decline
                  end;
               end;
            end;  % overall decrease

            if ((point(end)<(point(1)+3))&&(point(end)>(point(1)-3))),  % overall neutral

               % for repeating occurrences of ZERO
               if (length(minidx)>1),
                  if (nanmin(point)==0),
                     if (ismember(1,minidx)),
                        minidx = 1; % use start point as minimum
                     end;
                  end;
               end;
               if (length(maxidx)>1),
                  if (nanmax(point)==0),
                     if (ismember(1,maxidx)),
                        maxidx = 1; % use start point as maximum
                     end;
                  else %% rarely a non-zero duplicate maximum
                  %% but in these cases use the first as the max
                     maxidx = maxidx(1);
                  end;
               end;

               if ((maxidx==1)&&(minidx==7)),
                  categories(jj,ii) = 4;  % neutral + linear
               end;
               if ((minidx==1)&&(maxidx==7)),
                  categories(jj,ii) = 4;  % neutral + linear
               end;
               if ((maxidx~=1)&&(maxidx~=7)),
                  if ((nanmax(point)>(point(1)+3))||(nanmax(point)>(point(end)+3))),
                     categories(jj,ii) = 5;  % flat + peak and decline
                  else
                     categories(jj,ii) = 4;  % neutral + linear
                  end;
               end;
               if ((minidx~=1)&&(minidx~=7)),
                  if ((nanmin(point)<(point(1)-3))||(nanmin(point)<(point(end)-3))),
                     categories(jj,ii) = 6;  % flat + dip and rise
                  else
                     categories(jj,ii) = 4;  % neutral + linear
                  end;
               end;
               if (sum(point)==0),  % all zeroes
                  categories(jj,ii) = 4;  % neutral + linear
               end;
            end;  % overall neutral

         end; % length
      end;  % jj
   end;   % ii

   f = figure; colormap(patterns);
   acr_pcolormapr5(categories,lat,lon,[-60 90],[-180 180],[0.5 9.5],['Median GCM/GGCM SSP585 ' thiscropnice ' Production Change Response Pattern By GWL']);
   ax2 = gca;
   colorbar('delete');
   cb2 = colorbar('peer',ax2,'southoutside');
   set(cb2,'xtick',[1:9],'xticklabel',{'Linear \newline{Increase}','Peak \newline{Increase}','Dip \newline{Increase}','Linear \newline{Neutral}','Peak \newline{Neutral}','Dip \newline{Neutral}','Linear \newline{Decrease}','Peak \newline{Decrease}','Dip \newline{Decrease}'});
   print(f,'-dpng',['/home/mmphill2/figures/GGCMI/median_' thiscrop '_ssp585_gwl_response_pattern_default.png']);
   print(f,'-depsc',['/home/mmphill2/figures/GGCMI/median_' thiscrop '_ssp585_gwl_response_pattern_default.eps']);
   close all

end;

%%%%%%% central Julian Date of growing seasons
%%%%%%% use same processed crop masks

croptypes = {'mai','ri1','ri2','soy','swh','wwh'};
croptypesnice = {'Maize','Rice1','Rice2','Soybeans','SpringWheat','WinterWheat'};
irrtypes = {'rf','ir'};

cd /css/giss-cig/AgGRID/CropModels/inputs/growing_season/

pday = ones(size(lat,1),size(lon,2),length(croptypes),length(irrtypes))*NaN;
pday(:,:,1,1) = permute(ncread('mai_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','planting_day'),[2 1]);
pday(:,:,1,2) = permute(ncread('mai_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','planting_day'),[2 1]);
pday(:,:,2,1) = permute(ncread('ri1_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','planting_day'),[2 1]);
pday(:,:,2,2) = permute(ncread('ri1_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','planting_day'),[2 1]);
pday(:,:,3,1) = permute(ncread('ri2_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','planting_day'),[2 1]);
pday(:,:,3,2) = permute(ncread('ri2_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','planting_day'),[2 1]);
pday(:,:,4,1) = permute(ncread('soy_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','planting_day'),[2 1]);
pday(:,:,4,2) = permute(ncread('soy_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','planting_day'),[2 1]);
pday(:,:,5,1) = permute(ncread('swh_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','planting_day'),[2 1]);
pday(:,:,5,2) = permute(ncread('swh_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','planting_day'),[2 1]);
pday(:,:,6,1) = permute(ncread('wwh_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','planting_day'),[2 1]);
pday(:,:,6,2) = permute(ncread('wwh_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','planting_day'),[2 1]);

gslength = ones(size(lat,1),size(lon,2),length(croptypes),length(irrtypes))*NaN;
gslength(:,:,1,1) = permute(ncread('mai_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','growing_season_length'),[2 1]);
gslength(:,:,1,2) = permute(ncread('mai_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','growing_season_length'),[2 1]);
gslength(:,:,2,1) = permute(ncread('ri1_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','growing_season_length'),[2 1]);
gslength(:,:,2,2) = permute(ncread('ri1_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','growing_season_length'),[2 1]);
gslength(:,:,3,1) = permute(ncread('ri2_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','growing_season_length'),[2 1]);
gslength(:,:,3,2) = permute(ncread('ri2_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','growing_season_length'),[2 1]);
gslength(:,:,4,1) = permute(ncread('soy_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','growing_season_length'),[2 1]);
gslength(:,:,4,2) = permute(ncread('soy_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','growing_season_length'),[2 1]);
gslength(:,:,5,1) = permute(ncread('swh_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','growing_season_length'),[2 1]);
gslength(:,:,5,2) = permute(ncread('swh_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','growing_season_length'),[2 1]);
gslength(:,:,6,1) = permute(ncread('wwh_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','growing_season_length'),[2 1]);
gslength(:,:,6,2) = permute(ncread('wwh_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','growing_season_length'),[2 1]);

hday = ones(size(lat,1),size(lon,2),length(croptypes),length(irrtypes))*NaN;
hday(:,:,1,1) = permute(ncread('mai_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','maturity_day'),[2 1]);
hday(:,:,1,2) = permute(ncread('mai_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','maturity_day'),[2 1]);
hday(:,:,2,1) = permute(ncread('ri1_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','maturity_day'),[2 1]);
hday(:,:,2,2) = permute(ncread('ri1_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','maturity_day'),[2 1]);
hday(:,:,3,1) = permute(ncread('ri2_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','maturity_day'),[2 1]);
hday(:,:,3,2) = permute(ncread('ri2_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','maturity_day'),[2 1]);
hday(:,:,4,1) = permute(ncread('soy_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','maturity_day'),[2 1]);
hday(:,:,4,2) = permute(ncread('soy_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','maturity_day'),[2 1]);
hday(:,:,5,1) = permute(ncread('swh_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','maturity_day'),[2 1]);
hday(:,:,5,2) = permute(ncread('swh_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','maturity_day'),[2 1]);
hday(:,:,6,1) = permute(ncread('wwh_rf_ggcmi_crop_calendar_phase3_v1.01.nc4','maturity_day'),[2 1]);
hday(:,:,6,2) = permute(ncread('wwh_ir_ggcmi_crop_calendar_phase3_v1.01.nc4','maturity_day'),[2 1]);

rfpday = squeeze(pday(:,:,:,1));
rfpday(rfpday <= 0) = NaN;
rfhday = squeeze(hday(:,:,:,1));
rfhday(rfhday <= 0) = NaN;
rfgslength = squeeze(gslength(:,:,:,1));
rfgslength(rfgslength <= 0) = NaN;

cd /discover/nobackup/projects/giss_ana/users/mmphill2/GWLPaper/

oldhsv = hsv;
close all
newhsv = oldhsv*NaN;
colorscale = [48:64 1:47];
for cc=1:length(colorscale),
   newhsv(cc,:) = oldhsv(colorscale(cc),:);
end;

newcropmasks = importdata('GGCMI_CropMasks_processed.mat');
swhmask = permute(ncread('winter_and_spring_wheat_areas_phase3.nc4','swh_mask'),[2 1]);
wwhmask = permute(ncread('winter_and_spring_wheat_areas_phase3.nc4','wwh_mask'),[2 1]);
swhmask(swhmask==0) = NaN;
wwhmask(wwhmask==0) = NaN;

for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   thiscropnice = croptypesnice{xcrop};

   if (strcmp(thiscrop,'mai')),
      thiscropmask = squeeze(newcropmasks(:,:,1));
   end;
   if (strcmp(thiscrop,'swh')),
      thiscropmask = swhmask;
   end;
   if (strcmp(thiscrop,'wwh')),
      thiscropmask = wwhmask;
   end;
   if (strcmp(thiscrop,'soy')),
      thiscropmask = squeeze(newcropmasks(:,:,3));
   end;
   if (strcmp(thiscrop,'ri1')),
      thiscropmask = squeeze(newcropmasks(:,:,4));
   end;
   if (strcmp(thiscrop,'ri2')),
      thiscropmask = squeeze(newcropmasks(:,:,4));
   end;

   gsmap = lat*NaN;

   for ii=1:size(lon,2),
      for jj=1:size(lat,1),

         if (~isnan(rfgslength(jj,ii,xcrop)) && (rfgslength(jj,ii,xcrop)>0) && (rfpday(jj,ii,xcrop)>0)),

            thispday = ceil(squeeze(rfpday(jj,ii,xcrop)));
            thishday = ceil(squeeze(rfhday(jj,ii,xcrop)));

	          if(thishday<thispday),
	             thisgslength = ((365-thispday)+(thishday));
	          end;
	          if(thishday>thispday),
	             thisgslength = thishday-thispday;
	          end;

	          midseason = (thispday+(thisgslength/2));
	          if(midseason > 365),
	             gsmap(jj,ii) = midseason-365;
	          elseif(midseason < 1),
	             gsmap(jj,ii) = 365+midseason;
	          else
	             gsmap(jj,ii) = midseason;
	          end;
         end;

      end;
   end;

   f = figure; colormap(newhsv);
   acr_pcolormapr5(gsmap.*thiscropmask,lat,lon,[-60 90],[-180 180],[0 365],['Middle Day of Growing Season for Rainfed ' thiscropnice],['Julian date']);
   print(f, '-dpng', ['/home/mmphill2/midday_growing_season_rf_' thiscrop '_phase3.png']);
   print(f, '-depsc', ['/home/mmphill2/midday_growing_season_rf_' thiscrop '_phase3.eps']);
   close all

end;
