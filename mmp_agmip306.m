%		               mmp_agmip306
%
%  This script creates crop season temperature and
%  precipitation change maps for each GCM, SSP-RCP
%  scenario, and global warming level (time period)
%  for all crops and irrigation types using ISIMIP3
%  CMIP6 data. This is an updated version of
%  mmp_agmip206.m.
%
%  -- Creates a map of each temperature and
%  precipitation growing season change between the
%  1850-1900 mean and each GWL time period mean
%  for all combinations of GCM, crop, scenario and GWL.
%
%  -- Creates a map of each temperature and
%  precipitation growing season change between the
%  1850-1900 mean and each GWL time period mean
%  showing the mean across GCMs for all combinations
%  of crop, scenario and GWL.
%
%  INPUTS: none
%  OUTPUTS: none
%
%
%  This is written to be run on Discover.
%
%
%		                             author: Meridel Phillips
%                                mmp2192@columbia.edu
%				                         date:	8/18/2022
%
function mmp_agmip306();
%--------------------------------------------------
%--------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Begin Debug
%% End Debug

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /home/mmphill2/matlab/

modeltypes = {'acea','crover','cygma1p74','dssat-pythia','epic-iiasa','isam','ldndc','lpjml','pdssat','pepic','promet','simplace-lintul5'};
modeltypesnice = {'ACEA','CROVER','CYGMA','DSSAT-Pythia','EPIC-IIASA','ISAM','LandscapeDNDC','LPJmL','pDSSAT','PEPIC','PROMET','SIMPLACE-LINTUL5'};
scentypes = {'ssp126','ssp370','ssp585'};
scentypesnice = {'SSP126','SSP370','SSP585'};
gcmtypes = {'gfdl-esm4','ipsl-cm6a-lr','mpi-esm1-2-hr','mri-esm2-0','ukesm1-0-ll'};
gcmtypesnice = {'GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL'};
croptypes = {'mai','ri1','ri2','soy','swh','wwh'};
croptypesnice = {'Maize','Rice1','Rice2','Soybeans','SpringWheat','WinterWheat'};
irrtypes = {'rf','ir'};

cd /css/giss-cig/AgGRID/CropModels/phase3/matfiles/

lat = importdata('AgGRIDlat.mat');
lon = importdata('AgGRIDlon.mat');

cd /css/giss-cig/AgGRID/CropModels/inputs/growing_season/

%%% growing seasons

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

cd /discover/nobackup/projects/giss_ana/users/mmphill2/GWLPaper/

allyears = [1850:2100];

gwlevels = [0.69 1 1.5 2 2.5 3 3.5 4];

startyy = importdata('GWL_StartYear.mat');
endyy = importdata('GWL_EndYear.mat');

%%%%%% calculate growing season temp / precip
%%%%%% for each GCM, GWL, scenario and crop

ripfs = {'r1i1p1f1','r1i1p1f1','r1i1p1f1','r1i1p1f1','r1i1p1f2'};

histyears = [1850:2014];
subhistyears = [1970:2014];
histstartyears = [1850 1851 1861 1871 1881 1891 1901 1911 1921 1931 1941 1951 1961 1971 1981 1991 2001 2011];
histendyears = [1850 1860 1870 1880 1890 1900 1910 1920 1930 1940 1950 1960 1970 1980 1990 2000 2010 2014];
histfilelength = [1 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 10 4];

yearlength = repmat(365,1,length(histyears));
leapidx = [3:4:length(histyears)];
yearlength(leapidx) = 366;
yearlength(find(histyears == 1900)) = 365;

%%% pull in annual climate data and save
%%% as individual .mat files for faster use

%% temperature
for xgcm=1:length(gcmtypes),
   thisgcm = gcmtypes{xgcm};
   thisens = ripfs{xgcm};
   for yy=1:length(histstartyears),
      fileyrs = histstartyears(yy):histendyears(yy);
      filestr = [thisgcm '_' thisens '_w5e5_historical_tas_global_daily_' num2str(histstartyears(yy)) '_' num2str(histendyears(yy)) '.nc'];
      data = (ncread(filestr,'tas'))-273.15;
      startidx = find(histyears == histstartyears(yy));
      thisfileyylengths = yearlength(startidx:(startidx+histfilelength(yy)-1));
      yyend = cumsum(thisfileyylengths);
      yystart = [1 (yyend+1)];
      yystart(end) = [];
      for insideyear = 1:length(yystart),
         thisyear = fileyrs(insideyear);
         yydata = permute(squeeze(data(:,:,yystart(insideyear):yyend(insideyear))),[2 1 3]);
         if (size(yydata,3)>365),
            yydata(:,:,60) = [];
         end;
         disp([thisgcm ' ' num2str(thisyear)])
         save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/' thisgcm '_daily_tas_historical_' num2str(thisyear) '.mat'],'yydata');
         clear yydata
      end;
      clear data
   end;
end;

%% precipitation
for xgcm=1:length(gcmtypes),
   thisgcm = gcmtypes{xgcm};
   thisens = ripfs{xgcm};
   for yy=1:length(histstartyears),
      fileyrs = histstartyears(yy):histendyears(yy);
      filestr = [thisgcm '_' thisens '_w5e5_historical_pr_global_daily_' num2str(histstartyears(yy)) '_' num2str(histendyears(yy)) '.nc'];
      data = (ncread(filestr,'pr'))*86400;
      startidx = find(histyears == histstartyears(yy));
      thisfileyylengths = yearlength(startidx:(startidx+histfilelength(yy)-1));
      yyend = cumsum(thisfileyylengths);
      yystart = [1 (yyend+1)];
      yystart(end) = [];
      for insideyear = 1:length(yystart),
         thisyear = fileyrs(insideyear);
         yydata = permute(squeeze(data(:,:,yystart(insideyear):yyend(insideyear))),[2 1 3]);
         if (size(yydata,3)>365),
            yydata(:,:,60) = [];
         end;
         disp([thisgcm ' ' num2str(thisyear)])
         save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/' thisgcm '_daily_pr_historical_' num2str(thisyear) '.mat'],'yydata');
         clear yydata
      end;
      clear data
   end;
end;

%%% save annual future climate data

histyears = [2015:2100];
subhistyears = histyears;
histstartyears = [2015 2021 2031 2041 2051 2061 2071 2081 2091];
histendyears = [2020 2030 2040 2050 2060 2070 2080 2090 2100];
histfilelength = [6 10 10 10 10 10 10 10 10];

yearlength = repmat(365,1,length(histyears));
leapidx = [2:4:length(histyears)];
yearlength(leapidx) = 366;
yearlength(find(histyears == 2100)) = 365;

% temperature
for xgcm=1:length(gcmtypes),
   thisgcm = gcmtypes{xgcm};
   thisens = ripfs{xgcm};
   for xscen=2:length(scentypes),
      thisscen = scentypes{xscen};
      for yy=1:length(histstartyears), %% start at 2015
         fileyrs = histstartyears(yy):histendyears(yy);
         filestr = [thisgcm '_' thisens '_w5e5_' thisscen '_tas_global_daily_' num2str(histstartyears(yy)) '_' num2str(histendyears(yy)) '.nc'];
         data = (ncread(filestr,'tas'))-273.15;
         startidx = find(histyears == histstartyears(yy));
         thisfileyylengths = yearlength(startidx:(startidx+histfilelength(yy)-1));
         yyend = cumsum(thisfileyylengths);
         yystart = [1 (yyend+1)];
         yystart(end) = [];
         for insideyear = 1:length(yystart),
            thisyear = fileyrs(insideyear);
            yydata = permute(squeeze(data(:,:,yystart(insideyear):yyend(insideyear))),[2 1 3]);
            if (size(yydata,3)>365),
               yydata(:,:,60) = [];
            end;
            disp([thisgcm ' ' thisscen ' ' num2str(thisyear)])
            save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/' thisgcm '_daily_tas_' thisscen '_' num2str(thisyear) '.mat'],'yydata');
            clear yydata
         end;
         clear data
      end;
   end;
end;

% precipitation
for xgcm=1:length(gcmtypes),
   thisgcm = gcmtypes{xgcm};
   thisens = ripfs{xgcm};
   for xscen=2:length(scentypes),
      thisscen = scentypes{xscen};
      for yy=1:length(histstartyears), %% start at 2015
         fileyrs = histstartyears(yy):histendyears(yy);
         filestr = [thisgcm '_' thisens '_w5e5_' thisscen '_pr_global_daily_' num2str(histstartyears(yy)) '_' num2str(histendyears(yy)) '.nc'];
         data = (ncread(filestr,'pr'))*86400;
         startidx = find(histyears == histstartyears(yy));
         thisfileyylengths = yearlength(startidx:(startidx+histfilelength(yy)-1));
         yyend = cumsum(thisfileyylengths);
         yystart = [1 (yyend+1)];
         yystart(end) = [];
         for insideyear = 1:length(yystart),
            thisyear = fileyrs(insideyear);
            yydata = permute(squeeze(data(:,:,yystart(insideyear):yyend(insideyear))),[2 1 3]);
            if (size(yydata,3)>365),
               yydata(:,:,60) = [];
            end;
            disp([thisgcm ' ' thisscen ' ' num2str(thisyear)])
            save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/' thisgcm '_daily_pr_' thisscen '_' num2str(thisyear) '.mat'],'yydata');
            clear yydata
         end;
         clear data
      end;
   end;
end;

%%%%%% save crop season means for each GCM/GWL/scenario

cd /css/giss-cig/ISIMIP3/matfiles/

% temperature
for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   for xirr=1:length(irrtypes),
      thisirr = irrtypes{xirr};
      for xgcm=1:length(gcmtypes),
         thisgcm = gcmtypes{xgcm};
         thisens = ripfs{xgcm};
         for xscen=1:length(scentypes),
            thisscen = scentypes{xscen};
            cropmeans = ones(size(lat,1),size(lon,2),length(gwlevels),20)*NaN;
            for xgwl=1:length(gwlevels),
               thisgwl = gwlevels(xgwl);
               gwstart = startyy(xgcm,xgwl,xscen);
               gwend = endyy(xgcm,xgwl,xscen);
               gwyears  = gwstart:gwend;
               if (~isnan(gwstart)),
                  for xyear=1:length(gwyears),
                     thisyear = gwyears(xyear);
                     if (thisyear<2015),
                        data = importdata([thisgcm '_daily_tas_historical_' num2str(thisyear) '.mat']);
                     else
                        data = importdata([thisgcm '_daily_tas_' thisscen '_' num2str(thisyear) '.mat']);
                     end;
                     for ii=1:size(lon,2),
                        for jj=1:size(lat,1),
                           if (~isnan(gslength(jj,ii,xcrop,xirr)) && (gslength(jj,ii,xcrop,xirr)>0) && (pday(jj,ii,xcrop,xirr)>0)),
                              seasonstart = ceil(pday(jj,ii,xcrop,xirr)); %% some= 0.5
                              seasonend = ceil(pday(jj,ii,xcrop,xirr)+gslength(jj,ii,xcrop,xirr));
                              if (seasonend>365),
               		               seasonend = seasonend-365;
               		            end;
                              if (seasonend<seasonstart),
                                 cropmeans(jj,ii,xgwl,xyear) = squeeze(nanmean(squeeze(data(jj,ii,[seasonstart:365 1:seasonend]))));
                              else
                                 cropmeans(jj,ii,xgwl,xyear) = squeeze(nanmean(squeeze(data(jj,ii,seasonstart:seasonend))));
                              end;
                           end;
                        end;
                     end;

                     disp([thisirr ' ' thiscrop ' ' thisgcm ' ' thisscen ' ' num2str(thisgwl) ' ' num2str(thisyear)])

                  end; % xyear
               end; % NaN
            end; % gwl
            meangwldata = squeeze(nanmean(cropmeans,4));
            save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/' thisgcm '_' thisirr '_' thiscrop '_tas_' thisscen '_gwlmeans.mat'],'meangwldata');
         end; % scen
      end; % gcm
   end; % crop
end; % irr

% temperature 1850-1900 mean
for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   for xirr=1:length(irrtypes),
      thisirr = irrtypes{xirr};
      for xgcm=1:length(gcmtypes),
         thisgcm = gcmtypes{xgcm};
         thisens = ripfs{xgcm};
         baseyears = [1850:1900];
         cropmeans = ones(size(lat,1),size(lon,2),length(baseyears))*NaN;
         gwstart = 1850;
         gwend = 1900;
         gwyears  = gwstart:gwend;
         for xyear=1:length(gwyears),
            thisyear = gwyears(xyear);
            data = importdata([thisgcm '_daily_tas_historical_' num2str(thisyear) '.mat']);
            for ii=1:size(lon,2),
               for jj=1:size(lat,1),
                  if (~isnan(gslength(jj,ii,xcrop,xirr)) && (gslength(jj,ii,xcrop,xirr)>0) && (pday(jj,ii,xcrop,xirr)>0)),
                        seasonstart = ceil(pday(jj,ii,xcrop,xirr)); %% some= 0.5
                        seasonend = ceil(pday(jj,ii,xcrop,xirr)+gslength(jj,ii,xcrop,xirr));
                        if (seasonend>365),
               		         seasonend = seasonend-365;
               		      end;
                        if (seasonend<seasonstart),
                           cropmeans(jj,ii,xyear) = squeeze(nanmean(squeeze(data(jj,ii,[seasonstart:365 1:seasonend]))));
                        else
                           cropmeans(jj,ii,xyear) = squeeze(nanmean(squeeze(data(jj,ii,seasonstart:seasonend))));
                        end;
                  end;
               end;
            end;

            disp([thisirr ' ' thiscrop ' ' thisgcm ' ' num2str(thisyear)])

         end; % xyear
         meanbasedata = squeeze(nanmean(cropmeans,3));
         save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/' thisgcm '_' thisirr '_' thiscrop '_tas_mean_1850-1900.mat'],'meanbasedata');
      end; % gcm
   end; % crop
end; % irr

% precipitation
for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   for xirr=1:length(irrtypes),
      thisirr = irrtypes{xirr};
      for xgcm=1:length(gcmtypes),
         thisgcm = gcmtypes{xgcm};
         thisens = ripfs{xgcm};
         for xscen=1:length(scentypes),
            thisscen = scentypes{xscen};
            cropmeans = ones(size(lat,1),size(lon,2),length(gwlevels),20)*NaN;
            for xgwl=1:length(gwlevels),
               thisgwl = gwlevels(xgwl);
               gwstart = startyy(xgcm,xgwl,xscen);
               gwend = endyy(xgcm,xgwl,xscen);
               gwyears  = gwstart:gwend;
               if (~isnan(gwstart)),
                  for xyear=1:length(gwyears),
                     thisyear = gwyears(xyear);
                     if (thisyear<2015),
                        data = importdata([thisgcm '_daily_pr_historical_' num2str(thisyear) '.mat']);
                     else
                        data = importdata([thisgcm '_daily_pr_' thisscen '_' num2str(thisyear) '.mat']);
                     end;
                     for ii=1:size(lon,2),
                        for jj=1:size(lat,1),
                           if (~isnan(gslength(jj,ii,xcrop,xirr)) && (gslength(jj,ii,xcrop,xirr)>0) && (pday(jj,ii,xcrop,xirr)>0)),
                              seasonstart = ceil(pday(jj,ii,xcrop,xirr)); %% some= 0.5
                              seasonend = ceil(pday(jj,ii,xcrop,xirr)+gslength(jj,ii,xcrop,xirr));
                              if (seasonend>365),
               		               seasonend = seasonend-365;
               		            end;
                              if (seasonend<seasonstart),
                                 cropmeans(jj,ii,xgwl,xyear) = squeeze(nansum(squeeze(data(jj,ii,[seasonstart:365 1:seasonend]))));
                              else
                                 cropmeans(jj,ii,xgwl,xyear) = squeeze(nansum(squeeze(data(jj,ii,seasonstart:seasonend))));
                              end;
                           end;
                        end;
                     end;

                     disp([thisirr ' ' thiscrop ' ' thisgcm ' ' thisscen ' ' num2str(thisgwl) ' ' num2str(thisyear)])

                  end; % xyear
               end; % NaN
            end; % gwl
            meangwldata = squeeze(nanmean(cropmeans,4));
            save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/' thisgcm '_' thisirr '_' thiscrop '_pr_' thisscen '_gwlmeans.mat'],'meangwldata');
         end; % scen
      end; % gcm
   end; % crop
end; % irr

% precipitation 1850-1900 mean
for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   for xirr=1:length(irrtypes),
      thisirr = irrtypes{xirr};
      for xgcm=1:length(gcmtypes),
         thisgcm = gcmtypes{xgcm};
         thisens = ripfs{xgcm};
         baseyears = [1850:1900];
         cropmeans = ones(size(lat,1),size(lon,2),length(baseyears))*NaN;
         gwstart = 1850;
         gwend = 1900;
         gwyears  = gwstart:gwend;
         for xyear=1:length(gwyears),
            thisyear = gwyears(xyear);
            data = importdata([thisgcm '_daily_pr_historical_' num2str(thisyear) '.mat']);
            for ii=1:size(lon,2),
               for jj=1:size(lat,1),
                   if (~isnan(gslength(jj,ii,xcrop,xirr)) && (gslength(jj,ii,xcrop,xirr)>0) && (pday(jj,ii,xcrop,xirr)>0)),
                        seasonstart = ceil(pday(jj,ii,xcrop,xirr)); %% some= 0.5
                        seasonend = ceil(pday(jj,ii,xcrop,xirr)+gslength(jj,ii,xcrop,xirr));
                        if (seasonend>365),
               		         seasonend = seasonend-365;
               		      end;
                        if (seasonend<seasonstart),
                           cropmeans(jj,ii,xyear) = squeeze(nansum(squeeze(data(jj,ii,[seasonstart:365 1:seasonend]))));
                        else
                           cropmeans(jj,ii,xyear) = squeeze(nansum(squeeze(data(jj,ii,seasonstart:seasonend))));
                        end;
                  end;
               end;
            end;

            disp([thisirr ' ' thiscrop ' ' thisgcm ' ' num2str(thisyear)])

         end; % xyear
         meanbasedata = squeeze(nanmean(cropmeans,3));
         save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/' thisgcm '_' thisirr '_' thiscrop '_pr_mean_1850-1900.mat'],'meanbasedata');
      end; % gcm
   end; % crop
end; % irr

%%%%%% map growing season temp / precip change
%%%%%% for each GCM, GWL, scenario and crop
%%%%%% and mean across GCMs for each GWL/scenario

cd /discover/nobackup/projects/giss_ana/users/mmphill2/GWLPaper/

newcropmasks = importdata('GGCMI_CropMasks_processed.mat');
swhmask = permute(ncread('winter_and_spring_wheat_areas_phase3.nc4','swh_mask'),[2 1]);
wwhmask = permute(ncread('winter_and_spring_wheat_areas_phase3.nc4','wwh_mask'),[2 1]);
swhmask(swhmask==0) = NaN;
wwhmask(wwhmask==0) = NaN;

cd /css/giss-cig/ISIMIP3/matfiles/
% temperature
for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   thiscropnice = croptypesnice{xcrop};
   for xirr=1:length(irrtypes),
      thisirr = irrtypes{xirr};
      meandata = ones(size(lat,1),size(lon,2),length(gcmtypes),length(gwlevels),length(scentypes))*NaN;
      meanbasedata = ones(size(lat,1),size(lon,2),length(gcmtypes))*NaN;
      %thiscroparea = squeeze(allcroparea(:,:,xcrop,xirr));
      %thiscropmask = lat*NaN;
      %for ii=1:720,
      %   for jj=1:360,
      %      if (thiscroparea(jj,ii)>10),
      %         thiscropmask(jj,ii) = 1;
      %      end;
      %   end;
      %end;
      if (strcmp(thiscrop,'mai')),
         thiscropmask = squeeze(newcropmasks(:,:,1));
      end;
      if (strcmp(thiscrop,'swh')),
         %thiscropmask = squeeze(newcropmasks(:,:,2));
         thiscropmask = swhmask;
      end;
      if (strcmp(thiscrop,'wwh')),
         %thiscropmask = squeeze(newcropmasks(:,:,2));
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

      for xgcm=1:length(gcmtypes),
         thisgcm = gcmtypes{xgcm};
         thisgcmnice = gcmtypesnice{xgcm};
         for xscen=1:length(scentypes),
            thisscen = scentypes{xscen};
            thisscennice = scentypesnice{xscen};
            gcmdata = importdata([thisgcm '_' thisirr '_' thiscrop '_tas_' thisscen '_gwlmeans.mat']);
            meandata(:,:,xgcm,:,xscen) = gcmdata;
            gcmbase = importdata([thisgcm '_' thisirr '_' thiscrop '_tas_mean_1850-1900.mat']);
            meanbasedata(:,:,xgcm) = gcmbase;
            for xgwl=2:length(gwlevels),
               if (sum(sum((~isnan(squeeze(gcmdata(:,:,xgwl))))))>0),
                  gcmfut = squeeze(gcmdata(:,:,xgwl));
                  deltagcm = gcmfut-gcmbase;
                  p = importdata('/home/mmphill2/matlab/matfiles/cmap_bluered_zerogrey.mat');
                  f = figure; colormap(p);
                  acr_pcolormapr5(deltagcm.*thiscropmask,lat,lon,[-60 90],[-180 180],[(gwlevels(xgwl)-2) (gwlevels(xgwl)+2)],['Change in ' thisgcmnice ' ' thisirr thiscropnice ' Seasonal Temperature for ' thisscennice ' GWL' num2str(gwlevels(xgwl))],['degree change from 1850-1900 ' thisirr thiscropnice ' season mean']);
                  print(f,'-dpng',['/home/mmphill2/figures/GGCMI/GWL' num2str(gwlevels(xgwl)) '_' thisscennice '_' thisirr thiscropnice '_' thisgcmnice '_TemperatureChange.png']);
                  print(f,'-depsc',['/home/mmphill2/figures/GGCMI/GWL' num2str(gwlevels(xgwl)) '_' thisscennice '_' thisirr thiscropnice '_' thisgcmnice '_TemperatureChange.eps']);
                  close all
               end;
            end; % xgwl
        end; % xscen
     end; % xgcm

     meangcmbase = squeeze(nanmean(meanbasedata,3));
     for xgwl=2:length(gwlevels),
        meangcmfutallscen = squeeze(meandata(:,:,:,xgwl,:));
        meangcmfutallscen = squeeze(nanmean(squeeze(meangcmfutallscen(:,:,:)),3));
        meandeltagcmallscen = meangcmfutallscen-meangcmbase;
        if (sum(sum((~isnan(meandeltagcmallscen))))>0),
           %p = importdata('/home/mmphill2/matlab/matfiles/cmap_bluered_zerogrey.mat');
           p = importdata('/Users/mmphill2/Documents/MATLAB/cmap_bluered_zerogrey.mat');
           f = figure; colormap(p);
           acr_pcolormapr5(meandeltagcmallscen.*thiscropmask,lat,lon,[-60 90],[-180 180],[(gwlevels(xgwl)-2) (gwlevels(xgwl)+2)],['Mean Change in ' thisirr thiscropnice ' Seasonal Temperature for AllScenario GWL' num2str(gwlevels(xgwl))],['degree change from 1850-1900 ' thisirr thiscropnice ' season mean']);
           %print(f,'-dpng',['/home/mmphill2/GWL' num2str(gwlevels(xgwl)) '_allscen_' thisirr thiscropnice '_MeanTemperatureChange.png']);
           %print(f,'-depsc',['/home/mmphill2/GWL' num2str(gwlevels(xgwl)) '_allscen_' thisirr thiscropnice '_MeanTemperatureChange.eps']);
           print(f,'-dpng',['GWL' num2str(gwlevels(xgwl)) '_allscen_' thisirr thiscropnice '_MeanTemperatureChange.png']);
           print(f,'-depsc',['GWL' num2str(gwlevels(xgwl)) '_allscen_' thisirr thiscropnice '_MeanTemperatureChange.eps']);
           close all
        end;

        for xscen=1:length(scentypes),
           thisscen = scentypes{xscen};
           thisscennice = scentypesnice{xscen};

           meangcmfut = squeeze(nanmean(squeeze(meandata(:,:,:,xgwl,xscen)),3));
           meandeltagcm = meangcmfut-meangcmbase;
           if (sum(sum((~isnan(meandeltagcm))))>0),
             %p = importdata('/home/mmphill2/matlab/matfiles/cmap_bluered_zerogrey.mat');
              p = importdata('/Users/mmphill2/Documents/MATLAB/cmap_bluered_zerogrey.mat');
              f = figure; colormap(p);
              acr_pcolormapr5(meandeltagcm.*thiscropmask,lat,lon,[-60 90],[-180 180],[(gwlevels(xgwl)-2) (gwlevels(xgwl)+2)],['Mean Change in ' thisirr thiscropnice ' Seasonal Temperature for ' thisscennice ' GWL' num2str(gwlevels(xgwl))],['degree change from 1850-1900 ' thisirr thiscropnice ' season mean']);
              %print(f,'-dpng',['/home/mmphill2/GWL' num2str(gwlevels(xgwl)) '_' thisscennice '_' thisirr thiscropnice '_MeanTemperatureChange.png']);
              %print(f,'-depsc',['/home/mmphill2/GWL' num2str(gwlevels(xgwl)) '_' thisscennice '_' thisirr thiscropnice '_MeanTemperatureChange.eps']);
              print(f,'-dpng',['GWL' num2str(gwlevels(xgwl)) '_' thisscennice '_' thisirr thiscropnice '_MeanTemperatureChange.png']);
              print(f,'-depsc',['GWL' num2str(gwlevels(xgwl)) '_' thisscennice '_' thisirr thiscropnice '_MeanTemperatureChange.eps']);

              close all
           end;
        end; % xgwl
     end; % xscen

   end;
end;

cd /css/giss-cig/ISIMIP3/matfiles/

% precipitation
for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   thiscropnice = croptypesnice{xcrop};
   for xirr=1:length(irrtypes),
      thisirr = irrtypes{xirr};
      meandata = ones(size(lat,1),size(lon,2),length(gcmtypes),length(gwlevels),length(scentypes))*NaN;
      meanbasedata = ones(size(lat,1),size(lon,2),length(gcmtypes))*NaN;
      %thiscroparea = squeeze(allcroparea(:,:,xcrop,xirr));
      %thiscropmask = lat*NaN;
      %for ii=1:720,
      %   for jj=1:360,
      %      if (thiscroparea(jj,ii)>10),
      %         thiscropmask(jj,ii) = 1;
      %      end;
      %   end;
      %end;
      if (strcmp(thiscrop,'mai')),
         thiscropmask = squeeze(newcropmasks(:,:,1));
      end;
      if (strcmp(thiscrop,'swh')),
         %thiscropmask = squeeze(newcropmasks(:,:,2));
         thiscropmask = swhmask;
      end;
      if (strcmp(thiscrop,'wwh')),
         %thiscropmask = squeeze(newcropmasks(:,:,2));
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

      for xgcm=1:length(gcmtypes),
         thisgcm = gcmtypes{xgcm};
         thisgcmnice = gcmtypesnice{xgcm};
         for xscen=1:length(scentypes),
            thisscen = scentypes{xscen};
            thisscennice = scentypesnice{xscen};
            gcmdata = importdata([thisgcm '_' thisirr '_' thiscrop '_pr_' thisscen '_gwlmeans.mat']);
            meandata(:,:,xgcm,:,xscen) = gcmdata;
            gcmbase = importdata([thisgcm '_' thisirr '_' thiscrop '_pr_mean_1850-1900.mat']);
            meanbasedata(:,:,xgcm) = gcmbase;
            for xgwl=2:length(gwlevels),
               if (sum(sum((~isnan(squeeze(gcmdata(:,:,xgwl))))))>0),
                  gcmfut = squeeze(gcmdata(:,:,xgwl));
                  deltagcm = (((gcmfut./gcmbase)-1)*100);
                  p = importdata('/home/mmphill2/matlab/matfiles/cmap_percentchange.mat');
                  f = figure; colormap(p);
                  acr_pcolormapr5(deltagcm.*thiscropmask,lat,lon,[-60 90],[-180 180],[-30 30],['Change in ' thisgcmnice ' ' thisirr thiscropnice ' Seasonal Precipitation for ' thisscennice ' GWL' num2str(gwlevels(xgwl))],['percent change from 1850-1900 ' thisirr thiscropnice ' season mean']);
                  print(f,'-dpng',['/home/mmphill2/GWL' num2str(gwlevels(xgwl)) '_' thisscennice '_' thisirr thiscropnice '_' thisgcmnice '_PrecipitationChange.png']);
                  print(f,'-depsc',['/home/mmphill2/GWL' num2str(gwlevels(xgwl)) '_' thisscennice '_' thisirr thiscropnice '_' thisgcmnice '_PrecipitationChange.eps']);
                  close all
               end;
            end; % xgwl
        end; % xscen
     end; % xgcm

     meangcmbase = squeeze(nanmean(meanbasedata,3));
     for xgwl=2:length(gwlevels),
        meangcmfutallscen = squeeze(meandata(:,:,:,xgwl,:));
        meangcmfutallscen = squeeze(nanmean(squeeze(meangcmfutallscen(:,:,:)),3));
        meandeltagcmallscen = (((meangcmfutallscen./meangcmbase)-1)*100);
        if (sum(sum((~isnan(meandeltagcmallscen))))>0),
           %p = importdata('/home/mmphill2/matlab/matfiles/cmap_percentchange.mat');
           p = importdata('/Users/mmphill2/Documents/MATLAB/cmap_percentchange.mat');
           f = figure; colormap(p);
           acr_pcolormapr5(meandeltagcmallscen.*thiscropmask,lat,lon,[-60 90],[-180 180],[-30 30],['Mean Change in ' thisirr thiscropnice ' Seasonal Precipitation for AllScenario GWL' num2str(gwlevels(xgwl))],['percent change from 1850-1900 ' thisirr thiscropnice ' season mean']);
           %print(f,'-dpng',['/home/mmphill2/GWL' num2str(gwlevels(xgwl)) '_allscen_' thisirr thiscropnice '_MeanPrecipitationChange.png']);
           %print(f,'-depsc',['/home/mmphill2/GWL' num2str(gwlevels(xgwl)) '_allscen_' thisirr thiscropnice '_MeanPrecipitationChange.eps']);
           print(f,'-dpng',['GWL' num2str(gwlevels(xgwl)) '_allscen_' thisirr thiscropnice '_MeanPrecipitationChange.png']);
           print(f,'-depsc',['GWL' num2str(gwlevels(xgwl)) '_allscen_' thisirr thiscropnice '_MeanPrecipitationChange.eps']);

           close all
        end;

        for xscen=1:length(scentypes),
           thisscen = scentypes{xscen};
           thisscennice = scentypesnice{xscen};

           meangcmfut = squeeze(nanmean(squeeze(meandata(:,:,:,xgwl,xscen)),3));
           meandeltagcm = (((meangcmfut./meangcmbase)-1)*100);
           if (sum(sum((~isnan(meandeltagcm))))>0),
              %p = importdata('/home/mmphill2/matlab/matfiles/cmap_percentchange.mat');
              p = importdata('/Users/mmphill2/Documents/MATLAB/cmap_percentchange.mat');

              f = figure; colormap(p);
              acr_pcolormapr5(meandeltagcm.*thiscropmask,lat,lon,[-60 90],[-180 180],[-30 30],['Mean Change in ' thisirr thiscropnice ' Seasonal Precipitation for ' thisscennice ' GWL' num2str(gwlevels(xgwl))],['percent change from 1850-1900 ' thisirr thiscropnice ' season mean']);
              %print(f,'-dpng',['/home/mmphill2/GWL' num2str(gwlevels(xgwl)) '_' thisscennice '_' thisirr thiscropnice '_MeanPrecipitationChange.png']);
              %print(f,'-depsc',['/home/mmphill2/GWL' num2str(gwlevels(xgwl)) '_' thisscennice '_' thisirr thiscropnice '_MeanPrecipitationChange.eps']);
              print(f,'-dpng',['GWL' num2str(gwlevels(xgwl)) '_' thisscennice '_' thisirr thiscropnice '_MeanPrecipitationChange.png']);
              print(f,'-depsc',['GWL' num2str(gwlevels(xgwl)) '_' thisscennice '_' thisirr thiscropnice '_MeanPrecipitationChange.eps']);

              close all
           end;
        end; % xgwl
     end; % xscen

   end;
end;
