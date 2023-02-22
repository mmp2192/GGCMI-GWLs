%		                  mmp_agmip303
%
%  This script creates production change maps from GGCMI
%  Phase 3 simulations for each GCM, global warming level,
%  crop and irrigation type. This version uses pre-processed
%  GGCMI Phase 3 outputs. This is an updated version of
%  mmp_agmip203.m.
%
%  -- Creates a single production change map for each unique
%  combination of GGCM, GCM, crop, irrigation type, SSP
%  scenario, CO2 type and global warming level as available
%  from GGCMI Phase 3.
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
%				                         date:	8/16/2022
%
function mmp_agmip303();
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

cd /css/giss-cig/AgGRID/CropModels/phase3/matfiles/

lat = importdata('AgGRIDlat.mat');
lon = importdata('AgGRIDlon.mat');

cd /css/giss-cig/AgGRID/CropModels/inputs/MIRCA/

allcroparea = ones(360,720,length(croptypes))*NaN;
area = ncread('maize.nc4','sum');
allcroparea(:,:,1) = area';
area = ncread('wheat.nc4','sum');
allcroparea(:,:,2) = area';
area = ncread('soy.nc4','sum');
allcroparea(:,:,3) = area';
area = ncread('rice.nc4','sum');
allcroparea(:,:,4) = area';

cd /css/giss-cig/AgGRID/CropModels/phase3/matfiles/

allyears = [1979:2099];

for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   for xggcm=1:length(modeltypes),
      thisggcm = modeltypes{xggcm};
      thisggcmnice = modeltypesnice{xggcm};
      for xgcm=1:length(gcmtypes),
         thisgcm = gcmtypes{xgcm};
	       thisgcmnice = gcmtypesnice{xgcm};

	       for xscen=1:length(scentypes),
	          thisscen = scentypes{xscen};

            futfile = ['yield_ref_' thisggcm '_' thiscrop '_total_' thisgcm '_' thisscen '_processed.mat'];

		        if (exist(futfile)),
		           futyield = importdata([futfile]);
		           futproduct = futyield*NaN;
	             for yy=1:size(futyield,3),
	                futproduct(:,:,yy) = futyield(:,:,yy).*squeeze(allcroparea(:,:,xcrop));
	             end;
	             clear futyield
		           save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/' thisggcm '_' thisgcm '_total_' thiscrop '_production_' thisscen '_processed.mat'],'futproduct');

		         end; % fut file exist
	        end;  % scenario
       end;  % gcm
   end;  % ggcm
end;  % crop

gwlevels = [0.69 1 1.5 2 2.5 3 3.5 4];

cd /discover/nobackup/projects/giss_ana/users/mmphill2/GWLPaper/

startyy = importdata('GWL_StartYear.mat');
endyy = importdata('GWL_EndYear.mat');

allcropmask = ones(360,720,length(croptypes))*NaN;

cd /css/giss-cig/AgGRID/CropModels/inputs/MIRCA/

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

cd /discover/nobackup/projects/giss_ana/users/mmphill2/GWLPaper/

newcropmasks = importdata('GGCMI_CropMasks_processed.mat');

%%% actual maps of total production change
%%% by GWL compared to GWL0.69

hexcolors = importdata('/home/mmphill2/matlab/matfiles/cmap_productionchange.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;
p(19:end,3) = p(19:end,3)+0.25;
p(18,:) = [0.9529 0.9687 0.8093];

cd /css/giss-cig/AgGRID/CropModels/phase3/matfiles/

for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   thiscropnice = croptypesnice{xcrop};

   %cropmask = squeeze(allcropmask(:,:,xcrop));
   cropmask = squeeze(newcropmasks(:,:,xcrop));

   for xggcm=1:length(modeltypes),
      thisggcm = modeltypes{xggcm};
      thisggcmnice = modeltypesnice{xggcm};
      for xgcm=1:length(gcmtypes),
         thisgcm = gcmtypes{xgcm};
	       thisgcmnice = gcmtypesnice{xgcm};
	       for xscen = 1:length(scentypes),
	          thisscen = scentypes{xscen};
	          thisscennice = scentypesnice{xscen};

	          filename = [thisggcm '_' thisgcm '_total_' thiscrop '_production_' thisscen '_processed.mat'];
		        if (exist(filename)),
		           data = importdata([filename]);
		           yystart = startyy(xgcm,1,(xscen));
		           startidx = find(allyears == yystart);
		           yyend = endyy(xgcm,1,(xscen));
		           endidx = find(allyears == yyend);
		           basedata = squeeze(nanmean(squeeze(data(:,:,startidx:endidx)),3));
		           for gw=2:length(gwlevels),
		              yystart = startyy(xgcm,gw,(xscen));
		              startidx = find(allyears == yystart);
		              yyend = endyy(xgcm,gw,(xscen));
		              endidx = find(allyears == yyend);
			            futdata = squeeze(nanmean(squeeze(data(:,:,startidx:endidx)),3));
			            deltaproduct = (((futdata./basedata)-1)*100);
			            if (sum(sum(~isnan(deltaproduct)))>0),
		 	               f = figure; colormap(p);
			               acr_pcolormapr5(deltaproduct.*cropmask,lat,lon,[-60 90],[-180 180],[-50 50],[thisggcmnice ' ' thisgcmnice ' ' thiscropnice ' ' thisscennice ' GWL:' num2str(gwlevels(gw)) ' (Default)'],['% change in production']);
			               print(f,'-dpng',['/home/mmphill2/figures/GGCMI/' thisggcm '_' thisgcm '_' thiscrop '_' thisscen '_GWL' num2str(gwlevels(gw)) '_default_deltaproduction.png']);
			            end;
			            close all

		           end; % gw level
		        end; % file exists
         end; % scen
      end;  % gcm
   end;  % ggcm
end;  % crop

cd /css/giss-cig/AgGRID/CropModels/phase3/matfiles/

%%%%%%%%%%%%%%%%%%%% ensemble MEDIAN % CHANGE across GGCMs/GCMs
%%%%%%%%%%%%%%%%%%%% for each GWL x crop

hexcolors = importdata('/home/mmphill2/matlab/matfiles/cmap_productionchange.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;
p(19:end,3) = p(19:end,3)+0.25;
p(18,:) = [0.9529 0.9687 0.8093];

for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   thiscropnice = croptypesnice{xcrop};

   %cropmask = squeeze(allcropmask(:,:,xcrop));
   cropmask = squeeze(newcropmasks(:,:,xcrop));

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
   allcropdata = ones(size(lat,1),size(lat,2),length(modeltypes2),length(gcmtypes),length(scentypes),(length(gwlevels)-1))*NaN;
   for xggcm=1:length(modeltypes2),
      thisggcm = modeltypes2{xggcm};
      for xgcm=1:length(gcmtypes),
         thisgcm = gcmtypes{xgcm};
	       thisgcmnice = gcmtypesnice{xgcm};
	       for xscen = 1:length(scentypes),
	          thisscen = scentypes{xscen};
	          thisscennice = scentypesnice{xscen};

            filename = [thisggcm '_' thisgcm '_total_' thiscrop '_production_' thisscen '_processed.mat'];
		        if (exist(filename)),
		           data = importdata([filename]);
		           yystart = startyy(xgcm,1,(xscen));
		           startidx = find(allyears == yystart);
		           yyend = endyy(xgcm,1,(xscen));
		           endidx = find(allyears == yyend);
		           basedata = squeeze(nanmean(squeeze(data(:,:,startidx:endidx)),3));
		           for gw=2:length(gwlevels),
		              yystart = startyy(xgcm,gw,(xscen));
		              startidx = find(allyears == yystart);
		              yyend = endyy(xgcm,gw,(xscen));
		              endidx = find(allyears == yyend);
			            futdata = squeeze(nanmean(squeeze(data(:,:,startidx:endidx)),3));
			            deltaproduct = (((futdata./basedata)-1)*100);
                  allcropdata(:,:,xggcm,xgcm,(xscen),(gw-1)) = deltaproduct;
		           end; % gw level
		        end; % file exists
         end; % scen
      end;  % gcm
   end;  % ggcm

   meancropdata = permute(allcropdata,[1 2 5 6 3 4]);
   clear allcropdata
   meancropdata = squeeze(meancropdata(:,:,:,:,:));
   meancropdata = squeeze(nanmedian(meancropdata,5));
   for xscen = 1:length(scentypes),
      thisscen = scentypes{xscen};
      thisscennice = scentypesnice{xscen};
      for gw=2:length(gwlevels),
         thisgw = gwlevels(gw);
         deltaproduct = squeeze(meancropdata(:,:,(xscen),(gw-1)));

         if (sum(sum(~isnan(deltaproduct)))>0),
            f = figure; colormap(p);
            acr_pcolormapr5(deltaproduct.*cropmask,lat,lon,[-60 90],[-180 180],[-50 50],['Median Across GGCMs/GCMs ' thiscropnice ' ' thisscennice ' GWL:' num2str(gwlevels(gw)) ' (Default)'],['% change in production']);
            print(f,'-dpng',['/home/mmphill2/median_' thiscrop '_' thisscen '_GWL' num2str(gwlevels(gw)) '_default_deltaproduction.png']);
            print(f,'-depsc',['/home/mmphill2/median_' thiscrop '_' thisscen '_GWL' num2str(gwlevels(gw)) '_default_deltaproduction.eps']);
            close all
         end;
      end; % gwl
   end; % scen
end;  % crop

%%%%%%%%%%%%%%%%%%%% ensemble MEDIAN % CHANGE across GGCMs/GCMs
%%%%%%%%%%%%%%%%%%%% for each GWL x crop
%%%%%%%%%%%%%%%%%%%% UNIFIED SCENARIOS

hexcolors = importdata('/home/mmphill2/matlab/matfiles/cmap_productionchange.mat');
p = ones(length(hexcolors),3)*NaN;
for rgb=1:size(p,1),
   p(rgb,:) = hex2rgb(hexcolors(rgb));
end;
p(19:end,3) = p(19:end,3)+0.25;
p(18,:) = [0.9529 0.9687 0.8093];

for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   thiscropnice = croptypesnice{xcrop};

   %cropmask = squeeze(allcropmask(:,:,xcrop));
   cropmask = squeeze(newcropmasks(:,:,xcrop));

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
   allcropdata = ones(size(lat,1),size(lat,2),length(modeltypes2),length(gcmtypes),length(scentypes),(length(gwlevels)-1))*NaN;

   for xggcm=1:length(modeltypes2),
      thisggcm = modeltypes2{xggcm};
      for xgcm=1:length(gcmtypes),
         thisgcm = gcmtypes{xgcm};
	       thisgcmnice = gcmtypesnice{xgcm};
	       for xscen = 1:length(scentypes),
	          thisscen = scentypes{xscen};
	          thisscennice = scentypesnice{xscen};

            filename = [thisggcm '_' thisgcm '_total_' thiscrop '_production_' thisscen '_processed.mat'];
		        if (exist(filename)),
		           data = importdata([filename]);
		           yystart = startyy(xgcm,1,(xscen));
		           startidx = find(allyears == yystart);
		           yyend = endyy(xgcm,1,(xscen));
		           endidx = find(allyears == yyend);
		           basedata = squeeze(nanmean(squeeze(data(:,:,startidx:endidx)),3));
		           for gw=2:length(gwlevels),
		              yystart = startyy(xgcm,gw,(xscen));
		              startidx = find(allyears == yystart);
		              yyend = endyy(xgcm,gw,(xscen));
		              endidx = find(allyears == yyend);
			            futdata = squeeze(nanmean(squeeze(data(:,:,startidx:endidx)),3));
			            deltaproduct = (((futdata./basedata)-1)*100);
                  allcropdata(:,:,xggcm,xgcm,(xscen),(gw-1)) = deltaproduct;
		           end; % gw level
		        end; % file exists
         end; % scen
      end;  % gcm
   end;  % ggcm

   meancropdata = permute(allcropdata,[1 2 6 3 4 5]);
   clear allcropdata
   meancropdata = squeeze(meancropdata(:,:,:,:));
   meancropdata = squeeze(nanmedian(meancropdata,4));
   for gw=2:length(gwlevels),
      thisgw = gwlevels(gw);
      deltaproduct = squeeze(meancropdata(:,:,(gw-1)));
      if (sum(sum(~isnan(deltaproduct)))>0),
         f = figure; colormap(p);
         acr_pcolormapr5(deltaproduct.*cropmask,lat,lon,[-60 90],[-180 180],[-50 50],['Median Across GGCMs/GCMs ' thiscropnice ' AllScenario GWL:' num2str(gwlevels(gw)) ' (Default)'],['% change in production']);
         print(f,'-dpng',['/home/mmphill2/median_' thiscrop '_allscen_GWL' num2str(gwlevels(gw)) '_default_deltaproduction.png']);
         print(f,'-depsc',['/home/mmphill2/median_' thiscrop '_allscen_GWL' num2str(gwlevels(gw)) '_default_deltaproduction.eps']);
         close all
      end;
   end; % gwl
end;  % crop
