%		       mmp_agmip302
%
%  This script calculates absolute yields from GGCMI
%  Phase 3 pre-processed outputs (adjusted harvest dates,
%  combined crop seasons, combined rainfed and irrigated
%  production, etc.) and aggregates to national and global
%  production levels. This is an updated version of
%  mmp_agmip202.m.
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
%				             date:  8/15/2022
%
function mmp_agmip302();
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

cd /discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/

allcroparea = ones(360,720,length(croptypes))*NaN;
area = ncread('maize.nc4','sum');
allcroparea(:,:,1) = area';
area = ncread('wheat.nc4','sum');
allcroparea(:,:,2) = area';
area = ncread('soy.nc4','sum');
allcroparea(:,:,3) = area';
area = ncread('rice.nc4','sum');
allcroparea(:,:,4) = area';

lat = importdata('AgGRIDlat.mat');
lon = importdata('AgGRIDlon.mat');
gadm0countrymap = importdata(['gadm0countrymap.mat']);

gwlevels = [0.69 1 1.5 2 2.5 3 3.5 4];
startyy = ones(length(gcmtypes),length(gwlevels),length(scentypes))*NaN;
endyy = ones(length(gcmtypes),length(gwlevels),length(scentypes))*NaN;

%% annual global mean temperature time series

allyears = [1850:2100];
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

%% 20-year centered average
%% to replicate Hauser et al. GWLs
%% and calculate our own GWL2.5

temp20 = annualgmtanom*NaN;
for xgcm=1:length(gcmtypes),
   for xscen=1:length(scentypes),
      for yy=11:(length(fullyears)-9),
         temp20(xgcm,xscen,yy) = squeeze(nanmean(squeeze(annualgmtanom(xgcm,xscen,((yy-10):(yy+9))))));
      end;
   end;
end;

for xgcm=1:length(gcmtypes),
   for xscen=1:length(scentypes),
      thistemp20 = squeeze(temp20(xgcm,xscen,:));
      for xgwl=1:length(gwlevels),
         thisgwl = gwlevels(xgwl);
         centralyear = find(thistemp20>thisgwl,1,'first');
         if (length(centralyear)>0)
         startyy(xgcm,xgwl,xscen) = allyears(centralyear-10);
         endyy(xgcm,xgwl,xscen) = startyy(xgcm,xgwl,xscen)+19;
         end;
      end;
   end;
end;

save('GWL_StartYear.mat','startyy');
save('GWL_EndYear.mat','endyy');

allyears = [1979:2099];

refyield = ncread('GGCMI_observational_reference_yield_map.nc4','yield');

%% Loop through crops
for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   areamap = squeeze(allcroparea(:,:,xcrop));
      for xmodel = 1:length(modeltypes),
         thismodel = modeltypes{xmodel};
	       for xscen = 1:length(scentypes),
            thisscen = scentypes{xscen};
            globalproduct = ones(length(gcmtypes),length(allyears))*NaN;

            %% Loop through gcms
            for xgcm = 1:length(gcmtypes),
               thisgcm = gcmtypes{xgcm};
               filename = ['deltayield_' thismodel '_' thiscrop '_total_' thisgcm '_' thisscen '_processed.mat'];

               if(exist(filename)),
                  yield = importdata([filename]);

                  %%% calculate absolute yields
                  %%% using GGCMI observational yield reference map

                  absyield = yield*NaN;
                  for xyear=1:size(yield,3),
                     %absyield(:,:,xyear) = (((squeeze(yield(:,:,xyear)))/100)+1).*(squeeze(refyield(:,:,xcrop,3))'); %'
                     absyield(:,:,xyear) = (squeeze(refyield(:,:,xcrop,3))') + (((squeeze(yield(:,:,xyear)))/100).*(squeeze(refyield(:,:,xcrop,3))'));
                  end;
                  clear yield

                  save(['yield_ref_' thismodel '_' thiscrop '_total_' thisgcm '_' thisscen '_processed.mat'],'absyield');

                  %%%%%%% global aggregation
                  for xyear=1:size(absyield,3),
                     product = squeeze(absyield(:,:,xyear)).*areamap;
                     globalproduct(xgcm,xyear) = nansum(nansum(product));
                  end;  %% xyear

               end; % exist
            end; % gcm

            if (sum(sum(~isnan(globalproduct)))>0),
	             save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/global_production_' thismodel '_' thiscrop '_total_' thisscen '_processed.mat'],'globalproduct');
	          end;
         end;  % scen
      end;  % model
end;  % crop

%%%%%%%% country level aggregation

%% Loop through crops
for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
      areamap = squeeze(allcroparea(:,:,xcrop));
      %% Loop through models
      for xmodel = 1:length(modeltypes),
         thismodel = modeltypes{xmodel};
         %% Loop through gcms
         for xgcm = 1:length(gcmtypes),
            thisgcm = gcmtypes{xgcm};
	          %%% Loop through scenarios
	          for xscen = 1:length(scentypes),
               thisscen = scentypes{xscen};

		           filename = ['yield_ref_' thismodel '_' thiscrop '_total_' thisgcm '_' thisscen '_processed.mat'];

               if(exist(filename)),
		              yield = importdata([filename]);
		              countryproduct = ones(max(max(gadm0countrymap)),size(yield,3))*NaN;
		              countryarea = ones(max(max(gadm0countrymap)),size(yield,3))*NaN;
		              for yy=1:size(yield,3),
		                 yieldmask = ~isnan(squeeze(yield(:,:,yy)));
			               product = squeeze(yield(:,:,yy)).*areamap;
			               for thiscountry=1:max(max(gadm0countrymap)),
			                  areamask = areamap*NaN;
			                  productmask = product*NaN;
			                  if(sum(sum(((gadm0countrymap==thiscountry)&(yieldmask==1))))>0),
			                     areamask(gadm0countrymap==thiscountry) = areamap(gadm0countrymap==thiscountry);
			                     productmask(gadm0countrymap==thiscountry) = product(gadm0countrymap==thiscountry);

			                     countryproduct(thiscountry,yy) = nansum(nansum(productmask));
			                     countryarea(thiscountry,yy) = nansum(nansum(areamask));
			                  end;  %% yieldmask
			              end; %% country
		           end;  %% yy
               save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/country_production_' thismodel '_' thiscrop '_total_' thisgcm '_' thisscen '_processed.mat'],'countryproduct');

            end;     %% if filename exists
         end;        %% scenario
       end;           %% gcm
   end;                 %% models

   save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/country_area_total_' thiscrop '_MIRCA.mat'],'countryarea');
end;                       %% crops
