%		      mmp_agmip304
%
%  This script creates line plots of global warming levels
%  against global production from GGCMI-3 simulations. This
%  version uses pre-processed GGCMI Phase 3 outputs. This
%  is an updated version of mmp_agmip204.m.
%
%  -- Creates a line plot of global warming level (x-axis)
%  versus global production anomaly showing all GCMs and SSP
%  scenarios for each GGCM and crop.
%
%  -- Creates a line plot of global warming level (x-axis)
%  versus global production anomaly showing all GCMs for
%  SSP585 only for the median across GGCMs for each crop.
%
%  -- Creates a line plot of global warming level (x-axis)
%  versus global production anomaly showing all GCMs and SSP
%  scenarios for the median across GGCMs for each crop.
%
%  -- Creates a line plot of global warming level (x-axis)
%  versus global production anomaly showing all GCMs and SSP
%  scenarios for the median across only the GGCMs that ran
%  SSP370 for each crop.
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
function mmp_agmip304();
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

cd /discover/nobackup/projects/giss_ana/users/mmphill2/GWLPaper/

allyears = [1979:2099];

gwlevels = [0.69 1 1.5 2 2.5 3 3.5 4];

startyy = importdata('GWL_StartYear.mat');
endyy = importdata('GWL_EndYear.mat');

for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};

   for xmodel=1:length(modeltypes),
      thisggcm = modeltypes{xmodel};

      allscendata = ones(length(gcmtypes),length(allyears),length(scentypes))*NaN;

      for xscen=1:length(scentypes),
         thisscen = scentypes{xscen};

         filename = ['global_production_' thisggcm '_' thiscrop '_total_' thisscen '_processed.mat'];
         if(exist(filename)),
            data = importdata(['global_production_' thisggcm '_' thiscrop '_total_' thisscen '_processed.mat']);
            allscendata(:,:,xscen) = data;
	          clear data
	       end; % exist
       end; % scen

       if (sum(sum(sum(~isnan(allscendata))))>0),
          save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/total_global_production_' thisggcm '_' thiscrop '_processed.mat'],'allscendata');
       end;
   end; % model
end; % crop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% line plots of GWL x Global Production

cd /css/giss-cig/AgGRID/CropModels/phase3/matfiles/

for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   thiscropnice = croptypesnice{xcrop};
   for xmodel=1:length(modeltypes),
      thisggcm = modeltypes{xmodel};
      thisggcmnice = modeltypesnice{xmodel};

      scendata = ['total_global_production_' thisggcm '_' thiscrop '_processed.mat'];
      if(exist(scendata)),
         allscendata = importdata(['total_global_production_' thisggcm '_' thiscrop '_processed.mat']);
         allgwlevel = ones(length(gcmtypes),length(gwlevels),length(scentypes))*NaN;
         for gw=1:length(gwlevels),
            for xgcm=1:length(gcmtypes),
	             for xscen=1:length(scentypes),
                  gwstart = startyy(xgcm,gw,xscen);
	                startidx = find(allyears == gwstart);
	                gwend = endyy(xgcm,gw,xscen);
	                endidx = find(allyears == gwend);
	                allgwlevel(xgcm,gw,xscen) = squeeze(nanmean(squeeze(allscendata(xgcm,startidx:endidx,xscen))));
               end;
	          end;
         end;
         allgwlevelanom = ones(length(gcmtypes),length(gwlevels),length(scentypes))*NaN;
         for gw=1:length(gwlevels),
            for xgcm=1:length(gcmtypes),
	             for xscen=1:length(scentypes),
	                base = squeeze(allgwlevel(xgcm,1,xscen));
	                allgwlevelanom(xgcm,gw,xscen) = (((squeeze(allgwlevel(xgcm,gw,xscen))/base)-1)*100);
               end;
	          end;
         end;

         gcmcols = ones(5,3)*NaN;
         gcmcols(1,:) = [0 0 128];
         gcmcols(2,:) = [245 130 48];
         gcmcols(3,:) = [0 130 200];
         gcmcols(4,:) = [202 202 43];
         gcmcols(5,:) = [128 0 0];
         gcmcols = gcmcols/255;
         clear h
         f = figure;
         h(1) = plot(gwlevels,squeeze(allgwlevelanom(1,:,1)),'-p','color',squeeze(gcmcols(1,:)),'markerfacecolor',squeeze(gcmcols(1,:)),'markersize',4,'linewidth',2); hold on;
         plot(gwlevels,squeeze(allgwlevelanom(1,:,2)),'--p','color',squeeze(gcmcols(1,:)),'markerfacecolor',squeeze(gcmcols(1,:)),'markersize',4,'linewidth',2); hold on;
         plot(gwlevels,squeeze(allgwlevelanom(1,:,3)),':p','color',squeeze(gcmcols(1,:)),'markerfacecolor',squeeze(gcmcols(1,:)),'markersize',4,'linewidth',2); hold on;
         h(2) = plot(gwlevels,squeeze(allgwlevelanom(2,:,1)),'-p','color',squeeze(gcmcols(2,:)),'markerfacecolor',squeeze(gcmcols(2,:)),'markersize',4,'linewidth',2); hold on;
         plot(gwlevels,squeeze(allgwlevelanom(2,:,2)),'--p','color',squeeze(gcmcols(2,:)),'markerfacecolor',squeeze(gcmcols(2,:)),'markersize',4,'linewidth',2); hold on;
         plot(gwlevels,squeeze(allgwlevelanom(2,:,3)),':p','color',squeeze(gcmcols(2,:)),'markerfacecolor',squeeze(gcmcols(2,:)),'markersize',4,'linewidth',2); hold on;
         h(3) = plot(gwlevels,squeeze(allgwlevelanom(3,:,1)),'-p','color',squeeze(gcmcols(3,:)),'markerfacecolor',squeeze(gcmcols(3,:)),'markersize',4,'linewidth',2); hold on;
         plot(gwlevels,squeeze(allgwlevelanom(3,:,2)),'--p','color',squeeze(gcmcols(3,:)),'markerfacecolor',squeeze(gcmcols(3,:)),'markersize',4,'linewidth',2); hold on;
         plot(gwlevels,squeeze(allgwlevelanom(3,:,3)),':p','color',squeeze(gcmcols(3,:)),'markerfacecolor',squeeze(gcmcols(3,:)),'markersize',4,'linewidth',2); hold on;
         h(4) = plot(gwlevels,squeeze(allgwlevelanom(4,:,1)),'-p','color',squeeze(gcmcols(4,:)),'markerfacecolor',squeeze(gcmcols(4,:)),'markersize',4,'linewidth',2); hold on;
         plot(gwlevels,squeeze(allgwlevelanom(4,:,2)),'--p','color',squeeze(gcmcols(4,:)),'markerfacecolor',squeeze(gcmcols(4,:)),'markersize',4,'linewidth',2); hold on;
         plot(gwlevels,squeeze(allgwlevelanom(4,:,3)),':p','color',squeeze(gcmcols(4,:)),'markerfacecolor',squeeze(gcmcols(4,:)),'markersize',4,'linewidth',2); hold on;
         h(5) = plot(gwlevels,squeeze(allgwlevelanom(5,:,1)),'-p','color',squeeze(gcmcols(5,:)),'markerfacecolor',squeeze(gcmcols(5,:)),'markersize',4,'linewidth',2); hold on;
         plot(gwlevels,squeeze(allgwlevelanom(5,:,2)),'--p','color',squeeze(gcmcols(5,:)),'markerfacecolor',squeeze(gcmcols(5,:)),'markersize',4,'linewidth',2); hold on;
         plot(gwlevels,squeeze(allgwlevelanom(5,:,3)),':p','color',squeeze(gcmcols(5,:)),'markerfacecolor',squeeze(gcmcols(5,:)),'markersize',4,'linewidth',2); hold on;
         h(6) = plot([1 2],[NaN NaN],'-p','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',4,'linewidth',2);
         h(7) = plot([1 2],[NaN NaN],'--p','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',4,'linewidth',2);
         h(8) = plot([1 2],[NaN NaN],':p','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',4,'linewidth',2);
         xlim([0.69 4])
         xlabel(['Global Warming Level']);
         ylabel(['% Change in Global ' thiscropnice ' Production']);
         title({[thiscropnice ' Global Production Change by'], ['Global Warming Level from ' thisggcmnice]});
         legendnames = {'GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL','SSP1-2.6','SSP3-7.0','SSP5-8.5'};
         legend(h(:),legendnames,'location','eastoutside');
         print(f,'-dpng',['/home/mmphill2/figures/GGCMI/' thisggcm '_' thiscrop '_global_deltaproduction_gwl_defaultco2.png']);

         close all
      end;
   end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SSP585 only

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

   allgwlevelanom = ones(length(modeltypes2),length(gcmtypes),length(gwlevels),length(scentypes))*NaN;

   for xmodel=1:length(modeltypes2),
      thisggcm = modeltypes2{xmodel};

      scendata = ['total_global_production_' thisggcm '_' thiscrop '_processed.mat'];
      if(exist(scendata)),
         allscendata = importdata(['total_global_production_' thisggcm '_' thiscrop '_processed.mat']);
         allgwlevel = ones(length(gcmtypes),length(gwlevels),length(scentypes))*NaN;
         for gw=1:length(gwlevels),
            for xgcm=1:length(gcmtypes),
	             for xscen=1:length(scentypes),
                  gwstart = startyy(xgcm,gw,xscen);
	                startidx = find(allyears == gwstart);
	                gwend = endyy(xgcm,gw,xscen);
	                endidx = find(allyears == gwend);
	                allgwlevel(xgcm,gw,xscen) = squeeze(nanmean(squeeze(allscendata(xgcm,startidx:endidx,xscen))));
               end;
	          end;
         end;

         for gw=1:length(gwlevels),
            for xgcm=1:length(gcmtypes),
	             for xscen=1:length(scentypes),
	                base = squeeze(allgwlevel(xgcm,1,xscen));
	                allgwlevelanom(xmodel,xgcm,gw,xscen) = (((squeeze(allgwlevel(xgcm,gw,xscen))/base)-1)*100);
               end;
	          end;
         end;
      end;
   end;
   mediangwlevel = squeeze(nanmedian(squeeze(allgwlevelanom(:,:,:,3))));

   gcmcols = ones(5,3)*NaN;
   gcmcols(1,:) = [0 0 128];
   gcmcols(2,:) = [245 130 48];
   gcmcols(3,:) = [0 130 200];
   gcmcols(4,:) = [202 202 43];
   gcmcols(5,:) = [128 0 0];
   gcmcols = gcmcols/255;

   gcmgwlmax = ones(length(gcmtypes),2)*NaN;
   for xgcm=1:length(gcmtypes),
      maxgwlidx = sum(~isnan(squeeze(mediangwlevel(xgcm,:))));
      gcmgwlmax(xgcm,1) = gwlevels(maxgwlidx);
      gcmgwlmax(xgcm,2) = maxgwlidx;
   end;

   clear h
   f = figure;
   h(1) = plot(gwlevels,squeeze(mediangwlevel(1,:)),':p','color',squeeze(gcmcols(1,:)),'markerfacecolor',squeeze(gcmcols(1,:)),'markersize',4,'linewidth',2); hold on;
   h(2) = plot(gwlevels,squeeze(mediangwlevel(2,:)),':p','color',squeeze(gcmcols(2,:)),'markerfacecolor',squeeze(gcmcols(2,:)),'markersize',4,'linewidth',2); hold on;
   h(3) = plot(gwlevels,squeeze(mediangwlevel(3,:)),':p','color',squeeze(gcmcols(3,:)),'markerfacecolor',squeeze(gcmcols(3,:)),'markersize',4,'linewidth',2); hold on;
   h(4) = plot(gwlevels,squeeze(mediangwlevel(4,:)),':p','color',squeeze(gcmcols(4,:)),'markerfacecolor',squeeze(gcmcols(4,:)),'markersize',4,'linewidth',2); hold on;
   h(5) = plot(gwlevels,squeeze(mediangwlevel(5,:)),':p','color',squeeze(gcmcols(5,:)),'markerfacecolor',squeeze(gcmcols(5,:)),'markersize',4,'linewidth',2); hold on;
   %for xgcm=1:length(gcmtypes),
   %plot(gcmgwlmax(xgcm,1),mediangwlevel(xgcm,gcmgwlmax(xgcm,2)),'.','color',squeeze(gcmcols(xgcm,:)),'markerfacecolor',squeeze(gcmcols(xgcm,:)),'markersize',40); hold on;
   %end;
   xlim([0.69 4])
   xlabel(['Global Warming Level']);
   ylabel(['% Change in Global ' thiscropnice ' Production']);
   title({[thiscropnice ' Global Production Change'], ['by GWL for SSP5-8.5 (Median GGCM)']});
   legendnames = {'GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL'};
   if (strcmp(thiscrop,'mai')),
      legend(h(:),legendnames,'location','southwest');
   end;
   print(f,'-dpng',['/home/mmphill2/figures/GGCMI/median_' thiscrop '_global_deltaproduction_ssp585_gwl_defaultco2.png']);
   print(f,'-depsc',['/home/mmphill2/figures/GGCMI/median_' thiscrop '_global_deltaproduction_ssp585_gwl_defaultco2.eps']);

   close all

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% allscen -- separated by GCM

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

   allgwlevelanom = ones(length(modeltypes2),length(gcmtypes),length(gwlevels),length(scentypes))*NaN;

   for xmodel=1:length(modeltypes2),
      thisggcm = modeltypes2{xmodel};

      scendata = ['total_global_production_' thisggcm '_' thiscrop '_processed.mat'];
      if(exist(scendata)),

         allscendata = importdata(['total_global_production_' thisggcm '_' thiscrop '_processed.mat']);
         allgwlevel = ones(length(gcmtypes),length(gwlevels),length(scentypes))*NaN;
         for gw=1:length(gwlevels),
            for xgcm=1:length(gcmtypes),
	             for xscen=1:length(scentypes),
                  gwstart = startyy(xgcm,gw,xscen);
	                startidx = find(allyears == gwstart);
	                gwend = endyy(xgcm,gw,xscen);
	                endidx = find(allyears == gwend);
	                allgwlevel(xgcm,gw,xscen) = squeeze(nanmean(squeeze(allscendata(xgcm,startidx:endidx,xscen))));
              end;
	         end;
        end;

        for gw=1:length(gwlevels),
           for xgcm=1:length(gcmtypes),
	            for xscen=1:length(scentypes),
	               base = squeeze(allgwlevel(xgcm,1,xscen));
	               allgwlevelanom(xmodel,xgcm,gw,xscen) = (((squeeze(allgwlevel(xgcm,gw,xscen))/base)-1)*100);
               end;
	          end;
         end;
      end;
   end;
   allscengwlevelanom = permute(allgwlevelanom,[2 3 4 1]);
   allscengwlevelanom = squeeze(allscengwlevelanom(:,:,:));

   mediangwlevel = squeeze(nanmedian(allscengwlevelanom,3));

   gcmcols = ones(5,3)*NaN;
   gcmcols(1,:) = [0 0 128];
   gcmcols(2,:) = [245 130 48];
   gcmcols(3,:) = [0 130 200];
   gcmcols(4,:) = [202 202 43];
   gcmcols(5,:) = [128 0 0];
   gcmcols = gcmcols/255;

   gcmgwlmax = ones(length(gcmtypes),2)*NaN;
   for xgcm=1:length(gcmtypes),
      maxgwlidx = sum(~isnan(squeeze(mediangwlevel(xgcm,:))));
      gcmgwlmax(xgcm,1) = gwlevels(maxgwlidx);
      gcmgwlmax(xgcm,2) = maxgwlidx;
   end;

   clear h
   f = figure;
   h(1) = plot(gwlevels,squeeze(mediangwlevel(1,:)),':p','color',squeeze(gcmcols(1,:)),'markerfacecolor',squeeze(gcmcols(1,:)),'markersize',4,'linewidth',2); hold on;
   h(2) = plot(gwlevels,squeeze(mediangwlevel(2,:)),':p','color',squeeze(gcmcols(2,:)),'markerfacecolor',squeeze(gcmcols(2,:)),'markersize',4,'linewidth',2); hold on;
   h(3) = plot(gwlevels,squeeze(mediangwlevel(3,:)),':p','color',squeeze(gcmcols(3,:)),'markerfacecolor',squeeze(gcmcols(3,:)),'markersize',4,'linewidth',2); hold on;
   h(4) = plot(gwlevels,squeeze(mediangwlevel(4,:)),':p','color',squeeze(gcmcols(4,:)),'markerfacecolor',squeeze(gcmcols(4,:)),'markersize',4,'linewidth',2); hold on;
   h(5) = plot(gwlevels,squeeze(mediangwlevel(5,:)),':p','color',squeeze(gcmcols(5,:)),'markerfacecolor',squeeze(gcmcols(5,:)),'markersize',4,'linewidth',2); hold on;
   %for xgcm=1:length(gcmtypes),
   %plot(gcmgwlmax(xgcm,1),mediangwlevel(xgcm,gcmgwlmax(xgcm,2)),'.','color',squeeze(gcmcols(xgcm,:)),'markerfacecolor',squeeze(gcmcols(xgcm,:)),'markersize',40); hold on;
   %end;
   xlim([0.69 4])
   xlabel(['Global Warming Level']);
   ylabel(['% Change in Global ' thiscropnice ' Production']);
   title({[thiscropnice ' Global Production Change'], ['by GWL for ALLSCEN (Median GGCM/SSP)']});
   legendnames = {'GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL'};
   if (strcmp(thiscrop,'mai')),
      legend(h(:),legendnames,'location','southwest');
   end;
   print(f,'-dpng',['/home/mmphill2/figures/GGCMI/median_' thiscrop '_global_deltaproduction_allscen_byGCM_gwl_defaultco2.png']);
   print(f,'-depsc',['/home/mmphill2/figures/GGCMI/median_' thiscrop '_global_deltaproduction_allscen_byGCM_gwl_defaultco2.eps']);
   close all
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% allscen -- separated by SSP-RCP
%%%% INCLUSIVE version (more GGCMs)

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

   allgwlevelanom = ones(length(modeltypes2),length(gcmtypes),length(gwlevels),length(scentypes))*NaN;

   for xmodel=1:length(modeltypes2),
      thisggcm = modeltypes2{xmodel};

      scendata = ['total_global_production_' thisggcm '_' thiscrop '_processed.mat'];
      if(exist(scendata)),
         allscendata = importdata(['total_global_production_' thisggcm '_' thiscrop '_processed.mat']);
         allgwlevel = ones(length(gcmtypes),length(gwlevels),length(scentypes))*NaN;
         for gw=1:length(gwlevels),
            for xgcm=1:length(gcmtypes),
	             for xscen=1:length(scentypes),
                  gwstart = startyy(xgcm,gw,xscen);
	                startidx = find(allyears == gwstart);
	                gwend = endyy(xgcm,gw,xscen);
	                endidx = find(allyears == gwend);
	                allgwlevel(xgcm,gw,xscen) = squeeze(nanmean(squeeze(allscendata(xgcm,startidx:endidx,xscen))));
              end;
	         end;
        end;

        for gw=1:length(gwlevels),
           for xgcm=1:length(gcmtypes),
	            for xscen=1:length(scentypes),
	               base = squeeze(allgwlevel(xgcm,1,xscen));
	               allgwlevelanom(xmodel,xgcm,gw,xscen) = (((squeeze(allgwlevel(xgcm,gw,xscen))/base)-1)*100);
               end;
	          end;
         end;
      end;
   end;
   mediangwlevel = permute(allgwlevelanom,[3 4 1 2]);
   mediangwlevel = squeeze(nanmedian(squeeze(mediangwlevel(:,:,:)),3));

   gcmcols = ones(5,3)*NaN;
   gcmcols(1,:) = [0 0 128];
   gcmcols(2,:) = [245 130 48];
   gcmcols(3,:) = [0 130 200];
   gcmcols(4,:) = [202 202 43];
   gcmcols(5,:) = [128 0 0];
   gcmcols = gcmcols/255;

   scengwlmax = ones(length(scentypes),2)*NaN;
   for xscen=1:length(scentypes),
      maxgwlidx = sum(~isnan(squeeze(mediangwlevel(:,xscen))));
      scengwlmax(xscen,1) = gwlevels(maxgwlidx);
      scengwlmax(xscen,2) = maxgwlidx;
   end;

   clear h
   f = figure;
   h(1) = plot(gwlevels,squeeze(mediangwlevel(:,1)),'-p','color','k','markerfacecolor','k','markersize',4,'linewidth',2); hold on;
   h(2) = plot(gwlevels,squeeze(mediangwlevel(:,2)),'--p','color','k','markerfacecolor','k','markersize',4,'linewidth',2); hold on;
   h(3) = plot(gwlevels,squeeze(mediangwlevel(:,3)),':p','color','k','markerfacecolor','k','markersize',4,'linewidth',2); hold on;
   plot(1,squeeze(mediangwlevel(2,1)),'p','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',14); hold on;
   plot(3,squeeze(mediangwlevel(6,2)),'p','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',14); hold on;
   plot(3.5,squeeze(mediangwlevel(7,3)),'p','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',14); hold on;
   %for xscen=1:length(scentypes),
   %plot(gscengwlmax(xscen,1),mediangwlevel(scengwlmax(xscen,2),xscen),'p','color','k','markerfacecolor','k','markersize',6); hold on;
   %end;
   xlim([0.69 4])
   xlabel(['Global Warming Level']);
   ylabel(['% Change in Global ' thiscropnice ' Production']);
   title({[thiscropnice ' Global Production Change by GWL'], ['(Median Across GGCMs/ESMs)']});
   legendnames = {'SSP1-2.6','SSP3-7.0','SSP5-8.5'};
   if (strcmp(thiscrop,'mai')),
      legend(h(:),legendnames,'location','southwest');
   end;
   print(f,'-dpng',['/home/mmphill2/median_' thiscrop '_global_deltaproduction_allscen_bySSP_gwl_defaultco2.png']);
   print(f,'-depsc',['/home/mmphill2/median_' thiscrop '_global_deltaproduction_allscen_bySSP_gwl_defaultco2.eps']);
   close all

end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% allscen -- separated by SSP-RCP
%%%% CONSISTENT version (fewer GGCMs)

for xcrop=1:length(croptypes),
   thiscrop = croptypes{xcrop};
   thiscropnice = croptypesnice{xcrop};

   if (strcmp(thiscrop,'mai')),
      modeltypes2 = {'crover','cygma1p74','epic-iiasa','isam','ldndc','lpjml','pepic','promet'};
   end;
   if (strcmp(thiscrop,'ric')),
      modeltypes2 = {'cygma1p74','epic-iiasa','isam','ldndc','lpjml','pepic','promet'};
   end;
   if (strcmp(thiscrop,'soy')),
      modeltypes2 = {'crover','cygma1p74','epic-iiasa','isam','ldndc','lpjml','pepic','promet'};
   end;
   if (strcmp(thiscrop,'whe')),
      modeltypes2 = {'crover','epic-iiasa','isam','ldndc','lpjml','pepic','promet'};
   end;

   allgwlevelanom = ones(length(modeltypes2),length(gcmtypes),length(gwlevels),length(scentypes))*NaN;

   for xmodel=1:length(modeltypes2),
      thisggcm = modeltypes2{xmodel};

      scendata = ['total_global_production_' thisggcm '_' thiscrop '_processed.mat'];
      if(exist(scendata)),
         allscendata = importdata(['total_global_production_' thisggcm '_' thiscrop '_processed.mat']);
         allgwlevel = ones(length(gcmtypes),length(gwlevels),length(scentypes))*NaN;
         for gw=1:length(gwlevels),
            for xgcm=1:length(gcmtypes),
	             for xscen=1:length(scentypes),
                  gwstart = startyy(xgcm,gw,xscen);
	                startidx = find(allyears == gwstart);
	                gwend = endyy(xgcm,gw,xscen);
	                endidx = find(allyears == gwend);
	                allgwlevel(xgcm,gw,xscen) = squeeze(nanmean(squeeze(allscendata(xgcm,startidx:endidx,xscen))));
              end;
	         end;
        end;

        for gw=1:length(gwlevels),
           for xgcm=1:length(gcmtypes),
	            for xscen=1:length(scentypes),
	               base = squeeze(allgwlevel(xgcm,1,xscen));
	               allgwlevelanom(xmodel,xgcm,gw,xscen) = (((squeeze(allgwlevel(xgcm,gw,xscen))/base)-1)*100);
               end;
	          end;
         end;
      end;
   end;
   mediangwlevel = permute(allgwlevelanom,[3 4 1 2]);
   mediangwlevel = squeeze(nanmedian(squeeze(mediangwlevel(:,:,:)),3));

   gcmcols = ones(5,3)*NaN;
   gcmcols(1,:) = [0 0 128];
   gcmcols(2,:) = [245 130 48];
   gcmcols(3,:) = [0 130 200];
   gcmcols(4,:) = [202 202 43];
   gcmcols(5,:) = [128 0 0];
   gcmcols = gcmcols/255;

   scengwlmax = ones(length(scentypes),2)*NaN;
   for xscen=1:length(scentypes),
      maxgwlidx = sum(~isnan(squeeze(mediangwlevel(:,xscen))));
      scengwlmax(xscen,1) = gwlevels(maxgwlidx);
      scengwlmax(xscen,2) = maxgwlidx;
   end;

   clear h
   f = figure;
   h(1) = plot(gwlevels,squeeze(mediangwlevel(:,1)),'-p','color','k','markerfacecolor','k','markersize',4,'linewidth',2); hold on;
   h(2) = plot(gwlevels,squeeze(mediangwlevel(:,2)),'--p','color','k','markerfacecolor','k','markersize',4,'linewidth',2); hold on;
   h(3) = plot(gwlevels,squeeze(mediangwlevel(:,3)),':p','color','k','markerfacecolor','k','markersize',4,'linewidth',2); hold on;
   plot(1,squeeze(mediangwlevel(2,1)),'p','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',14); hold on;
   plot(3,squeeze(mediangwlevel(6,2)),'p','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',14); hold on;
   plot(3.5,squeeze(mediangwlevel(7,3)),'p','color',[0.5 0.5 0.5],'markerfacecolor',[0.5 0.5 0.5],'markersize',14); hold on;
   %for xscen=1:length(scentypes),
   %plot(gscengwlmax(xscen,1),mediangwlevel(scengwlmax(xscen,2),xscen),'p','color','k','markerfacecolor','k','markersize',6); hold on;
   %end;
   xlim([0.69 4])
   xlabel(['Global Warming Level']);
   ylabel(['% Change in Global ' thiscropnice ' Production']);
   title({[thiscropnice ' Global Production Change by GWL'], ['(Median Across GGCMs/ESMs)']});
   legendnames = {'SSP1-2.6','SSP3-7.0','SSP5-8.5'};
   if (strcmp(thiscrop,'mai')),
      legend(h(:),legendnames,'location','southwest');
   end;
   print(f,'-dpng',['/home/mmphill2/median_' thiscrop '_global_deltaproduction_allscen_bySSP_gwl_defaultco2_consistent.png']);
   print(f,'-depsc',['/home/mmphill2/median_' thiscrop '_global_deltaproduction_allscen_bySSP_gwl_defaultco2_consistent.eps']);
   close all

end;
