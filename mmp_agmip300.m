%		                  mmp_agmip300
%
%  This script converts GGCMI Phase 3 pre-processed
%  outputs (adjusted harvest dates, combined crop seasons,
%  combined rainfed and irrigated production, etc.) to
%  .mat files for plotting and other analyses. Based on
%  acr_agmip503.m. This is an updated version of
%  mmp_agmip200.m.
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
%				                         date:	8/11/2022
%
function mmp_agmip300();
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
gcmtypes = {'gfdl-esm4','ipsl-cm6a-lr','mpi-esm1-2-hr','mri-esm2-0','ukesm1-0-ll'};
gcmtypesnice = {'GFDL-ESM4','IPSL-CM6A-LR','MPI-ESM1-2-HR','MRI-ESM2-0','UKESM1-0-LL'};
croptypes = {'mai','whe','soy','ric'};
allyears = [1979:2099];

cd /discover/nobackup/mmphill2/GGCMI3_Processed/

%% Loop through models
for thismodel = 1:length(modeltypes),
   for thisscen = 1:length(scentypes),
      for thisgcm = 1:length(gcmtypes),

         filename = ['GGCMI_Phase3_annual_' scentypes{thisscen} '_' gcmtypes{thisgcm} '_' modeltypes{thismodel} '.nc4'];
         if(exist(filename)),
            deltayield = ncread([filename],'yield change');

            %% Loop through crops
            for thiscrop=1:length(croptypes),
               %% Prepare variable
               outvar = permute(squeeze(deltayield(:,:,79:end,thiscrop)),[2 1 3]);
               if (sum(sum(sum(~isnan(outvar))))>0),
               save(['/discover/nobackup/projects/giss_ana/users/mmphill2/AgGRID/matfiles/deltayield_' modeltypes{thismodel} '_' croptypes{thiscrop} '_total_' gcmtypes{thisgcm} '_' scentypes{thisscen} '_processed.mat'],'outvar');
               end;
               clear outvar

            end;           %% crop
         end;              %% if exists
     end;                  %% gcm
   end;                    %% scen
end;                       %% ggcm
