#!/bin/sh
# cdo -b F32 -f nc4 copy 

domain=SO
source=erai
source_long=era-interim
datapath=/local/data/${source_long}/
outpath=/local/projects/so_decadal_variability/flux/
inperiod=197901-201812 
outperiod=1979-2018

rm ${outpath}${domain}_flux_sst_${source}_${outperiod}.nc
cdo selname,sst ${datapath}erai_sst_SO_monthly_197901-201812.nc ${outpath}tmp.nc
cdo setname,sst ${outpath}tmp.nc ${outpath}${domain}_flux_sst_${source}_${outperiod}.nc
rm ${outpath}tmp.nc
