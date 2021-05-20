#!/bin/sh
# cdo -b F32 -f nc4 copy 

domain=SO
source=merra2
datapath=/local/data/${source}/
outpath=/local/projects/so_decadal_variability/flux/
inperiod=198001-201912 
outperiod=1980-2019

rm ${outpath}${domain}_flux_sst_${source}_${outperiod}.nc
cdo selname,TSKINWTR ${datapath}merra2_ocn_mon_SO_198001-201912.nc4 ${outpath}tmp.nc
cdo setname,sst ${outpath}tmp.nc ${outpath}${domain}_flux_sst_${source}_${outperiod}.nc
rm ${outpath}tmp.nc
