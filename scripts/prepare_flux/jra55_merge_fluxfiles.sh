#!/bin/sh
# cdo -b F32 -f nc4 copy 

domain=SO
source=jra55
datapath=/local/data/${source}/
outpath=/local/projects/so_decadal_variability/flux/
inperiod=1979-2020 
outperiod=1979-2020

rm ${outpath}${domain}_flux_sst_${source}_${outperiod}.nc
cdo selname,var118 ${datapath}${source}_brtmp_${inperiod}.nc ${outpath}tmp.nc
cdo sellonlatbox,-180,180,-90,-30 ${outpath}tmp.nc ${outpath}${domain}_tmp.nc
cdo setname,sst ${outpath}${domain}_tmp.nc ${outpath}${domain}_flux_sst_${source}_${outperiod}.nc
rm ${outpath}*tmp.nc
