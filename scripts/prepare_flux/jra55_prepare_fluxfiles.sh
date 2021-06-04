#!/bin/sh
# cdo -b F32 -f nc4 copy 

domain=SO
source=jra55
datapath=/local/data/${source}_SO/
outpath=/local/projects/so_decadal_variability/flux/
inperiod=1979-2020
outperiod=1979-2020

rm ${datapath}${domain}_flux_ht_${source}_${outperiod}.nc
rm ${datapath}${domain}_flux_sr_${source}_${outperiod}.nc
cdo selname,var122 ${datapath}${domain}_jra55_shtfl_${inperiod}.nc ${datapath}tmp.nc
cdo setname,ht ${datapath}tmp.nc ${datapath}${domain}_flux_ht_${source}_${outperiod}.nc
cdo setname,sr ${datapath}tmp.nc ${datapath}${domain}_flux_sr_${source}_${outperiod}.nc
rm ${datapath}tmp.nc

rm ${datapath}${domain}_flux_fw_${source}_${outperiod}.nc
cdo selname,var61 ${datapath}${domain}_jra55_tprat_${inperiod}.nc ${datapath}tmp.nc
cdo setname,fw ${datapath}tmp.nc ${datapath}${domain}_flux_fw_${source}_${outperiod}.nc
rm ${datapath}tmp.nc
