#!/bin/sh

domain=SO
source=en4
gridtype=jra55
# datapath=/local/data/${source}/v4.2.1/
# correction=g10
datapath=/local/projecta/so_decadal_variability/ocean/
outpath=/local/projects/so_decadal_variability/ocean/
inperiod=1979-2019 
outperiod=197901-201812
gridfile=/local/projects/so_decadal_variability/grid/SO_grid_des_${gridtype}
method=con2

CDO_RESET_HISTORY=1 ; export CDO_RESET_HISTORY

#echo cdo -P 4 remap${method},${gridfile} ${datapath}${correction}/EN.4.2.1.f.analysis.${correction}.197901-201812.nc ${outpath}${domain}_ocean_sa_${source}_${gridtype}_${outperiod}.nc
rm ${outpath}${domain}_ocean_sa_${source}_${gridtype}_${outperiod}.nc
cdo -P 4 remap${method},${gridfile} ${outpath}${domain}_ocean_sa_${source}_${outperiod}.nc ${outpath}${domain}_ocean_sa_${source}_${gridtype}_${outperiod}.nc
rm ${outpath}${domain}_ocean_alpha_${source}_${gridtype}_${outperiod}.nc
cdo -P 4 remap${method},${gridfile} ${outpath}${domain}_ocean_alpha_${source}_${outperiod}.nc ${outpath}${domain}_ocean_alpha_${source}_${gridtype}_${outperiod}.nc
rm ${outpath}${domain}_ocean_beta_${source}_${gridtype}_${outperiod}.nc
cdo -P 4 remap${method},${gridfile} ${outpath}${domain}_ocean_beta_${source}_${outperiod}.nc ${outpath}${domain}_ocean_beta_${source}_${gridtype}_${outperiod}.nc
