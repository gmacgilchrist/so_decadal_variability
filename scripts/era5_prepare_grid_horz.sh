#!/bin/sh

inpath=/local/projects/so_decadal_variability/flux/
outpath=/local/projects/so_decadal_variability/grid/
source=era5
period=1979-2019

rm ${outpath}SO_grid_des_${source}
cdo griddes ${inpath}SO_flux_sst_${source}_${period}.nc > ${outpath}SO_grid_des_${source}
rm ${outpath}SO_grid_area_${source}.nc
cdo gridarea ${inpath}SO_flux_sst_${source}_${period}.nc ${outpath}SO_grid_garea_${source}.nc
rm ${outpath}SO_grid_dx_${source}.nc
cdo griddx ${inpath}SO_flux_sst_${source}_${period}.nc ${outpath}SO_grid_dx_${source}.nc
rm ${outpath}SO_grid_dy_${source}.nc
cdo griddy ${inpath}SO_flux_sst_${source}_${period}.nc ${outpath}SO_grid_dy_${source}.nc
rm ${outpath}SO_grid_mask_${source}.nc
cdo setname,mask ${outpath}SO_grid_garea_${source}.nc ${outpath}SO_grid_mask_${source}.nc
