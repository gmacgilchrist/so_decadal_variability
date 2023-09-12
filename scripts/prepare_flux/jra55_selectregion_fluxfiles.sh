#!/bin/sh
# cdo -b F32 -f nc4 copy 

domain=SO
source=jra55
datapath=/local/data/${source}/
outpath=/local/data/jra55_SO/
inperiod=1979-2020 
outperiod=1979-2020
declare -a vars=("brtmp" "dlwrf" "dswrf" "evp" "icec" "lhtfl" "shtfl" "tprat" "ulwrf" "uswrf")
for v in "${vars[@]}"; do
	rm ${outpath}${domain}_${source}_${v}_${outperiod}.nc
	cdo sellonlatbox,-180,180,-90,-30 ${datapath}${source}_${v}_${inperiod}.nc ${outpath}${domain}_${source}_${v}_${outperiod}.nc
done
