%%
clear all
clc
close all

datafile1='/remote/oahu/data/jra55_SO/SO_jra55_';
suffix = '_1979-2020.nc'
icefile= '/remote/oahu/data/jra55/seaice_SO_monthly_clim_perm2_1992-2008.nc';
outfile1='/remote/oahu/data/jra55_SO/SO_flux_ht_jra55_1979-2019.nc';
outfile2='/remote/oahu/data/jra55_SO/SO_flux_sr_jra55_1979-2019.nc';
outfile3='/remote/oahu/data/jra55_SO/SO_flux_fw_jra55_1979-2019.nc';
landfile='/remote/oahu/data/jra55/TL319_SO.nc';
%%
%cdo remapcon2,griddes_SO ../seaice_freshwater/prod_monthly_clim_perm2_1992-2008.nc seaice_SO_monthly_clim_perm2_1992-2008.nc
%%
fact=1/(60*60*24).*1000;%to convert from m per day to kg per s
ncid=netcdf.open(icefile,'NOWRITE');
ice=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'prod')));
netcdf.close(ncid)
clear ncid
ice(ice<=-9.9900e+20)=0;
ice(isnan(ice))=0;
ice=repmat(ice,[1 1 41]);
rhoice=925; %kg/m3 eg http://www.sciencedirect.com/science/article/pii/0165232X9500007X
rhowater=1000; %kg/m3 density of freshwater
s_ice=6/1000; %frac of salt in ice...see eg vancoppenolle
s_seawater=34.7/1000; %kg/m3 salinity of CDW
Li= 334*10^3;
conversion=(rhoice*(1-(s_ice/s_seawater)))/rhowater;
conversion=1./100./conversion./(24*60*60).*rhoice.*Li; %convert to heat flux associated with melting cm/day -> W/m2
ice_ht=ice.*100.*conversion;
ice_fw=-ice.*fact;
clear fact ice

% convention: positive=heat into ocean, W*s* m**-2, accumulated per day
fact=1;%1/(60*60*24);%to convert from W*s per day to W
ncid=netcdf.open([datafile1,'shtfl',suffix],'NOWRITE');
sshf=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'var122')));
netcdf.close(ncid)
ncid=netcdf.open([datafile1,'lhtfl',suffix],'NOWRITE');
slhf=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'var121')));
netcdf.close(ncid)
ncid=netcdf.open([datafile1,'dlwrf',suffix],'NOWRITE');
strd=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'var205')));
netcdf.close(ncid)
ncid=netcdf.open([datafile1,'dswrf',suffix],'NOWRITE');
ssrd=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'var204')));
netcdf.close(ncid)
ncid=netcdf.open([datafile1,'ulwrf',suffix],'NOWRITE');
stru=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'var212')));
netcdf.close(ncid)
ncid=netcdf.open([datafile1,'uswrf',suffix],'NOWRITE');
ssru=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'var211')));
netcdf.close(ncid)
clear ncid
sshf(sshf>1.e+19)=NaN;
slhf(slhf>1.e+19)=NaN;
strd(strd>1.e+19)=NaN;
ssrd(ssrd>1.e+19)=NaN;
stru(stru>1.e+19)=NaN;
ssru(ssru>1.e+19)=NaN;
ht=(-sshf-slhf+strd-stru).*fact;
sr=(ssrd-ssru).*fact;
clear ssrd ssru sshf slhf strd stru
ht(isnan(ht))=0;
sr(isnan(sr))=0;


ncid = netcdf.open(landfile,'NC_NOWRITE');
mask = squeeze(double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'LAND_GDS4_SFC'))));
netcdf.close(ncid);
mask(mask>1.e+19)=NaN;
mask(mask<0)=NaN;
mask=repmat(mask,[1 1 size(ht,3)]);
ht(mask==1)=0;
sr(mask==1)=0;
clear fact
ht=ht+ice_ht;
clear ice_ht


% convention: positive=freshwater into ocean, mm m^-2 accumulated per day
fact=1/(60*60*24).*1000./1000;%to convert from mm per day to kg per s

ncid=netcdf.open([datafile1,'evp',suffix],'NOWRITE');
e=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'var57')));
netcdf.close(ncid)
ncid=netcdf.open([datafile1,'tprat',suffix],'NOWRITE');
tp=double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,'var61')));
netcdf.close(ncid)
clear ncid
e(e>1.e+19)=NaN;
tp(tp>1.e+19)=NaN;
fw=(-e+tp).*fact;
clear e tp
fw(mask==1)=0;
clear fact
fw(isnan(fw))=0;

fw=fw+ice_fw;
clear ice_fw

%%
ncid=netcdf.open(outfile1,'WRITE');
varid=netcdf.inqVarID(ncid,'ht');
netcdf.putVar(ncid,varid,ht);
netcdf.reDef(ncid)
try
netcdf.delAtt(ncid,varid,'standard_name')
catch
end
try
netcdf.delAtt(ncid,varid,'long_name')
catch
end
netcdf.delAtt(ncid,varid,'units')
netcdf.putAtt(ncid,varid,'units','W m^-2');
netcdf.endDef(ncid)
netcdf.close(ncid)
clear ncid varid

ncid=netcdf.open(outfile2,'WRITE');
varid=netcdf.inqVarID(ncid,'sr');
netcdf.putVar(ncid,varid,sr);
netcdf.reDef(ncid)
try
netcdf.delAtt(ncid,varid,'standard_name')
catch
end
try
netcdf.delAtt(ncid,varid,'long_name')
catch
end
netcdf.delAtt(ncid,varid,'units')
netcdf.putAtt(ncid,varid,'units','W m^-2');
netcdf.endDef(ncid)
netcdf.close(ncid)
clear ncid varid

ncid=netcdf.open(outfile3,'WRITE');
varid=netcdf.inqVarID(ncid,'fw');
netcdf.putVar(ncid,varid,fw);
netcdf.reDef(ncid)
try
netcdf.delAtt(ncid,varid,'standard_name')
catch
end
try
netcdf.delAtt(ncid,varid,'long_name')
catch
end
netcdf.delAtt(ncid,varid,'units')
netcdf.putAtt(ncid,varid,'units','kg m^-2 s^-1');
netcdf.endDef(ncid)
netcdf.close(ncid)
clear ncid varid
