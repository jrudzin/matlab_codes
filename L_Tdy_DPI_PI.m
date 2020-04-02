function [alpha,ustar,L,Tdy,DPI,PI]=L_Tdy_DPI(ofname,fin_track,stormday,c_d,T0,Ta,qa)

%Created by Johna Rudzin (rudzinjohna@gmail.com), 1/2/2019

%Code created using equations from Balaguru et al., 2015, Emaunel 1999

%estimating mixing length scale, mixing length temperature,and DPI from T/S profiles

%%inputs

%ofname :: string, file name for T/S profile. Needs to be pre-storm T/S profiles. File format should be columns of : Z, T, S
%fin_track :: string, file name for best track file
%stormday :: integer, day of year of storm passage for area of T/S profile
%c_d :: integer, value of momentum drag coefficient to be used
%T0 :: integer, outflow temperature, in Kelvin
%qa :: integer, value of specific humidity of air at 10m

%%outputs

%alpha :: stratification term (kgm^-3/m)
%ustar :: surface friction velocity (m/s)
%L :: mixing length (m)
%Tdy :: depth-integrated (over L) temperature (degC)
%DPI :: dynamic potential intensity (m/s)
%PI :: potential intensity (m/s)

%constants
kappa=0.4; %von Karman unitless
rho0=1025; %density of seawater kg/m^-3
rhoair=1.2; %density of air kg/m^3
g=9.81; %gravity m/s^2
Lv=2.5*10^6; %J/kg
cp=1004.7; %J/K/kg

addpath /programs/toolboxes/sw

odata=load(ofname); %load T/S profiles from .mat file

Z=odata(:,1);T=odata(:,2);S=odata(:,3);
sst(ii)=T(1);

%%%get constant density mld as in de Boyer Montegut et al. (2007)

st=sw_pden(S,T,Z,Z(1)); %estimate sigma-t

deltasig=sw_pden(S(1),T(1)-0.5,Z(1),Z(1))-sw_pden(S(1),T(1),Z(1),Z(1));
ind=find(round(st*10)./10-0 == round((st(1)+deltasig)*10)./10); 
ind=ind(1); %incase theres more than one index, get top index

mld=Z(ind);

%constant temperature mixed layer
[ild,jild]=get_mld(Z,T); 

%get rate of change with depth of potential density from ML to IL

mld50=mld+50; %add 50m depth to MLD to get bottom depth limit for estimate as suggested in correspondence from K. Balaguru

ind50=find(ppint == mld50(ii));

%stratification term
alpha(ii)=(st(ind50)-st(ind))/(mld50(ii)-mld(ii));

%estimate ustar

%read in best track data file
[nday0 lon0 lat0 wspd p0 rmax0 cat0 ]=textread(fin_track,'%f %f %f %f %f %f %f');

nind=find(rmax0 == -99);rmax0(nind)=NaN; %
rmax=rmax0*1.852*1000; %nm to m
wspd=wspd*0.514; %convert to m/s from knots

bind=find(nday0 == stormday); get index of storm time passage

%%% Uh - calculates the distance between the two track points and then
%%% calculates translation speed
for i=1:length(lat0)-1;
    dist(i)=deg2km(distance(lon0(i),lat0(i),lon0(i+1),lat0(i+1)));
    Uh(i+1)=dist(i)/6; %km/hr
    Uh(i+1)=Uh(i+1)*(1000/3600); %m/s
end

t=rmax(bind)/Uh(bind); %time period of mixing

tau=rhoair*c_d*wspd(bind)^2; %momentum flux
ustar=sqrt(tau/rho0); %sfc friction velocity

%%%%Estimate L and Tdy

L(i) = mld + ((2*rho0*ustar.^3*t)/(kappa*g*alpha2))^(1/3);

cumsum=0;
for i=1:round(L)

    cumsum=cumsum+ T(i);
end
Tdy=(1/round(L))*cumsum;

%%%%%%Estimate PI & DPI

c_k=c_d;
sst_kelvin=sst+273.15;
Tdy_kelvin=Tdy+273.15;
qsst=qsat(sst);
qTdy=qsat(Tdy);

%enthalpy terms
ksst=(Lv*qsst+cp*sst_kelvin);
kTdy=(Lv*qTdy+cp*Tdy_kelvin);
k=(Lv*qair+cp*Ta);

PI=sqrt(((sst_kelvin-T0)/T0)*(c_k/c_d)*(ksst-k)); %m/s
DPI=sqrt(((Tdy_kelvin-T0)/T0)*(c_k/c_d)*(kTdy-k)); %m/s

