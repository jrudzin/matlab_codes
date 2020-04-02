function [mld,jmld]=get_mld(z,T)
mld=0;
jref=0;

%zref=24;  % usual case
zref=12;  % use this value when shallow OML are present

for j=1:length(z);
 if z(j) == zref;
  jref=j;
 end;
end;
if jref > 0;
 Tave=nanmean(T(1:jref));
 for j=jref:length(z);
  if (Tave-T(j)) > 0.5 & mld == 0;
   mld=z(j);
   mlT=nanmean(T(1:j));
   jmld=j;
  end;
 end;
end;
if mld > 0;
 mld=mld-2;
 jmld=jmld-1;
end;
