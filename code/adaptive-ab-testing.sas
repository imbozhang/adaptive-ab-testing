
proc power; 
 twosamplefreq TEST=pchq 
 proportiondiff = 0.006
 refproportion = 0.03
 npergroup = .
 power = 0.85
 alpha=0.10;
 title "Sample Size Calculation for Traditional A/B Testing";
run;

/*Case 1: two-sided alpha=0.1, power=0.85, interim at 0.5;*/
proc iml;
  cp=t((10:90)/100);
  alpha=0.1;
  beta=0.85;
  zalpha=quantile("Normal",1-alpha/2,0,1);
  zbeta=quantile("Normal",beta,0,1);
  n=12744;
  Np=12744*2;
  t=n/Np;
  *interim observed zt;
  zt=(sqrt(t))*(quantile("Normal",cp,0,1)*sqrt(1-t)+zalpha);
  *sample size re-calculation;
  n2_0=n+(n/zt##2)#((zalpha*sqrt(Np)-zt#sqrt(n))/sqrt(Np-n)+zbeta)##2;
  n2=n2_0;
  n2=ceil(n2_0);
  n2[loc(n2_0<np)]=np;
  n2star=n2-n;
  *lower bound b(zt,n2);
  b=(sqrt(n2star/(np-n))#(zalpha#sqrt(np)-zt#sqrt(n))+zt#sqrt(n))/sqrt(n2);
  create out var {cp zt n2 b zalpha beta};
  append;
  close out;
quit;

proc sql noprint;
 select zalpha into: zalpha from out
 where zalpha >.
 ;
 select beta into: beta from out
 where beta >.
 ;
 select min(cp) into: CPlower from out
 where b<=&zalpha
 ;
quit;

/*CP = 0.31 is CPlower */
ods graphics on;
symbol1 interpol=join value=dot;
proc gplot data=out;
 plot b*cp/vref=&zalpha href=&CPlower 0.85;
 title "b versus CP";
 title2 "Promissing Zone Requires CP Such That b<=z_alpha"; 
 title3 "CPlower = %cmpres(&CPlower), Target Power = %cmpres(&beta)";
 run;

/*Test Case on Two-Stage A/B Testing*/
data test;
 input group $ response count;
 cards;
 A  1 191 
 A  0 6181
 B  1 216
 B  0 6156
 ;
run;
ods output ChiSq=ChiSq;
proc freq data=test;
 tables group*response/chisq;
 weight count;
 title "Observed Test Statistics at End of First Stage";
run;

proc sql noprint;
 select sqrt(Value) into :z
 from ChiSq
 where Statistic="Chi-Square"
 ;
 select sum(count) into :n
 from test
 ;
quit;

proc iml;
 alpha=0.1;
 beta=0.85;
 zalpha=quantile("Normal",1-alpha/2,0,1);
 zbeta=quantile("Normal",beta,0,1);
 n=&n;z=&z;
 Np=12744*2;
 t=n/Np;
 title "Conditional Power Based on Observed Test Statistics";
 *calculate conditional power;
 zz=(&z/sqrt(t)-zalpha)/sqrt(1-t);
 cp=round(cdf("Normal",zz,0,1),.001);
 print cp;
 *recalculate sample size;
 NZn=n+(n/(z**2))*(((zalpha*sqrt(Np)-z*sqrt(N))/sqrt(Np-n)+zbeta)**2);
 *Round up re-calculated sample size and make it an even number;
 N2=ceil(NZn/2)*2;
 create ReCal var {alpha beta zalpha zbeta z cp n2};
 append;
 close ReCal;
quit;

data decision;
 set ReCal;
 length decision $100.;
 if cp<&CPlower then decision="Stop A/B Testing";
 else if &CPlower<=cp<&beta then decision="Increase Total Sample Size to "||n2;
 else decision="Continue A/B Testing to Originally Planned Sample Size";
 label alpha="Significance Level"
       beta="Planned Power"
	   z="Test Statistics at End of First Stage"
	   cp="Conditional Power"
	   n2="Re-calculated Sample Size"
	   decision="Decision"
 ;
run;
proc print data=decision label noobs;
 var alpha beta z cp decision;
 title "Decision Made at the End of First Stage";
run;
