{Unit to implement some probability functions, v. 1.0.
Copyright (C) 2005 Franco Milani
email: nadir@solaris.polarhome.com or astlab@inwind.it
web sites:
http://www.polarhome.com/~franco
http://www.polarhome.com:723/~nadir/software
http://spazioinwind.libero.it/frm/software
-----------------------------------------------------------------------------
This software is released in the hope that it will be useful but on an as-is
basis, so without any warranty, either expressed or implied. It can be used
and modified for personal purposes, or linked to distributed programs,
in which case an acknowledgment will be appreciated. This software can also
be redistributed, but only in the original form. Possible changes in the
distributed versions will be made only by the author, also according to
hints, requests or bug reports.
-----------------------------------------------------------------------------
This unit contains some functions useful in probability and statistics.
To achieve the highest possible precision with extended real variables,
the routines have been translated from the Cephes Mathematical Library
written in C by Stephen L. Moshier [1].
The accuracy in the usual ranges, not too near to overflow or underflow,
is of 14 or more significant decimal digits.
The following is a listing of the available functions and some of their
characteristics ("inf" means infinity).

gamma(x)       gamma function, defined by

                                     +inf
                                       - 
                                      | |   x-1
               if x > 0  gamma(x) =   |    t   exp(-t)dt
                                    | |
                                     -
                                     0
                                                pi
               if x < 0  gamma(x) = - ------------------------
                                      gamma(1+|x|)*sin(pi*|x|)

               Domain: ]0,+inf[ and x<0 except negative integers.
               Co-domain: ]-inf,0[ U ]0,+inf[.


lngamma(x)     natural logarithm of the absolute value of the gamma function:
               lngamma(x) = ln(|gamma(x)|). Domain: ]0,+inf[ and x<0 except
               negative integers. Co-domain: ]-inf,+inf[.

            
gami(a,x)      incomplete gamma function, defined by

                                          x
                                          -
                                1        | |  -t  a-1
               gami(a,x)  =  --------    |   e   t   dt
                             gamma(a)  | |
                                        -
                                        0

               with a>0 and x>0. Domain: [0,+inf[. Co-domain: [0,1[. 


cgami(a,x)     complemented incomplete gamma function, defined by
 
                                                         inf
                                                          -
                                                1        | |  -t  a-1
               cgami(a,x) = 1 - gami(a,x) =  --------    |   e   t   dt
                                             gamma(a)  | |
                                                        -
                                                        x
 
               with a>0 and x>0. Domain: [0,+inf[. Co-domain: ]0,1].


icgami(a,x)    inverse of complemented incomplete gamma function.
               If y = cgami(a,x) then x = icgami(a,y). Positive arguments.
               Domain: ]0,1]. Co-domain: [0,+inf[.

betai(a,b,x)   incomplete beta integral, defined by

                                                   x
                                                   -
                                 gamma(a+b)       | |  a-1     b-1
               betai(a,b,x) = -----------------   |   t   (1-t)   dt.
                              gamma(a)*gamma(b) | |
                                                 -
                                                 0
 
               with a>0, b>0 and 0<=x<=1.


ibetai(a,b,x)  inverse incomplete beta integral.
               If y = betai(a,b,x) then x = ibetai(a,b,y).
              

erf(a)         error function, defined by

                                       a 
                                       -
                             2        | |        2
               erf(a)  =  --------    |    exp(-t )dt
                          sqrt(pi)  | |
                                     -
                                     0
             
               Domain: ]-inf,+inf[. Co-domain: ]-1,1[.
               Underflow with result -1 if a < UERF (about -5.93022),
               overflow with result 1 if a > OERF (about 5.93027).


ierf(a)        inverse error function. If y = erf(a) then a = ierf(y).
               Domain: ]-1,1[. Co-domain: ]-inf,+inf[.


erfc(a)        complementary error function, defined by

                                                     +inf
                                                       -
                                             2        | |        2
               erfc(a)  =  1 - erf(a)  =  --------    |    exp(-t )dt
                                          sqrt(pi)  | |
                                                     -
                                                     a
  
               Domain: ]-inf,+inf[. Co-domain: ]0,2[. 
               Underflow with result 0 if a > URFC (about 106.7296),
               overflow with result 2 if a < UERF (about -5.93022).


ierfc(a)       inverse complementary error function. If y = erfc(a) then
               a = ierfc(y). Domain: ]0,2[. Co-domain: ]-inf,+inf[.


phi(a)         cumulative probability of the standard normal distribution:
               is the probability p(0<=x<=a) defined by

                                        a
                                        -
                              1        | |        2
               phi(a)  =  ---------    |    exp(-t /2)dt
                          sqrt(2pi)  | |
                                      -
                                      0

               = 0.5*erf(a/sqrt(2)).
               Domain: [0,+inf[. Co-domain: [0,1/2[. Overflow with result
               0.5 if a > OPHI (about 8.38667).


iphi(a)        inverse cumulative probability of the standard normal
               distribution. If y = phi(x) then x = iphi(y).
               Domain: [0,1/2[. Co-domain: [0,+inf[.
               

lnor(a)        lower tail cumulative probability of the standard normal
               distribution: is the probability p(x<=a) defined by

                                         a
                                         -
                               1        | |        2
               lnor(a)  =  ---------    |    exp(-t /2)dt
                           sqrt(2pi)  | |
                                       -
                                     -inf
  
               = (1 + erf(a/sqrt(2)))/2 = 1 - erfc(a/sqrt(2))/2.
               Domain: ]-inf,+inf[. Co-domain: ]0,1[.
               Underflow with result 0 if a < ULNR (about -150.9385),
               overflow with result 1 if a > OLNR (about 8.30476).


ilnor(a)       inverse lower tail cumulative probability of the standard
               normal distribution. If y = lnor(a) then a = ilnor(y).
               Domain: ]0,1[. Co-domain: ]-inf,+inf[.


unor(a)        upper tail cumulative probability of the standard normal
               distribution. Complementary of lnor, is the probability
               p(x>=a) defined by

                                                   +inf
                                                     -
                                           1        | |        2
               unor(a) = 1 - lnor(a) = ---------    |    exp(-t /2)dt
                                       sqrt(2pi)  | |
                                                   -
                                                   a

               = (1 - erf(a/sqrt(2)))/2 = erfc(a/sqrt(2))/2.
               Domain: ]-inf,+inf[. Co-domain: ]0,1[. Underflow with
               result 0 if a > UUNR (about 9.15529), overflow with
               result 1 if a < OUNR (about -8.30476).


iunor(a)       inverse upper tail cumulative probability of the standard
               normal distribution. If y = unor(a) then a = iunor(y).
               Domain: ]0,1[. Co-domain: ]-inf,+inf[.


chsq(a,v)      lower cumulative probability of the chi-square distribution,
               with the number of degrees of freedom v as parameter.
               Domain: [0,+inf[. Co-domain: ]0,1]. Returns the probability
               p(x<=a) defined by

                                                a
                                                -
                                   1           | |  v/2-1 -t/2
               chsq(a,v)  =  --------------    |   t     e    dt
                              v/2            | |
                             2   gamma(v/2)   -
                                              0

               Overflow with result 1 if a > about 70 and v=1. This limit
               increases if v increases.


ichsq(a,v)     inverse lower cumulative probability of the chi-square
               distribution with v degrees of freedom. If y = chsq(a,v)
               then a = ichsq(y,v). Domain: ]0,1]. Co-domain: [0,+inf[.


chsqc(a,v)     upper tail cumulative probability of the chi-square
               distribution, with the number of degrees of freedom v
               as parameter. Complementary of chsq. Domain: [0,+inf[.
               Co-domain: ]0,1]. Returns the probability P(x>=a) defined by

                                               +inf
                                                 -
                                    1           | |  v/2-1 -t/2
               chsqc(a,v)  =  --------------    |   t     e    dt
                               v/2            | |
                              2   gamma(v/2)   -
                                               a
  
               Underflow with result 0 if a > about 22788 and v=1.
               This limit increases if v increases.


ichsqc(a,v)    inverse upper tail cumulative probability of the chi-square
               distribution with v degrees of freedom. If y = chsqc(a,v)
               then a = ichsqc(y,v). Domain: ]0,1] with the limitation
               a > MICG (1e-4932) to avoid overflow. Co-domain: [0,+inf[.

-----------------------------------------------------------------------------

References

[1] Cephes Mathematical Library, available at the Moshier's site
    http://www.moshier.net
[2] M. Abramowitz and I. A. Stegun, "Handbook of Mathematical Functions",
    National Bureau of Standards, 1972
[3] H. B. Dwight, "Tables of Integrals and Other Mathematical Data",
    MacMillan Company, 1961

-----------------------------------------------------------------------------}
UNIT ustat10;

INTERFACE

{gamma function}
FUNCTION gamma(x: extended): extended;

{natural logarithm of the absolute value of the gamma function}
FUNCTION lngamma(x: extended): extended;

{incomplete gamma integral}
FUNCTION gami(a: extended; x: extended): extended;

{complemented incomplete gamma integral}
FUNCTION cgami(a: extended; x: extended): extended;

{inverse of complemented incomplete gamma integral}
FUNCTION icgami(a: extended; x: extended): extended;

{incomplete beta integral}
FUNCTION betai(a,b,x: extended): extended;

{inverse incomplete beta integral}
FUNCTION ibetai(a,b,x: extended): extended;

{error function}
FUNCTION erf(a: extended): extended;

{inverse error function}
FUNCTION ierf(a: extended): extended;

{complementary error function}
FUNCTION erfc(a: extended): extended;

{inverse complementary error function}
FUNCTION ierfc(a: extended): extended;

{cumulative probability of the standard normal distribution}
FUNCTION phi(a: extended): extended;

{inverse cumulative probability of the standard normal distribution}
FUNCTION iphi(a: extended): extended;

{lower tail cumulative probability of the standard normal distribution}
FUNCTION lnor(a: extended): extended;

{inverse lower tail cumulative probability of the standard normal distribution}
FUNCTION ilnor(a: extended): extended;

{upper tail cumulative probability of the standard normal distribution}
FUNCTION unor(a: extended): extended;

{inverse upper tail cumulative probability of the standard normal distribution}
FUNCTION iunor(a: extended): extended;

{lower cumulative probability of the chi-square distribution with v degrees
of freedom}
FUNCTION chsq(a: extended; v: longword): extended;

{inverse lower cumulative probability of the chi-square distribution with
v degrees of freedom}
FUNCTION ichsq(a: extended; v: longword): extended;

{upper tail cumulative probability of the chi-square distribution with v
degrees of freedom}
FUNCTION chsqc(a: extended; v: longword): extended;

{inverse upper tail cumulative probability of the chi-square distribution
with v degrees of freedom}
FUNCTION ichsqc(a: extended; v: longword): extended;

IMPLEMENTATION

CONST

CSPI: extended = 0.31415926535897932384626E+1; {pi}
SQT2: extended = 1.4142135623730950488;        {sqrt(2)}
S2PI: extended = 2.506628274631000502416;      {sqrt(2*pi)}
L2PI: extended = 0.91893853320467274178;       {ln(sqrt(2*pi))}
XPM2: extended = 1.35335283236612691894E-1;    {exp(-2)}
ISQ2: extended = 7.07106781186547524401E-1;    {1/sqrt(2)}
XPML: extended = 8.6466471676338730811E-1;     {1-exp(-2)}
MXLL: extended =-1.1355830259113584004E+4;
LLIM: extended = 0.03125;
MXLG: extended = 1.04848146839019521116E+4928;
MXGL: extended = 1.1356523406294143949492E+4;
MIGL: extended =-1.1355137111933024058873E+4;
MXNL: extended = 1.189731495357231765021E+4932;
MCHP: extended = 5.421010862427522170037E-20;
NBIG: extended = 9.223372036854775808E+18;
BINV: extended = 1.084202172485504434007E-19;
UERF: extended =-5.9302291691084904282;        {erf and erfc limit}
OERF: extended = 5.9302742314456255421;        {erf overflow point}
URFC: extended = 1.0672964589947618630E+2;     {erfc underflow point}
OPHI: extended = 8.386674246702086461;         {phi overflow point}
ULNR: extended =-1.50938512738317213808E+2;    {lnor underflow point}
OLNR: extended = 8.3047644693506286744;        {lnor overflow point}
UUNR: extended = 9.1552937726860725459;        {unor underflow point}
OUNR: extended =-8.3047644693506286744;        {unor overflow point}
OGAM: extended = 1755.455;                     {gamma overflow point}
MICG: extended = 1e-4932;                      {icgami overflow point}

{coefficients used to approximate by rational functions}

Aj: array [0..7] of extended =
(3.387132872796366608,1.3314166789178437745E+2,1.9715909503065514427E+3,
 1.3731693765509461125E+4,4.5921953931549871457E+4,6.7265770927008700853E+4,
 3.3430575583588128105E+4,2.5090809287301226727E+3);
Bj: array [1..7] of extended =
(4.2313330701600911252E+1,6.871870074920579083E+2,5.3941960214247511077E+3,
 2.1213794301586595867E+4,3.930789580009271061E+4,2.8729085735721942674E+4,
 5.226495278852854561E+3);
Cj: array [0..7] of extended =
(1.42343711074968357734,4.6303378461565452959,5.7694972214606914055,
 3.64784832476320460504,1.27045825245236838258,2.4178072517745061177E-1,
 2.27238449892691845833E-2,7.7454501427834140764E-4);
Dj: array [1..7] of extended =
(2.05319162663775882187,1.6763848301838038494,6.8976733498510000455E-1,
 1.4810397642748007459E-1,1.51986665636164571966E-2,5.475938084995344946E-4,
 1.05075007164441684324E-9);
Ej: array [0..7] of extended =
(6.6579046435011037772,5.4637849111641143699,1.7848265399172913358,
 2.9656057182850489123E-1,2.6532189526576123093E-2,1.2426609473880784386E-3,
 2.71155556874348757815E-5,2.01033439929228813265E-7);
Fj: array [1..7] of extended =
(5.9983220655588793769E-1,1.3692988092273580531E-1,1.48753612908506148525E-2,
 7.868691311456132591E-4,1.8463183175100546818E-5,1.4215117583164458887E-7,
 2.04426310338993978564E-15);
Ph: array [0..9] of extended =
(1.130609921802431462353E9,2.290171954844785638925E9,2.295563412811856278515E9,
 1.448651275892911637208E9,6.234814405521647580919E8,1.87009507112043671593E8,
 3.833161455208142870198E7,4.964439504376477951135E6,3.198859502299390825278E5,
-9.085943037416544232472E-6);
Qh: array [0..9] of extended =
(1.130609910594093747762E9,3.56592869656703138891E9,5.188672873106859049556E9,
 4.58801818891860972689E9,2.729005809811924550999E9,1.138778654945478547049E9,
 3.358653716579278063988E8,6.822450775590265689648E7,8.79923997735126107761E6,
 5.669830829076399819566E5);
Rh: array [0..4] of extended =
(3.621349282255624026891,7.173690522797138522298,3.445028155383625172464,
 5.537445669807799246891E-1,2.697535671015506686136E-2);
Sh: array [0..4] of extended =
(1.072884067182663823072E1,1.533713447609627196926E1,6.572990478128949439509,
 1.005392977603322982436,4.781257488046430019872E-2);
Tk: array [0..6] of extended =
(1.097496774521124996496E-1,5.402980370004774841217,2.871822526820825849235E2,
 2.677472796799053019985E3,4.825977363071025440855E4,1.549905740900882313773E5,
 1.104385395713178565288E6);
Uk: array [0..5] of extended =
(4.525777638142203713736E1,9.715333124857259246107E2,1.245905812306219011252E4,
 9.942956272177178491525E4,4.636021778692893773576E5,9.787360737578177599571E5);
U0: array [0..7] of extended =
(8.779679420055069160496E-3,-7.649544967784380691785E-1,2.971493676711545292135,
-4.144980036933753828858,2.765359913000830285937,-9.570456817794268907847E-1,
 1.659219375097958322098E-1,-1.140013969885358273307E-2);
V0: array [0..6] of extended =
(-5.303846964603721860329,9.908875375256718220854,-9.031318655459381388888,
  4.496118508523213950686,-1.250016921424819972516,1.823840725000038842075E-1,
 -1.088633151006419263153E-2);
U1: array [0..9] of extended =
(4.302849750435552180717,4.3602094518370966826E1,9.454613328844768318162E1,
 9.336735653151873871756E1,5.305046472191852391737E1,1.775851836288460008093E1,
 3.640308340137013109859,3.69135490017122412239E-1,1.403530274998072987187E-2,
 1.377145111380960566197E-4);
V1: array [0..8] of extended =
(2.001425109170530136741E1,7.079893963891488254284E1,8.033277265194672063478E1,
 5.034715121553662712917E1,1.779820137342627204153E1,3.845554944954699547539,
 3.993627390181238962857E-1,1.52687068952219119138E-2,1.4987006762866754669E-4);
U2: array [0..7] of extended =
(3.244525725312906932464,6.856256488128415760904,3.765479340423144482796,
 1.240893301734538935324,1.740282292791367834724E-1,9.08283420099310744175E-3,
 1.617870121822776093899E-4,7.377405643054504178605E-7);
V2: array [0..6] of extended =
(6.021509481727510630722,3.528463857156936773982,1.289185315656302878699,
 1.87429014261570360951E-1,9.867655920899636109122E-3,1.760452434084258930442E-4,
 8.028288500688538331773E-7);
U3: array [0..7] of extended =
(2.020331091302772535752,2.133020661587413053144,2.114822217898707063183E-1,
-6.500909615246067985872E-3,-7.279315200737344309241E-4,-1.275404675610280787619E-5,
-6.433966387613344714022E-8,-7.772828380948163386917E-11);
V3: array [0..6] of extended =
(2.278210997153449199574,2.345321838870438196534E-1,-6.916708899719964982855E-3,
-7.908542088737858288849E-4,-1.387652389480217178984E-5,-7.001476867559193780666E-8,
-8.458494263787680376729E-11);
Sr: array [0..8] of extended =
(-1.193945051381510095614E-3,7.220599478036909672331E-3,-9.622023360406271645744E-3,
 -4.219773360705915470089E-2,1.665386113720805206758E-1,-4.200263503403344054473E-2,
 -6.558780715202540684668E-1,5.772156649015328608253E-1,1.0);
Sn: array [0..8] of extended =
(1.13337416724389438201E-3,7.220837261893170325704E-3,9.621911155035976733706E-3,
-4.219773343731191721664E-2,-1.665386113944413519335E-1,-4.200263503402112910504E-2,
 6.558780715202536547116E-1,5.772156649015328608727E-1,-1.0);
Ar: array [0..6] of extended =
(4.885026142432270781165E-3,-1.880801938119376907179E-3,8.412723297322498080632E-4,
-5.952345851765688514613E-4,7.936507795855070755671E-4,-2.77777777775034960344E-3,
 8.333333333333331447505E-2);
Br: array [0..6] of extended =
(-2.16369082764381285764E3,-8.72387152284351145979E4,-1.104326814691464261197E6,
 -6.111225012005214299996E6,-1.625568062543700591014E7,-2.003937418103815175475E7,
 -8.875666783650703802159E6);
Cr: array [0..6] of extended =
(-5.139481484435370143617E2,-3.403570840534304670537E4,-6.227441164066219501697E5,
 -4.81494037941188218663E6,-1.785433287045078156959E7,-3.138646407656182662088E7,
 -2.099336717757895876142E7);
Tw: array [0..5] of extended =
(6.97281375836585777429E-5,7.84039221720066627474E-4,-2.29472093621399176955E-4,
-2.68132716049382716049E-3,3.47222222222222222222E-3,8.33333333333333333333E-2);
St: array [0..8] of extended =
(7.147391378143610789273E-4,-2.363848809501759061727E-5,-5.950237554056330156018E-4,
 6.98933226062319317187E-5,7.840334842744753003862E-4,-2.294719747873185405699E-4,
-2.681327161876304418288E-3,3.472222222230075327854E-3,8.333333333333331800504E-2);
Pz: array [0..7] of extended =
(4.212760487471622013093E-5,4.5429319606080091556E-4,4.092666828394035500949E-3,
 2.385363243461108252554E-2,1.113062816019361559013E-1,3.629515436640239168939E-1,
 8.378004301573126728826E-1,1.000000000000000000009);
Qz: array [0..8] of extended =
(-1.397148517476170440917E-5,2.346584059160635244282E-4,-1.237799246653152231188E-3,
 -7.955933682494738320586E-4,2.773706565840072979165E-2,-4.633887671244534213831E-2,
 -2.243510905670329164562E-1,4.150160950588455434583E-1,9.999999999999999999908E-1);
Sz: array [0..8] of extended =
(-1.193945051381510095614E-3,7.220599478036909672331E-3,-9.622023360406271645744E-3,
 -4.219773360705915470089E-2,1.665386113720805206758E-1,-4.200263503403344054473E-2,
 -6.558780715202540684668E-1,5.772156649015328608253E-1,1.0);
Tz: array [0..8] of extended =
(1.13337416724389438201E-3,7.220837261893170325704E-3,9.621911155035976733706E-3,
-4.219773343731191721664E-2,-1.665386113944413519335E-1,-4.200263503402112910504E-2,
 6.558780715202536547116E-1,5.772156649015328608727E-1,-1.0);

VAR
vx,pb,qb,xb,yb,zb,yu,zu,ya,xe,ye,ze,x2,y2,z2,ac,bc,xc,xv,t,w,y,r,m: extended;
ap,aq,au,az,fz,ans,ax,uc,ur,vc,yc,vr,vt,vy,vz,pv,qv,zv,uv,sa,ta,ua: extended;
pk,p1,p2,qk,q1,q2,qx,yl,yh,x0,x1,zy,zd,lg,dt,df,sy,sw,sv,va,na,t1: extended;
xk,pa,qa,m1,m2,m3,m4,rz,zk,tv,ad,ts,bp,bq,k1,k2,k3,k4,k5,k6,k7,k8: extended;
av,bv,y0,dz,yz,xz,xq,xs,lm,yq,di,dr,yr,yk,xt,zr,ai: extended;
ir,rf,nf,iv,ci,nz: longword;
dq,dir: longint;
sgngam: smallint;
fg,fl,nc,flag: byte;

{natural logarithm of the absolute value of the gamma function}
FUNCTION lngamma(x: extended): extended;
BEGIN
if x<13.0 then begin
az:=1.0;
ap:=int(x+0.5);
fz:=x-ap;
au:=x;
while au>=3.0 do begin
ap:=ap-1.0;
au:=ap+fz;
az:=az*au
end;
fl:=0;
while au<2.0 do begin
if abs(au)<=LLIM then begin
fl:=1;
break
end;
az:=az/(ap+fz);
ap:=ap+1.0;
au:=ap+fz
end;
if fl=1 then begin
if au=0.0 then lngamma:=MXNL
else begin
if au<0.0 then begin
au:=-au;
ap:=au*((((((((Sn[0]*au+Sn[1])*au+Sn[2])*au+Sn[3])*au+Sn[4])*au+Sn[5])*au
+Sn[6])*au+Sn[7])*au+Sn[8]);
ap:=az/ap
end
else begin
ap:=au*((((((((Sr[0]*au+Sr[1])*au+Sr[2])*au+Sr[3])*au+Sr[4])*au+Sr[5])*au
+Sr[6])*au+Sr[7])*au+Sr[8]);
ap:=az/ap
end;
if ap<0.0 then ap:=-ap;
lngamma:=ln(ap)
end
end
else begin
if az<0.0 then az:=-az;
if au=2.0 then lngamma:=ln(az)
else begin
au:=ap-2.0+fz;
ap:=au*((((((Br[0]*au+Br[1])*au+Br[2])*au+Br[3])*au+Br[4])*au+Br[5])*au+Br[6])/
(((((((au+Cr[0])*au+Cr[1])*au+Cr[2])*au+Cr[3])*au+Cr[4])*au+Cr[5])*au+Cr[6]);
ap:=ap+ln(az);
lngamma:=ap
end
end
end
else begin
if x>MXLG then lngamma:=MXNL
else begin
aq:=(x-0.5)*ln(x)-x+L2PI;
if x<=1.0e10 then begin
ap:=1.0/(x*x);
aq:=aq+((((((Ar[0]*ap+Ar[1])*ap+Ar[2])*ap+Ar[3])*ap+Ar[4])*ap+Ar[5])*ap+Ar[6])/x
end;
lngamma:=aq
end
end
END;

{complemented incomplete gamma integral}
FUNCTION cgami(a: extended; x: extended): extended;
BEGIN
if (a<=0.0) or (x<=0.0) then cgami:=1.0
else begin
ax:=a*ln(x)-x-lngamma(a);
ax:=exp(ax);
if (x<1.0) or (x<a) then begin
ur:=a;
uc:=1.0;
ans:=1.0;
repeat
ur:=ur+1.0;
uc:=uc*x/ur;
ans:=ans+uc;
until uc/ans<=MCHP;
cgami:=1.0-ans*ax/a
end
else begin
vy:=1.0-a;
vz:=x+vy+1.0;
vc:=0.0;
p2:=1.0;
q2:=x;
p1:=x+1.0;
q1:=vz*x;
ans:=p1/q1;
repeat
vc:=vc+1.0;
vy:=vy+1.0;
vz:=vz+2.0;
yc:=vy*vc;
pk:=p1*vz-p2*yc;
qk:=q1*vz-q2*yc;
if qk<>0 then begin
vr:=pk/qk;
vt:=abs((ans-vr)/vr);
ans:=vr
end
else vt:=1.0;
p2:=p1;
p1:=pk;
q2:=q1;
q1:=qk;
if abs(pk)>NBIG then begin
p2:=p2/NBIG;
p1:=p1/NBIG;
q2:=q2/NBIG;
q1:=q1/NBIG
end;
until vt<=MCHP;
cgami:=ans*ax
end
end
END;

{incomplete gamma integral}
FUNCTION gami(a: extended; x: extended): extended;
BEGIN
if (x<=0) or (a<=0) then gami:=0.0
else if (x>1.0) and (x>a) then gami:=1.0-cgami(a,x)
else begin
ax:=a*ln(x)-x-lngamma(a);
ax:=exp(ax);
ur:=a;
uc:=1.0;
ans:=1.0;
repeat
ur:=ur+1.0;
uc:=uc*x/ur;
ans:=ans+uc;
until uc/ans<=MCHP;
gami:=ans*ax/a
end
END;

{inverse of complemented incomplete gamma integral}
FUNCTION icgami(a: extended; x: extended): extended;
BEGIN
x0:=MXNL;
yl:=0.0;
x1:=0.0;
yh:=1.0;
dt:=4.0*MCHP;
zd:=1.0/(9.0*a);
zy:=1.0-zd-ilnor(x)*sqrt(zd);
qx:=a*zy*zy*zy;
lg:=lngamma(a);
nc:=1;
for ci:=1 to 10 do begin
if (qx>x0) or (qx<x1) then break;
zy:=cgami(a,qx);
if (zy<yl) or (zy>yh) then break;
if zy<x then begin
x0:=qx;
yl:=zy
end
else begin
x1:=qx;
yh:=zy
end;
zd:=(a-1.0)*ln(x0)-x0-lg;
if zd<MXLL then break;
zd:=(x-zy)/exp(zd);
qx:=qx-zd;
if ci>2 then
if abs(zd/qx)<dt then begin
nc:=0;
break
end
end;
if nc=1 then begin
if x0=MXNL then begin
zd:=0.0625;
if qx<=0.0 then qx:=1.0;
repeat
qx:=(zd+1.0)*qx;
zy:=cgami(a,qx);
if zy<x then break;
zd:=zd+zd;
until false;
x0:=qx;
yl:=zy
end;
zd:=0.5;
dir:=0;
for ci:=1 to 400 do begin
qx:=x1+zd*(x0-x1);
zy:=cgami(a,qx);
lg:=(x0-x1)/(x1+x0);
if abs(lg)<dt then break;
lg:=(zy-x)/x;
if abs(lg)<dt then break;
if qx<=0.0 then break;
if zy>x then begin
x1:=qx;
yh:=zy;
if dir<0 then begin
dir:=0;
zd:=0.5
end
else if dir>1 then zd:=0.5*zd+0.5
else zd:=(x-yl)/(yh-yl);
inc(dir)
end
else begin
x0:=qx;
yl:=zy;
if dir>0 then begin
dir:=0;
zd:=0.5
end
else if dir<-1 then zd:=0.5*zd
else zd:=(x-yl)/(yh-yl);
dec(dir)
end
end
end;
icgami:=qx;
END;

{Stirling's approximation}
FUNCTION stirling(x: extended): extended;
BEGIN
sw:=1.0/x;
if x>1024.0 then begin
sw:=(((((Tw[0]*sw+Tw[1])*sw+Tw[2])*sw+Tw[3])*sw+Tw[4])*sw+Tw[5])*sw+1.0;
sv:=exp((0.5*x-0.25)*ln(x));
sy:=sv*(sv/exp(x))
end
else begin
sw:=1.0+sw*((((((((St[0]*sw+St[1])*sw+St[2])*sw+St[3])*sw+St[4])*sw
+St[5])*sw+St[6])*sw+St[7])*sw+St[8]);
sy:=exp((x-0.5)*ln(x))/exp(x)
end;
stirling:=S2PI*sy*sw;
END;

{gamma function}
FUNCTION gamma(x: extended): extended;
BEGIN
sgngam:=1;
if x<0.0 then qv:=-x else qv:=x;
if qv>13.0 then begin
if qv>OGAM then begin
writeln('gamma overflow, the argument must be < ',OGAM:1:3);
gamma:=0.0;
exit
end;
if x<0.0 then begin
pv:=int(qv);
if pv=qv then begin
writeln('gamma overflow on integer negative argument');
gamma:=0.0;
exit
end;
iv:=trunc(pv);
if iv and 1 = 0 then sgngam:=-1;
zv:=qv-pv;
if zv>0.5 then begin
pv:=pv+1.0;
zv:=qv-pv
end;
zv:=qv*sin(CSPI*zv);
zv:=abs(zv)*stirling(qv);
if zv<=CSPI/MXNL then begin
writeln('gamma overflow');
gamma:=0.0;
exit
end;
zv:=CSPI/zv
end
else zv:=stirling(x);
qv:=sgngam*zv
end
else begin
zv:=1.0;
uv:=x;
while uv>=3.0 do begin
uv:=uv-1.0;
zv:=zv*uv
end;
while uv<-0.03125 do begin
zv:=zv/uv;
uv:=uv+1.0
end;
if uv<=0.03125 then begin
if uv=0.0 then begin
writeln('gamma overflow');
gamma:=0.0;
exit
end;
if uv<0.0 then begin
uv:=-uv;
pv:=uv*((((((((Tz[0]*uv+Tz[1])*uv+Tz[2])*uv+Tz[3])*uv+Tz[4])*uv+Tz[5])*uv
+Tz[6])*uv+Tz[7])*uv+Tz[8])
end
else
pv:=uv*((((((((Sz[0]*uv+Sz[1])*uv+Sz[2])*uv+Sz[3])*uv+Sz[4])*uv+Sz[5])*uv
+Sz[6])*uv+Sz[7])*uv+Sz[8]);
qv:=zv/pv
end
else begin
while uv<2.0 do begin
zv:=zv/uv;
uv:=uv+1.0
end;
if uv=2.0 then qv:=zv
else begin
uv:=uv-2.0;
pv:=((((((Pz[0]*uv+Pz[1])*uv+Pz[2])*uv+Pz[3])*uv+Pz[4])*uv+Pz[5])*uv
+Pz[6])*uv+Pz[7];
qv:=(((((((Qz[0]*uv+Qz[1])*uv+Qz[2])*uv+Qz[3])*uv+Qz[4])*uv+Qz[5])*uv
+Qz[6])*uv+Qz[7])*uv+Qz[8];
qv:=zv*pv/qv
end
end
end;
gamma:=qv
END;

{continued fraction expansion for incomplete beta integral (form 1)}
FUNCTION confr_a(a,b,x: extended): extended;
BEGIN
k1:=a;
k2:=a+b;
k3:=a;
k4:=a+1.0;
k5:=1.0;
k6:=b-1.0;
k7:=k4;
k8:=a+2.0;
m2:=0.0;
m4:=1.0;
m1:=1.0;
m3:=1.0;
ad:=1.0;
rz:=1.0;
ts:=3.0*MCHP;
fg:=0;
for nz:=1 to 434 do begin
xk:=-(x*k1*k2)/(k3*k4);
pa:=m1+m2*xk;
qa:=m3+m4*xk;
m2:=m1;
m1:=pa;
m4:=m3;
m3:=qa;
xk:=(x*k5*k6)/(k7*k8);
pa:=m1+m2*xk;
qa:=m3+m4*xk;
m2:=m1;
m1:=pa;
m4:=m3;
m3:=qa;
if qa<>0.0 then rz:=pa/qa;
if rz<>0.0 then begin
tv:=abs((ad-rz)/rz);
ad:=rz
end
else tv:=1.0;
if tv<ts then begin
fg:=1;
break
end;
k1:=k1+1.0;
k2:=k2+1.0;
k3:=k3+2.0;
k4:=k4+2.0;
k5:=k5+1.0;
k6:=k6-1.0;
k7:=k7+2.0;
k8:=k8+2.0;
bp:=abs(pa); bq:=abs(qa);
if bp+bq>NBIG then begin
m1:=m1*BINV;
m2:=m2*BINV;
m3:=m3*BINV;
m4:=m4*BINV
end;
if (bp<BINV) or (bq<BINV) then begin 
m1:=m1*NBIG;
m2:=m2*NBIG;
m3:=m3*NBIG;
m4:=m4*NBIG
end
end;
if fg=0 then writeln('betai confr_a: loss of precision');
confr_a:=ad
END;

{continued fraction expansion for incomplete beta integral (form 2)}
FUNCTION confr_b(a,b,x: extended): extended;
BEGIN
k1:=a;
k2:=b-1.0;
k3:=a;
k4:=a+1.0;
k5:=1.0;
k6:=a+b;  
k7:=a+1;
k8:=a+2.0;
m2:=0.0;
m4:=1.0;
m1:=1.0;
m3:=1.0;
ad:=1.0;
rz:=1.0;
zk:=x/(1.0-x);
ts:=3.0*MCHP;
fg:=0;
for nz:=1 to 434 do begin
xk:=-(zk*k1*k2)/(k3*k4);
pa:=m1+m2*xk;
qa:=m3+m4*xk;
m2:=m1;
m1:=pa;
m4:=m3;
m3:=qa;
xk:=(zk*k5*k6)/(k7*k8);
pa:=m1+m2*xk;
qa:=m3+m4*xk;
m2:=m1;
m1:=pa;
m4:=m3;
m3:=qa;
if qa<>0.0 then rz:=pa/qa;
if rz<>0.0 then begin
tv:=abs((ad-rz)/rz);
ad:=rz
end
else tv:=1.0;
if tv<ts then begin
fg:=1;
break
end;
k1:=k1+1.0;
k2:=k2-1.0;
k3:=k3+2.0;
k4:=k4+2.0;
k5:=k5+1.0;
k6:=k6+1.0;
k7:=k7+2.0;
k8:=k8+2.0;
bp:=abs(pa); bq:=abs(qa);
if bp+bq>NBIG then begin
m1:=m1*BINV;
m2:=m2*BINV;
m3:=m3*BINV;
m4:=m4*BINV
end;
if (bp<BINV) or (bq<BINV) then begin 
m1:=m1*NBIG;
m2:=m2*NBIG;
m3:=m3*NBIG;
m4:=m4*NBIG
end
end;
if fg=0 then writeln('betai confr_b: loss of precision');
confr_b:=ad
END;

{power series expansion for incomplete beta integral}
FUNCTION powser(a,b,x: extended): extended;
BEGIN
ai:=1.0/a;
ua:=(1.0-b)*x;
va:=ua/(a+1.0);
t1:=va;
ta:=ua;
na:=2.0;
sa:=0.0;
zr:=MCHP*ai;
while abs(va)>zr do begin
ua:=(na-b)*x/na;
ta:=ta*ua;
va:=ta/(a+na);
sa:=sa+va;
na:=na+1
end;
sa:=sa+t1+ai;
ua:=a*ln(x);
if ((a+b)<OGAM) and (abs(ua)<MXGL) then begin
ta:=gamma(a+b)/(gamma(a)*gamma(b));
sa:=sa*ta*exp(a*ln(x))
end
else begin
ta:=lngamma(a+b)-lngamma(a)-lngamma(b)+ua+ln(sa);
if ta<MIGL then sa:=0.0 else sa:=exp(ta)
end;
powser:=sa
END;

{incomplete beta integral}
FUNCTION betai(a,b,x: extended): extended;
BEGIN
if (a<=0.0) or (b<=0.0) then begin
writeln('wrong betai argument, a and b must be > 0');
betai:=0.0;
exit
end;
if (x<=0.0) or (x>=1.0) then begin
if x=0.0 then betai:=0.0
else if x=1.0 then betai:=1.0
else begin
writeln('wrong betai argument, x must be >= 0 and <= 1'); 
betai:=0.0
end;
exit
end;
flag:=0;
if (b*x<=1.0) and (x<=0.95) then t:=powser(a,b,x)
else begin
w:=1.0-x;
if x>a/(a+b) then begin
flag:=1;
ac:=b;
bc:=a;
xc:=x;
xv:=w
end
else begin
ac:=a;
bc:=b;
xv:=x;
xc:=w
end;
if (flag=1) and (bc*xv<=1.0) and (xv<=0.95) then t:=powser(ac,bc,xv)
else begin
y:=xv*(ac+bc-2.0)-(ac-1.0);
if y<0.0 then w:=confr_a(ac,bc,xv) else w:=confr_b(ac,bc,xv)/xc;
y:=ac*ln(xv);
t:=bc*ln(xc);
if (ac+bc<OGAM) and (abs(y)<MXGL) and (abs(t)<MXGL) then begin
t:=exp(bc*ln(xc)+ac*ln(xv));
t:=(t/ac)*w;
t:=t*gamma(ac+bc)/(gamma(ac)*gamma(bc))
end
else begin
y:=y+t+lngamma(ac+bc)-lngamma(ac)-lngamma(bc);
y:=y+ln(w/ac);
if y<MIGL then t:=0.0 else t:=exp(y)
end
end
end;
if flag=1 then
if t<=MCHP then t:=1.0-MCHP else t:=1.0-t;
betai:=t
END;

{inverse incomplete beta integral}
FUNCTION ibetai(a,b,x: extended): extended;
LABEL dict,cmpt,quit;
BEGIN
if (a<=0.0) or (b<=0.0) then begin
writeln('wrong ibetai argument, a and b must be > 0');
ibetai:=0.0;
exit
end;
if (x<0.0) or (x>1.0) then begin
writeln('wrong ibetai argument, x must be >= 0 and <= 1');
ibetai:=0.0;
exit
end;
if x=0.0 then begin
ibetai:=0.0;
exit
end;
if x=1.0 then begin
ibetai:=1.0;
exit
end;
xq:=0.0;
yr:=0.0;
xs:=1.0;
yk:=1.0;
if (a<=1.0) or (b<=1.0) then begin
dr:=1.0e-7;
rf:=0;
av:=a;
bv:=b;
y0:=x;
xz:=av/(av+bv);
yz:=betai(av,bv,xz);
nf:=0;
goto dict
end
else begin
nf:=0;
dr:=1.0e-4
end;
yq:=-ilnor(x);
if x>0.5 then begin
rf:=1;
av:=b;
bv:=a;
y0:=1.0-x;
yq:=-yq
end
else begin
rf:=0;
av:=a;
bv:=b;
y0:=x
end;
lm:=(yq*yq-3.0)/6.0;
xz:=4.0/(1.0/(av-0.5)+1.0/(bv-0.5));
dz:=2.0*(yq*sqrt(xz+lm)/xz-(0.5/(bv-0.5)-0.5/(av-0.5))
*(lm+(5.0/6.0)-2.0/(3.0*xz)));
if dz<MIGL then begin
writeln('ibetai underflow');
ibetai:=0.0;
exit
end;
xz:=av/(av+bv*exp(dz));
yz:=betai(av,bv,xz);
yq:=(yz-y0)/y0;
if abs(yq)<0.2 then goto cmpt;
dict:
dq:=0;
di:=0.5;
for ir:=0 to 400 do begin
if ir>0 then begin
xz:=xq+di*(xs-xq);
if xz=1.0 then xz:=1.0-MCHP;
if xz=0.0 then begin
di:=0.5;
xz:=xq+di*(xs-xq);
if xz=0.0 then begin
writeln('ibetai underflow');
ibetai:=0.0;
exit
end
end;
yz:=betai(av,bv,xz);
yq:=(xs-xq)/(xs+xq);
if abs(yq)<dr then goto cmpt;
yq:=(yz-y0)/y0;
if abs(yq)<dr then goto cmpt
end;
if yz<y0 then begin
xq:=xz;
yr:=yz;
if dq<0 then begin
dq:=0;
di:=0.5
end
else if dq>3 then di:=di*(2.0-di)
else if dq>1 then di:=0.5*di+0.5
else  di:=(y0-yz)/(yk-yr);
inc(dq);
if xq>0.95 then begin
if rf=1 then begin
rf:=0; 
av:=a;
bv:=b;
y0:=x
end
else begin
rf:=1;
av:=b;
bv:=a;
y0:=1.0-x
end;
xz:=1.0-xz;
yz:=betai(av,bv,xz);
xq:=0.0;
yr:=0.0;
xs:=1.0;
yk:=1.0;
goto dict
end
end
else begin
xs:=xz;
if (rf=1) and (xs<MCHP) then begin
xz:=0.0;
goto quit
end;
yk:=yz;
if dq>0 then begin
dq:=0;
di:=0.5
end
else if dq<-3 then di:=di*di
else if dq<-1 then di:=0.5*di
else di:=(yz-y0)/(yk-yr);
dec(dq)
end
end;
if xq>=1.0 then begin
xz:=1.0-MCHP;
goto quit
end;
if xz<=0.0 then begin
writeln('ibetai underflow');
ibetai:=0.0;
exit
end;
cmpt:
if nf>0 then goto quit;
nf:=1;
lm:=lngamma(av+bv)-lngamma(av)-lngamma(bv);
for ir:=0 to 14 do begin
if ir>0 then yz:=betai(av,bv,xz);
if yz<yr then begin
xz:=xq;
yz:=yr
end
else if yz>yk then begin
xz:=xs;
yz:=yk
end
else if yz<y0 then begin
xq:=xz;
yr:=yz
end
else begin
xs:=xz;
yk:=yz
end;
if (xz=1.0) or (xz=0.0) then break;
dz:=(av-1.0)*ln(xz)+(bv-1.0)*ln(1.0-xz)+lm;
if dz<MIGL then goto quit;
if dz>MXGL then break;
dz:=(yz-y0)/exp(dz);
xt:=xz-dz;
if xt<=xq then begin
yz:=(xz-xq)/(xs-xq);
xt:=xq+0.5*yz*(xz-xq);
if xt<=0.0 then break
end;
if xt>=xs then begin
yz:=(xs-xz)/(xs-xq);
xt:=xs-0.5*yz*(xs-xz);
if xt>=1.0 then break
end;
xz:=xt;
if abs(dz/xz)<128.0*MCHP then goto quit
end;
dr:=256.0*MCHP;
goto dict;
quit:
if rf>0 then
if xz<=MCHP then xz:=1.0-MCHP else xz:=1.0-xz;
ibetai:=xz
END;

{complementary error function}
FUNCTION erfc(a: extended): extended;
BEGIN
if a<UERF then yb:=2.0
else if a>URFC then yb:=0.0
else begin
xb:=abs(a);
if xb<1.0 then begin
xb:=a*a;
zb:=a*((((((Tk[0]*xb+Tk[1])*xb+Tk[2])*xb+Tk[3])*xb+Tk[4])*xb+Tk[5])*xb+Tk[6])/
((((((xb+Uk[0])*xb+Uk[1])*xb+Uk[2])*xb+Uk[3])*xb+Uk[4])*xb+Uk[5]);
yb:=1.0-zb
end
else begin
yb:=1.0/xb;
if xb<8.0 then begin
pb:=((((((((Ph[0]*yb+Ph[1])*yb+Ph[2])*yb+Ph[3])*yb+Ph[4])*yb+Ph[5])*yb
+Ph[6])*yb+Ph[7])*yb+Ph[8])*yb+Ph[9];
qb:=(((((((((yb+Qh[0])*yb+Qh[1])*yb+Qh[2])*yb+Qh[3])*yb+Qh[4])*yb+Qh[5])*yb
+Qh[6])*yb+Qh[7])*yb+Qh[8])*yb+Qh[9]
end
else begin
zb:=yb*yb;
pb:=yb*((((Rh[0]*zb+Rh[1])*zb+Rh[2])*zb+Rh[3])*zb+Rh[4]);
qb:=((((zb+Sh[0])*zb+Sh[1])*zb+Sh[2])*zb+Sh[3])*zb+Sh[4]
end;
yb:=exp(-a*a)*pb/qb;
if a<0 then yb:=2.0-yb
end
end;
erfc:=yb
END;

{inverse complementary error function}
FUNCTION ierfc(a: extended): extended;
BEGIN
if (a<=0.0) or (a>=2.0) then begin
writeln('wrong ierfc argument, must be > 0 and < 2');
ierfc:=0.0
end
else if a=1.0 then ierfc:=0.0
else ierfc:=ierf(1.0-a)
END;

{error function}
FUNCTION erf(a: extended): extended;
BEGIN
if a<UERF then yu:=-1.0
else if a>OERF then yu:=1.0
else if abs(a)>1.0 then yu:=1.0-erfc(a)
else begin
zu:=a*a;
yu:=a*((((((Tk[0]*zu+Tk[1])*zu+Tk[2])*zu+Tk[3])*zu+Tk[4])*zu+Tk[5])*zu+Tk[6])/
((((((zu+Uk[0])*zu+Uk[1])*zu+Uk[2])*zu+Uk[3])*zu+Uk[4])*zu+Uk[5])
end;
erf:=yu
END;

{inverse error function}
FUNCTION ierf(a: extended): extended;
BEGIN
if (a<=-1.0) or (a>=1.0) then begin
writeln('wrong ierf argument, must be > -1 and < 1');
ierf:=0.0
end
else if a=0.0 then ierf:=0.0
else ierf:=-ilnor((1.0-a)/2.0)/SQT2
END;

{cumulative probability of the standard normal distribution}
FUNCTION phi(a: extended): extended;
BEGIN
if a<0.0 then begin
writeln('wrong phi argument, must be >= 0');
phi:=0.0
end
else if a=0.0 then phi:=0.0
else if a>OPHI then phi:=0.5
else phi:=0.5*erf(a/SQT2)
END;

{inverse cumulative probability of the standard normal distribution}
FUNCTION iphi(a: extended): extended;
BEGIN
if (a<0.0) or (a>=0.5) then begin
writeln('wrong iphi argument, must be >= 0 and < 0.5');
iphi:=0.0;
exit
end;
if a<=0.425 then begin
r:=0.180625-a*a;
m:=a*(((((((Aj[7]*r+Aj[6])*r+Aj[5])*r+Aj[4])*r+Aj[3])*r+Aj[2])*r
+Aj[1])*r+Aj[0])/(((((((Bj[7]*r+Bj[6])*r+Bj[5])*r+Bj[4])*r+Bj[3])*r
+Bj[2])*r+Bj[1])*r+1.0)
end
else begin
r:=sqrt(-ln(0.5-a));
if r<=5.0 then begin
r:=r-1.6;
m:=(((((((Cj[7]*r+Cj[6])*r+Cj[5])*r+Cj[4])*r+Cj[3])*r+Cj[2])*r+
Cj[1])*r+Cj[0])/(((((((Dj[7]*r+Dj[6])*r+Dj[5])*r+Dj[4])*r+Dj[3])*r
+Dj[2])*r+Dj[1])*r+1.0)
end
else begin
r:=r-5.0;
m:=(((((((Ej[7]*r+Ej[6])*r+Ej[5])*r+Ej[4])*r+Ej[3])*r+Ej[2])*r
+Ej[1])*r+Ej[0])/(((((((Fj[7]*r+Fj[6])*r+Fj[5])*r+Fj[4])*r+Fj[3])*r+
Fj[2])*r+Fj[1])*r+1.0)
end
end;
if m<0 then iphi:=-m else iphi:=m
END;

{lower tail normal cumulative distribution}
FUNCTION lnor(a: extended): extended;
BEGIN
if a<ULNR then ya:=0.0
else if a>OLNR then ya:=1.0
else if abs(a)<1.0 then ya:=0.5+0.5*erf(a*ISQ2)
else begin
ya:=0.5*erfc(abs(a)*ISQ2);
if a>0.0 then ya:=1.0-ya
end;
lnor:=ya
END;

{inverse lower tail normal cumulative distribution}
FUNCTION ilnor(a: extended): extended;
BEGIN
if (a<=0.0) or (a>=1.0) then begin
writeln('wrong ilnor argument, must be > 0 and < 1');
xe:=0.0
end
else begin
fl:=1;
ye:=a;
if ye>XPML then begin
ye:=1.0-ye;
fl:=0
end;
if ye>XPM2 then begin
ye:=ye-0.5;
y2:=ye*ye;
xe:=(((((((U0[0]*y2+U0[1])*y2+U0[2])*y2+U0[3])*y2+U0[4])*y2+U0[5])*y2
+u0[6])*y2+U0[7])/
(((((((y2+V0[0])*y2+V0[1])*y2+V0[2])*y2+V0[3])*y2+V0[4])*y2+V0[5])*y2+
V0[6]);
xe:=ye*(1.0+y2*xe)*S2PI
end
else begin
xe:=sqrt(-2*ln(ye));
x2:=xe-ln(xe)/xe;
ze:=1.0/xe;
if xe<8.0 then
z2:=ze*(((((((((U1[0]*ze+U1[1])*ze+U1[2])*ze+U1[3])*ze+U1[4])*ze+U1[5])*ze+
U1[6])*ze+U1[7])*ze+U1[8])*ze+U1[9])/
(((((((((ze+V1[0])*ze+V1[1])*ze+V1[2])*ze+V1[3])*ze+V1[4])*ze+V1[5])*ze+
V1[6])*ze+V1[7])*ze+V1[8])
else if xe<32.0 then
z2:=ze*(((((((U2[0]*ze+U2[1])*ze+U2[2])*ze+U2[3])*ze+U2[4])*ze+U2[5])*ze+
U2[6])*ze+U2[7])/
(((((((ze+V2[0])*ze+V2[1])*ze+V2[2])*ze+V2[3])*ze+V2[4])*ze+V2[5])*ze+V2[6])
else
z2:=ze*(((((((U3[0]*ze+U3[1])*ze+U3[2])*ze+U3[3])*ze+U3[4])*ze+U3[5])*ze+
U3[6])*ze+U3[7])/
(((((((ze+V3[0])*ze+V3[1])*ze+V3[2])*ze+V3[3])*ze+V3[4])*ze+V3[5])*ze+V3[6]);
xe:=x2-z2;
if fl=1 then xe:=-xe
end
end;
ilnor:=xe
END;

{upper tail normal cumulative distribution}
FUNCTION unor(a: extended): extended;
BEGIN
if a<OUNR then unor:=1.0
else if a>UUNR then unor:=0.0
else unor:=1.0-lnor(a)
END;

{inverse upper tail normal cumulative distribution}
FUNCTION iunor(a: extended): extended;
BEGIN
if (a<=0.0) or (a>=1.0) then begin
writeln('wrong iunor argument, must be > 0 and < 1');
iunor:=0.0
end
else iunor:=ilnor(1.0-a)
END;

{upper tail cumulative probability of the chi-square distribution with v
degrees of freedom}
FUNCTION chsqc(a: extended; v: longword): extended;
BEGIN
if a<0 then begin
writeln('wrong chsqc argument, must be >= 0');
chsqc:=0.0
end
else if (a=0.0) or (v=0) then chsqc:=1.0
else if v=2 then chsqc:=exp(-0.5*a)
else begin
vx:=v;
chsqc:=cgami(vx/2.0,a/2.0)
end
END;

{inverse upper tail cumulative probability of the chi-square distribution
with v degrees of freedom}
FUNCTION ichsqc(a: extended; v: longword): extended;
BEGIN
if (a<=0.0) or (a>1.0) then begin
writeln('wrong ichsqc argument, must be > 0 and <= 1');
ichsqc:=0.0
end
else if a<MICG then begin
writeln('ichsqc overflow, a must be > ',MICG);
ichsqc:=0.0
end
else if v=0 then begin
writeln('wrong ichsqc n.d.f., must be >= 1');
ichsqc:=0.0
end
else if a=1.0 then ichsqc:=0.0 
else begin
df:=v;
ichsqc:=2.0*icgami(0.5*df,a)
end
END;

{lower cumulative probability of the chi-square distribution with v degrees
of freedom}
FUNCTION chsq(a: extended; v: longword): extended;
BEGIN
if a<0 then begin
writeln('wrong chsq argument, must be >= 0');
chsq:=0.0
end
else if (a=0.0) or (v=0) then chsq:=0.0
else if v=2 then chsq:=1.0-exp(-0.5*a)
else begin
vx:=v;
chsq:=gami(vx/2.0,a/2.0)
end
END;

{inverse lower cumulative probability of the chi-square distribution with v
degrees of freedom}
FUNCTION ichsq(a: extended; v: longword): extended;
BEGIN
if (a<0.0) or (a>=1.0) then begin
writeln('wrong ichsq argument, must be >= 0 and < 1');
ichsq:=0.0
end
else if v=0 then begin
writeln('wrong ichsq n.d.f., must be >= 1');
ichsq:=0.0
end
else if a=0.0 then ichsq:=0.0 
else begin
df:=v;
ichsq:=2.0*icgami(0.5*df,1.0-a)
end
END;

END.
