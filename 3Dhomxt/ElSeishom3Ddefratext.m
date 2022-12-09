% SPDX-License-Identifier: CC0-1.0

% defining all the variables
clear all;
clc;
close all;
% predetermine fontsize
fs=18;
% Setting the parameters for a source wave
nf              =   128;
mf              =   4; % multiplication factor to have sinc-interpolation in time-domain
% Define the center frequency of the source waveform
fc              =   40;
% Define frequency axis
freq            =   linspace(0,4*fc,nf);
freq(1)         =   1e-6; % To prevent unnecessary problems upon division by zero. Can't be too close to zero, though, 
                         % since some variables either blow up or vanish in comparison to the rest of the entries in the frequency vector.
ss              =   2i*pi*freq;
% compute the frequency step size
df              =   freq(3)-freq(2);
dt              =   1/(mf*nf*df);
t               =   linspace(-mf*nf/2,mf*nf/2-1,mf*nf)*dt;


% Define Medium Parameters (frequency-independent)
rhof            =   1.0e3; % fluid density
rho             =   2.7e3;
Gfr             =   9.0e9; % shear modulus of the framework of grains
eta             =   1.0e-3; % fluid viscosity
k0              =   1.3e-12; % medium permeability (static)
Kfr             =   4.0e9;
Ks              =   4.0e10;
Kf              =   2.2e9;
Concentr        =   1.0e-4;
bplus           =   3.0e11;
bmin            =   3.0e11;
porosity        =   0.3;
epsilonRF       =   80.0;
epsilonRS       =   4.0;
alpha_inf       =   3.0;
pH              =   7.0;
similaritypar   =   8;

c0              =   299792458; %// velocity of light in free-space
mu0             =   4.0e-7*pi; %// free-space magnetic permeability
epsilon0        =   1.0/(mu0*c0*c0); %// free-space electric permittivity
e               =   1.602e-19; %// elementary charge
z_1             =   1; %// ion valences
z_1c            =   -1; %// valency of the conjugate ion
NA              =   6.022e23; %// Avogadro's constant [mol^{-1}]
kb              =   1.381e-23;  %// Boltzmann constant
T               =   295.0; %// Temperature in Kelvin

epsilonR        =   (porosity/alpha_inf)*(epsilonRF-epsilonRS)+epsilonRS ;
omegac          =   (porosity*eta)/(alpha_inf*k0*rhof);% critical frequency
zetap           =   8e-3+26e-3*log10(Concentr);% zeta potential, empirical relation (experimental studies, Pride&Morgan 1991)
L0              =   -(porosity*epsilon0*epsilonRF*zetap)/(alpha_inf*eta); % static coupling coefficient
L               =   L0; %because omega << omega_c (critical frequency)
N               =   10e3*Concentr*NA*abs(z_1c); % bulk-ionic concentration (of species i, in this case only 2 species; binary symmetric electrolyte)
sigmaF          =   ((e*z_1)^2)*N*(bplus+bmin); %! conductivity of the pore-fluid phase (simplified cond. pure electrolyte)
sigmaE          =   (porosity*sigmaF)/(alpha_inf); %! bulk electric conductivity (freq-(in)dependent)
% %!sigma_E = (poros*sigmaF)/(alpha_inf)*(1+(2*(sigma_em+sigma_os))/(sigmaF*Lambda)); ! = complete version: bulk electric conductivity (freq-dependent */
rhoB            =   (1.0-porosity)*rho+porosity*rhof; % effective density of the fluid (in relative motion)
epsilon         =   epsilon0*epsilonR;
Delta           =   Kf*((1-porosity)*Ks-Kfr)/(porosity*(Ks)^2); %! combination of the frequency-independent compression moduli
Kg              =   (Kfr+porosity*Kf+(1+porosity)*Ks*Delta)/(1+Delta); %! Gassmann's bulk modulus
C               =   (Kf+Ks*Delta)/(1+Delta);
M               =    Kf/(porosity*(1+Delta)); %! the elastic media parameter S is defined as Kg-2/3Gfr-C**2/ M
H               =   Kg + 4.0*Gfr/3.0;
S               =   Kg-((2.0/3.0)*Gfr)-((C*C)/M);
sigmaM          =   0.0; % magnetic conductivity
Kc              =   (S+2.*Gfr);
% defining the grid
nx=128;
xx1=linspace(-350,350,nx);
xx2=xx1;
[x1,x2,s]=meshgrid(xx1,xx2,ss);
% fixed vertical source-receiver distance
x3=50;
% distance and unit vectors
R = sqrt(x1.^2 + x2.^2 + x3.^2);
th1=x1./R;
th2=x2./R;
th3=x3./R;
% frequency dependent parameters
k               =   k0./(sqrt(1+4.*s./(similaritypar.*omegac))+s./omegac);% frequency-dependent dynamic permeability
rhoE            =   eta./(s.*k); % effective density
rhoc            =   rhoB-((rhof.*rhof)./rhoE); % Complex Density rhoC
zeta            =   sigmaM+s.*mu0; 
etae            =   sigmaE+s.*epsilon;
varsigma        =   etae - s.*rhoE.*L.^2;
chi             =   s.*rhof*L;
% source signature
wav=-2*sqrt(1/pi)*(s/(2*pi*fc)).^2.*exp((s/(2*pi*fc)).^2)/fc;
% spherical wavenumbers of eqs B-9, B-10, B-12, B-13.
ypf     =   s./sqrt(2).*sqrt((rhoB.*M-2.*rhof.*C+rhoE.*H.*etae./varsigma)./(H.*M-C.^2) - sqrt(((rhoB.*M-2.*rhof.*C+rhoE.*H.*etae./varsigma)./(H.*M-C.^2)).^2 + 4.*(rhof.^2-rhoB.*rhoE.*etae./varsigma)./(H.*M-C.^2)));
yps     =   s./sqrt(2).*sqrt((rhoB.*M-2.*rhof.*C+rhoE.*H.*etae./varsigma)./(H.*M-C.^2) + sqrt(((rhoB.*M-2.*rhof.*C+rhoE.*H.*etae./varsigma)./(H.*M-C.^2)).^2 + 4.*(rhof.^2-rhoB.*rhoE.*etae./varsigma)./(H.*M-C.^2)));
ys      =   1/sqrt(2).*sqrt(s.^2.*rhoc./Gfr+zeta.*etae - sqrt((s.^2.*rhoc./Gfr-zeta.*etae).^2-4.*s.^3.*zeta.*(rhof.*L).^2./Gfr));
yem     =   1/sqrt(2).*sqrt(s.^2.*rhoc./Gfr+zeta.*etae + sqrt((s.^2.*rhoc./Gfr-zeta.*etae).^2-4.*s.^3.*zeta.*(rhof.*L).^2./Gfr));
% recurring coefficient factors
dems    =   1./(yem.^2-ys.^2);            
dpspf   =   1./(yps.^2-ypf.^2);     

% deformation rate, h_11 only

% Green's functions of eq C-10
GS = exp(-ys.*R)./(4*pi*R);
GEM = exp(-yem.*R)./(4*pi*R);
GPs = exp(-yps.*R)./(4*pi*R);
GPf = exp(-ypf.*R)./(4*pi*R);
% grad (equation C-15)
[GS1,GS2,GS3]=Greeni(ys,th1,th2,th3,R);
[GEM1,GEM2,GEM3]=Greeni(yem,th1,th2,th3,R);
[GPf1,GPf2,GPf3]=Greeni(ypf,th1,th2,th3,R);
[GPs1,GPs2,GPs3]=Greeni(yps,th1,th2,th3,R);
% grad-grad divided by gamma^2 (equation C-16)
[GS11,GS21,GS31,GS22,GS32,GS33]=Greenij(ys,th1,th2,th3,R);
[GEM11,GEM21,GEM31,GEM22,GEM32,GEM33]=Greenij(yem,th1,th2,th3,R);
[GPf11,GPf21,GPf31,GPf22,GPf32,GPf33]=Greenij(ypf,th1,th2,th3,R);
[GPs11,GPs21,GPs31,GPs22,GPs32,GPs33]=Greenij(yps,th1,th2,th3,R);
% grad-grad - delta_{ij} is coded on the fly
% grad-grad-grad divided by gamma^2 without delta_{ij}-grad in (equation C-18)
% because then it has high-degree of symmetry and only 10 out of 27 are
% needed
[GS111,GS211,GS311,GS221,GS321,GS331,GS222,GS322,GS332,GS333]=Greenijk(ys,th1,th2,th3,R);
[GEM111,GEM211,GEM311,GEM221,GEM321,GEM331,GEM222,GEM322,GEM332,GEM333]=Greenijk(yem,th1,th2,th3,R);
[GPf111,GPf211,GPf311,GPf221,GPf321,GPf331,GPf222,GPf322,GPf332,GPf333]=Greenijk(ypf,th1,th2,th3,R);
[GPs111,GPs211,GPs311,GPs221,GPs321,GPs331,GPs222,GPs322,GPs332,GPs333]=Greenijk(yps,th1,th2,th3,R);
% grad-grad-grad-grad/gam^4 again for symmetry reasons so now only 15 out
% of 81 are needed
[GS1111,GS2111,GS3111,GS2211,GS3211,GS3311,GS2221,GS3221,GS3321,...
    GS3331,GS2222,GS3222,GS3322,GS3332,GS3333]=Greenijkl(ys,th1,th2,th3,R);
[GEM1111,GEM2111,GEM3111,GEM2211,GEM3211,GEM3311,GEM2221,GEM3221,GEM3321,...
    GEM3331,GEM2222,GEM3222,GEM3322,GEM3332,GEM3333]=Greenijkl(ys,th1,th2,th3,R);
[GPf1111,GPf2111,GPf3111,GPf2211,GPf3211,GPf3311,GPf2221,GPf3221,GPf3321,...
    GPf3331,GPf2222,GPf3222,GPf3322,GPf3332,GPf3333]=Greenijkl(ys,th1,th2,th3,R);
[GPs1111,GPs2111,GPs3111,GPs2211,GPs3211,GPs3311,GPs2221,GPs3221,GPs3321,...
    GPs3331,GPs2222,GPs3222,GPs3322,GPs3332,GPs3333]=Greenijkl(ys,th1,th2,th3,R);

% particle velocity
KSvh = (zeta.*etae-ys.^2).*dems;
KEMvh = (zeta.*etae-yem.^2).*dems;
KPsvh = -2*Gfr*M*(s.^2.*rhoE.*etae./(M*varsigma)-yps.^2).*dpspf/(H*M-C^2);
KPfvh = -2*Gfr*M*(s.^2.*rhoE.*etae./(M*varsigma)-ypf.^2).*dpspf/(H*M-C^2);
NPsvh = -(s.^2.*(rhoE.*H.*etae./varsigma-rhof*C)/(H*M-C^2)-yps.^2).*dpspf;
NPfvh = -(s.^2.*(rhoE.*H.*etae./varsigma-rhof*C)/(H*M-C^2)-ypf.^2).*dpspf;
v111 = 2*(KSvh.*(GS111-GS1)-KEMvh.*(GEM111-GEM1))+KPfvh.*(GPf111-GPf1)...
    -KPsvh.*(GPs111-GPs1) + NPfvh.*GPf1-NPsvh.*GPs1;
v211 = 2*(KSvh.*GS211-KEMvh.*GEM211)+KPfvh.*(GPf211-GPf2)-KPsvh.*(GPs211-GPs2)...
    + NPfvh.*GPf2-NPsvh.*GPs2;
v311 = 2*(KSvh.*GS311-KEMvh.*GEM311)+KPfvh.*(GPf311-GPf3)-KPsvh.*(GPs311-GPs3)...
    + NPfvh.*GPf3-NPsvh.*GPs3;
v121 = KSvh.*(2*GS211-GS2)-KEMvh.*(2*GEM211-GEM2)+KPfvh.*GPf211-KPsvh.*GPs211;
v221 = KSvh.*(2*GS221-GS1)-KEMvh.*(2*GEM221-GEM1) + KPfvh.*GPf221-KPsvh.*GPs211;
v321 = 2*KSvh.*GS321-2*KEMvh.*GEM321+KPfvh.*GPf321-KPsvh.*GPs321;
v122 = 2*KSvh.*GS221-2*KEMvh.*GEM221+KPfvh.*(GPf221-GPf1)-KPsvh.*(GPs221-GPs1) ...
      + NPfvh.*GPf1-NPsvh.*GPs1;
v222 = 2*(KSvh.*(GS222-GS2)-KEMvh.*(GEM222-GEM2))+KPfvh.*(GPf222-GPf2)...
    -KPsvh.*(GPs222-GPs2) + NPfvh.*GPf2-NPsvh.*GPs2;
v322 = 2*KSvh.*GS322-2*KEMvh.*GEM322+KPfvh.*(GPf322-GPf3)-KPsvh.*(GPs322-GPs3)...
    + NPfvh.*GPf3-NPsvh.*GPs3;
v113 = KSvh.*(2*GS311-GS3)-KEMvh.*(2*GEM311-GEM3)+KPfvh.*GPf311-KPsvh.*GPs311;
% v213 = v321;
v313 = KSvh.*(2*GS331-GS1)-KEMvh.*(2*GEM331-GEM1) + KPfvh.*GPf331-KPsvh.*GPs331;
% v123 = v321;
v223 = KSvh.*(2*GS322-GS3)-KEMvh.*(2*GEM322-GEM3)+KPfvh.*GPf322-KPsvh.*GPs322;
v323 = KSvh.*(2*GS332-GS2)-KEMvh.*(2*GEM332-GEM2) + KPfvh.*GPf332-KPsvh.*GPs332;
v133 = 2*KSvh.*GS331-2*KEMvh.*GEM331+KPfvh.*(GPf331-GPf1)-KPsvh.*(GPs331-GPs1) ...
      + NPfvh.*GPf1-NPsvh.*GPs1;
v233 = 2*KSvh.*GS332-2*KEMvh.*GEM322+KPfvh.*(GPf332-GPf2)-KPsvh.*(GPs332-GPs2)...
    + NPfvh.*GPf-NPsvh.*GPs2;
v333 = 2*(KSvh.*(GS333-GS3)-KEMvh.*(GEM333-GEM3))+KPfvh.*(GPf333-GPf3)-KPsvh.*(GPs333-GPs3)...
    + NPfvh.*GPf3-NPsvh.*GPs3;

% electric field
KSeh = chi.*zeta.*dems;
KEMeh = KSeh;
KPseh = -2*s.*rhoE.*L.*C.*Gfr.*(s.^2.*rhof/C-yps.^2).*dpspf./((H*M-C^2).*varsigma);
KPfeh = -2*s.*rhoE.*L.*C.*Gfr.*(s.^2.*rhof/C-ypf.^2).*dpspf./((H*M-C^2).*varsigma);
NPseh = -s.^3.*rhoE.*L.*(rhof*H-rhoB*C).*dpspf./((H*M-C^2).*varsigma);
NPfeh = NPseh;
E111 = 2*(KSeh.*(GS111-GS1)-KEMeh.*(GEM111-GEM1))+KPfeh.*(GPf111-GPf1)...
    -KPseh.*(GPs111-GPs1) + NPfeh.*GPf1-NPseh.*GPs1;
E211 = 2*(KSeh.*GS211-KEMeh.*GEM211)+KPfeh.*(GPf211-GPf2)-KPseh.*(GPs211-GPs2)...
    + NPfeh.*GPf2-NPseh.*GPs2;
E311 = 2*(KSeh.*GS311-KEMeh.*GEM311)+KPfeh.*(GPf311-GPf3)-KPseh.*(GPs311-GPs3)...
    + NPfeh.*GPf3-NPseh.*GPs3;
E121 = KSeh.*(2*GS211-GS2)-KEMeh.*(2*GEM211-GEM2)+KPfeh.*GPf211-KPseh.*GPs211;
E221 = KSeh.*(2*GS221-GS1)-KEMeh.*(2*GEM221-GEM1) + KPfeh.*GPf221-KPseh.*GPs211;
E321 = 2*KSeh.*GS321-2*KEMeh.*GEM321+KPfeh.*GPf321-KPseh.*GPs321;
E122 = 2*KSeh.*GS221-2*KEMeh.*GEM221+KPfeh.*(GPf221-GPf1)-KPseh.*(GPs221-GPs1) ...
      + NPfeh.*GPf1-NPseh.*GPs1;
E222 = 2*(KSeh.*(GS222-GS2)-KEMeh.*(GEM222-GEM2))+KPfeh.*(GPf222-GPf2)...
    -KPseh.*(GPs222-GPs2) + NPfeh.*GPf2-NPseh.*GPs2;
E322 = 2*KSeh.*GS322-2*KEMeh.*GEM322+KPfeh.*(GPf322-GPf3)-KPseh.*(GPs322-GPs3)...
    + NPfeh.*GPf3-NPseh.*GPs3;
e113 = KSeh.*(2*GS311-GS3)-KEMeh.*(2*GEM311-GEM3)+KPfeh.*GPf311-KPseh.*GPs311;
% E213 = E321;
E313 = KSeh.*(2*GS331-GS1)-KEMeh.*(2*GEM331-GEM1) + KPfeh.*GPf331-KPseh.*GPs331;
% E123 = E321;
E223 = KSeh.*(2*GS322-GS3)-KEMeh.*(2*GEM322-GEM3)+KPfeh.*GPf322-KPseh.*GPs322;
E323 = KSeh.*(2*GS332-GS2)-KEMeh.*(2*GEM332-GEM2) + KPfeh.*GPf332-KPseh.*GPs332;
E133 = 2*KSeh.*GS331-2*KEMeh.*GEM331+KPfeh.*(GPf331-GPf1)-KPseh.*(GPs331-GPs1) ...
      + NPfeh.*GPf1-NPseh.*GPs1;
E233 = 2*KSeh.*GS332-2*KEMeh.*GEM322+KPfeh.*(GPf332-GPf2)-KPseh.*(GPs332-GPs2)...
    + NPfeh.*GPf-NPseh.*GPs2;
E333 = 2*(KSeh.*(GS333-GS3)-KEMeh.*(GEM333-GEM3))+KPfeh.*(GPf333-GPf3)-KPseh.*(GPs333-GPs3)...
    + NPfeh.*GPf3-NPseh.*GPs3;
% filtration velocity
KSwh = rhof.*(zeta.*varsigma-ys.^2).*dems./rhoE;
KEMwh = rhof.*(zeta.*varsigma-yem.^2).*dems./rhoE;
KPswh = -2*Gfr*C.*(s.^2*rhof/C-yps.^2).*dpspf./(H*M-C^2);
KPfwh = -2*Gfr*C.*(s.^2*rhof/C-ypf.^2).*dpspf./(H*M-C^2);
NPfwh = -s.^2*(rhof*H-rhoB*C).*dpspf./(H*M-C^2);
NPswh = NPfwh;
w111 = 2*(KSwh.*(GS111-GS1)-KEMwh.*(GEM111-GEM1))+KPfwh.*(GPf111-GPf1)...
    -KPswh.*(GPs111-GPs1) + NPfwh.*GPf1-NPswh.*GPs1;
w211 = 2*(KSwh.*GS211-KEMwh.*GEM211)+KPfwh.*(GPf211-GPf2)-KPswh.*(GPs211-GPs2)...
    + NPfwh.*GPf2-NPswh.*GPs2;
w311 = 2*(KSwh.*GS311-KEMwh.*GEM311)+KPfwh.*(GPf311-GPf3)-KPswh.*(GPs311-GPs3)...
    + NPfwh.*GPf3-NPswh.*GPs3;
w121 = KSwh.*(2*GS211-GS2)-KEMwh.*(2*GEM211-GEM2)+KPfwh.*GPf211-KPswh.*GPs211;
w221 = KSwh.*(2*GS221-GS1)-KEMwh.*(2*GEM221-GEM1) + KPfwh.*GPf221-KPswh.*GPs211;
w321 = 2*KSwh.*GS321-2*KEMwh.*GEM321+KPfwh.*GPf321-KPswh.*GPs321;
w122 = 2*KSwh.*GS221-2*KEMwh.*GEM221+KPfwh.*(GPf221-GPf1)-KPswh.*(GPs221-GPs1) ...
      + NPfwh.*GPf1-NPswh.*GPs1;
w222 = 2*(KSwh.*(GS222-GS2)-KEMwh.*(GEM222-GEM2))+KPfwh.*(GPf222-GPf2)...
    -KPswh.*(GPs222-GPs2) + NPfwh.*GPf2-NPswh.*GPs2;
w322 = 2*KSwh.*GS322-2*KEMwh.*GEM322+KPfwh.*(GPf322-GPf3)-KPswh.*(GPs322-GPs3)...
    + NPfwh.*GPf3-NPswh.*GPs3;
w113 = KSwh.*(2*GS311-GS3)-KEMwh.*(2*GEM311-GEM3)+KPfwh.*GPf311-KPswh.*GPs311;
% w213 = w321;
w313 = KSwh.*(2*GS331-GS1)-KEMwh.*(2*GEM331-GEM1) + KPfwh.*GPf331-KPswh.*GPs331;
% w123 = w321;
w223 = KSwh.*(2*GS322-GS3)-KEMwh.*(2*GEM322-GEM3)+KPfwh.*GPf322-KPswh.*GPs322;
w323 = KSwh.*(2*GS332-GS2)-KEMwh.*(2*GEM332-GEM2) + KPfwh.*GPf332-KPswh.*GPs332;
w133 = 2*KSwh.*GS331-2*KEMwh.*GEM331+KPfwh.*(GPf331-GPf1)-KPswh.*(GPs331-GPs1) ...
      + NPfwh.*GPf1-NPswh.*GPs1;
w233 = 2*KSwh.*GS332-2*KEMwh.*GEM322+KPfwh.*(GPf332-GPf2)-KPswh.*(GPs332-GPs2)...
    + NPfwh.*GPf-NPswh.*GPs2;
w333 = 2*(KSwh.*(GS333-GS3)-KEMwh.*(GEM333-GEM3))+KPfwh.*(GPf333-GPf3)-KPswh.*(GPs333-GPs3)...
    + NPfwh.*GPf3-NPswh.*GPs3;
% acoustic pressure
KPPph = 2*Gfr*s.*(rhoE.*etae*C./varsigma-rhof*M).*dpspf/(H*M-C^2);
NPsph = -s.*(s.^2*C.*(rhoB*rhoE.*etae./varsigma-rhof^2)...
        ./(H*M-C^2)-rhof.*yps.^2).*dpspf;
NPfph = -s.*(s.^2*C.*(rhoB*rhoE.*etae./varsigma-rhof^2)...
        ./(H*M-C^2)-rhof.*ypf.^2).*dpspf;
p11 = KPPph.*(ypf.^2.*(GPf11-GPf)-yps.^2.*(GPs11-GPs)) + (NPfph.*GPf-NPsph.*GPs);
p21 = KPPph.*(ypf.^2.*GPf21-yps.^2.*GPs21);
p31 = KPPph.*(ypf.^2.*GPf31-yps.^2.*GPs31);
p22 = KPPph.*(ypf.^2.*(GPf22-GPf)-yps.^2.*(GPs22-GPs)) + (NPfph.*GPf-NPsph.*GPs);
p32 = KPPph.*(ypf.^2.*GPf32-yps.^2.*GPs32);
p33 = KPPph.*(ypf.^2.*(GPf33-GPf)-yps.^2.*(GPs33-GPs)) + (NPfph.*GPf-NPsph.*GPs);
% magnetic field
KSSmh = chi.*dems;
H111 = 0;
H211 = 2*KSSmh.*(ys.^2.*GS31-yem.^2.*GEM31);
H311 = -2*KSSmh.*(ys.^2.*GS21-yem.^2.*GEM21);

H121 = -KSSmh.*(ys.^2.*GS31-yem.^2.*GEM31);
H221 = KSSmh.*(ys.^2.*GS32-yem.^2.*GEM32);
H321 = KSSmh.*(ys.^2.*(GS11-GS22)-yem.^2.*(GEM11-GEM22));
H131 = KSSmh.*(ys.^2.*GS21-yem.^2.*GEM21);
H231 = KSSmh.*(ys.^2.*(GS11-GS33)-yem.^2.*(GEM11-GEM33));
H331 = -KSSmh.*(ys.^2.*GS32-yem.^2.*GEM32);
H122 = -2*KSSmh.*(ys.^2.*GS32-yem.^2.*GEM32);
H222 = 0;
H322 = 2*KSSmh.*(ys.^2.*GS21-yem.^2.*GEM21);
H132 = KSSmh.*(ys.^2.*(GS22-GS33)-yem.^2.*(GEM22-GEM33));
H232 = -KSSmh.*(ys.^2.*GS21-yem.^2.*GEM21);
H332 = KSSmh.*(ys.^2.*GS31-yem.^2.*GEM31);
H133 = 2*KSSmh.*(ys.^2.*GS32-yem.^2.*GEM32);
H233 = -2*KSSmh.*(ys.^2.*GS31-yem.^2.*GEM31);
H333 = 0;
% bulk stress field
KSth = Gfr.*(zeta.*etae-ys.^2).*dems./s;
KEMth = Gfr.*(zeta.*etae-yem.^2).*dems./s;
KPfth = 2*Gfr.*ypf.^2.*(ypf.^2-s.^2.*(rhoE*H.*etae./varsigma-rhof*C)./(H*M-C^2)).*dpspf./s;
KPsth = 2*Gfr.*yps.^2.*(yps.^2-s.^2.*(rhoE*H.*etae./varsigma-rhof*C)./(H*M-C^2)).*dpspf./s;
NPfth = H*ypf.^2.*(ypf.^2-s.^2.*(rhoE*H.*etae./varsigma-2*rhof*C+rhoB*C^2/H)...
        ./(H*M-C^2)).*dpspf./s;
NPsth = H*yps.^2.*(yps.^2-s.^2.*(rhoE*H.*etae./varsigma-2*rhof*C+rhoB*C^2/H)...
        ./(H*M-C^2)).*dpspf./s;
QPfth = 4*Gfr^2*M*ypf.^2.*(ypf.^2-s.^2.*rhoE.*etae./(M*varsigma)).*dpspf./(s.*(H*M-C^2));
QPsth = 4*Gfr^2*M*yps.^2.*(yps.^2-s.^2.*rhoE.*etae./(M*varsigma)).*dpspf./(s.*(H*M-C^2));
% h_{11}-source 
tau1111 = 2*(KPfth.*(GPf11-GPf)-KPsth.*(GPs11-GPs)) ...
    + QPfth.*(GPf1111-2*GPf11+GPf)-QPsth.*(GPs1111-2*GPs11+GPs) ...
    + NPfth.*GPf-NPsth.*GPs + 4*(KSth.*ys.^2.*(GS1111-GS11)-KEMth.*yem.^2.*(GEM1111-GEM11));
tau2111 = KPfth.*GPf21-KPsth.*GPs21 + QPfth.*(GPf2111-GPf21)-QPsth.*(GPs2111-GPs21) ...
    + 2*(KSth.*ys.^2.*(2.*GS2111-GS21)-KEMth.*yem.^2.*(2.*GEM2111-GEM21));
tau3111 = KPfth.*GPf31-KPsth.*GPs31 + QPfth.*(GPf3111-GPf31)-QPsth.*(GPs3111-GPs31) ...
    + 2*(KSth.*ys.^2.*(2.*GS3111-GS31)-KEMth.*yem.^2.*(2*GEM3111-GEM31));
tau2211 = KPfth.*(GPf11+GPf22-2*GPf)-KPsth.*(GPs11+GPs22-2*GPs) ...
    + QPfth.*(GPf2211-(GPf11+GPf22)+GPf)-QPsth.*(GPs2211-(GPs11+GPs22)+GPs) ...
    + NPfth.*GPf-NPsth.*GPs + 4*(KSth.*ys.^2.*GS2211-KEMth.*yem.^2.*GEM2211);
tau3211 = KPfth.*GPf32-KPsth.*GPs32 + QPfth.*(GPf3211-GPf32)-QPsth.*(GPs3211-GPs32) ...
    + 4*(KSth.*ys.^2.*(GS3211-GS32)-KEMth.*yem.^2.*(GEM3211-GEM32));
tau3311 = KPfth.*(GPf11+GPf33-2*GPf)-KPsth.*(GPs11+GPs33-2*GPs) ...
    + QPfth.*(GPf3311-(GPf11+GPf33)+GPf)-QPsth.*(GPs3311-(GPs11+GPs33)+GPs) ...
    + NPfth.*GPf-NPsth.*GPs + 4*(KSth.*ys.^2.*GS3311-KEMth.*yem.^2.*GEM3311);
% h_{21}-source 
tau1121 = tau2111;
tau2121 = QPfth.*GPf2211-QPsth.*GPs2211 + KSth.*ys.^2.*(4*GS2211-GS11-GS22)...
    -KEMth.*yem.^2.*(4*GEM2211-GEM11-GEM22);
tau3121 = QPfth.*GPf3211-QPsth.*GPs3211 + KSth.*ys.^2.*(4*GS3211-GS32)...
    -KEMth.*yem.^2.*(4*GEM3211-GEM32);
tau2221 = KPfth.*GPf21-KPsth.*GPs21 + QPfth.*(GPf2221-GPf21)-QPsth.*(GPs2221-GPs21)...
    + 2*KSth.*ys.^2.*(2*GS2221-GS21)-2*KEMth.*yem.^2.*(2*GEM2221-GEM21);
tau3221 = QPfth.*GPf3221-QPsth.*GPs3221 + KSth.*ys.^2.*(4*GS3221-GS31)...
    -KEMth.*yem.^2.*(4*GEM3221-GEM31);
tau3321 = KPfth.*GPf21-KPsth.*GPs21 + QPfth.*(GPf3321-GPf21)-QPsth.*(GPs3321-GPs21)...
    + 4*KSth.*ys.^2.*GS3321-4*KEMth.*yem.^2.*GEM3321;
% h_{31}-source 
tau1131 = tau3111;
tau2131 = tau3121;
tau3131 = QPfth.*GPf3311-QPsth.*GPs3311 + KSth.*ys.^2.*(4*GS3311-GS11-GS33)...
    -KEMth.*yem.^2.*(4*GEM3311-GEM11-GEM33);
tau2231 = KPfth.*GPf31-KPsth.*GPs31 + QPfth.*(GPf3221-GPf31)-QPsth.*(GPs3221-GPs31)...
    + 4*KSth.*ys.^2.*GS3221-4*KEMth.*yem.^2.*GEM3221;
tau3231 = QPfth.*GPf3321-QPsth.*GPs3321 + KSth.*ys.^2.*(4*GS3321-GS21)...
    -KEMth.*yem.^2.*(4*GEM3321-GEM21);
tau3331 = KPfth.*GPf31-KPsth.*GPs31 + QPfth.*(GPf3331-GPf31)-QPsth.*(GPs3331-GPs31)...
    + 2*KSth.*ys.^2.*(2*GS3331-GS31)-2*KEMth.*yem.^2.*(2*GEM3331-GEM31);
% h_{22}-source 
tau1122 = tau2211;
tau2122 = tau2221;
tau3122 = tau2231;
tau2222 = 2*(KPfth.*(GPf11-GPf)-KPsth.*(GPs11-GPs)) ...
    + QPfth.*(GPf2222-2*GPf22+GPf)-QPsth.*(GPs2222-2*GPs22+GPs) ...
    + NPfth.*GPf-NPsth.*GPs + 4*(KSth.*ys.^2.*(GS2222-GS22)-KEMth.*yem.^2.*(GEM2222-GEM22));
tau3222 = KPfth.*GPf32-KPsth.*GPs32 + QPfth.*(GPf3222-GPf32)-QPsth.*(GPs3222-GPs32)...
    + 2*KSth.*ys.^2.*(2*GS3222-GS32)-2*KEMth.*yem.^2.*(2*GEM3222-GEM32);
tau3322 = KPfth.*(GPf22+GPf33-2*GPf)-KPsth.*(GPs22+GPs33-2*GPs) ...
    + QPfth.*(GPf3322-(GPf22+GPf33)+GPf)-QPsth.*(GPs3322-(GPs22+GPs33)+GPs) ...
    + NPfth.*GPf-NPsth.*GPs + 4*(KSth.*ys.^2.*GS3322-KEMth.*yem.^2.*GEM3322);
% h_{32}-source 
tau1132 = tau3211;
tau2132 = tau3221;
tau3132 = tau3231;
tau2232 = tau3222;
tau3232 = QPfth.*GPf3322-QPsth.*GPs3322 + KSth.*ys.^2.*(4*GS3322-GS22-GS33)...
    -KEMth.*yem.^2.*(4*GEM3322-GEM22-GEM33);
tau3332 = KPfth.*GPf32-KPsth.*GPs32 + QPfth.*(GPf3332-GPf32)-QPsth.*(GPs3332-GPs32)...
    + 2*KSth.*ys.^2.*(2*GS3332-GS32)-2*KEMth.*yem.^2.*(2*GEM3332-GEM32);
% h_{33}-source 
tau1133 = tau3311;
tau2133 = tau3321;
tau3133 = tau3331;
tau2233 = tau3322;
tau3233 = tau3332;
tau3333 = 2*(KPfth.*(GPf33-GPf)-KPsth.*(GPs33-GPs)) ...
    + QPfth.*(GPf3333-2*GPf33+GPf)-QPsth.*(GPs3333-2*GPs33+GPs) ...
    + NPfth.*GPf-NPsth.*GPs + 4*(KSth.*ys.^2.*(GS3333-GS33)-KEMth.*yem.^2.*(GEM3333-GEM33));
% compute and plot some time-fields
tv111=2*real(ifft(v111.*wav,mf*nf,3))/dt;
tE211=2*real(ifft(E211.*wav,mf*nf,3))/dt;
tw111=2*real(ifft(w111.*wav,mf*nf,3))/dt;
tp11=2*real(ifft(p11.*wav,mf*nf,3))/dt;
tH211=2*real(ifft(H211.*wav,mf*nf,3))/dt;
ttau1111=2*real(ifft(tau1111.*wav,mf*nf,3))/dt;
for it=1:nf
    imagesc(xx1,xx2,squeeze(tE211(:,:,it)));
    set(gca,'fontsize',fs)
    colormap(1-gray)
    colorbar;
    set(gca,'fontsize',18,'plotboxaspectratio',[1 1 1])
    xlabel('x distance from source [m]')
    ylabel('y distance from source [m]')
    pause(0.1)
end
