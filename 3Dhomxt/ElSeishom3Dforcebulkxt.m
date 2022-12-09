% SPDX-License-Identifier: CC0-1.0

% defining all the variables
clear all;
clc;
close all;
% predetermine fontsize
fs=14;
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

% force acting on the bulk

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
% for v, e, and w we need only the upper triangle of the 3x3 matrix
% particle velocity
KSvf = -s.*(zeta.*etae-ys.^2).*dems/Gfr;
KEMvf = -s.*(zeta.*etae-yem.^2).*dems/Gfr;
KPsvf =  s.*M.*(s.^2.*rhoE.*etae./(M.*varsigma)-yps.^2).*dpspf/(H*M-C^2);
KPfvf =  s.*M.*(s.^2.*rhoE.*etae./(M.*varsigma)-ypf.^2).*dpspf/(H*M-C^2);
v11 = KSvf.*(GS11-GS)-KEMvf.*(GEM11-GEM) + KPfvf.*GPf11-KPsvf.*GPs11;
v12 = KSvf.*GS21-KEMvf.*GEM31 + KPfvf.*GPf21-KPsvf.*GPs21;
v13 = KSvf.*GS31-KEMvf.*GEM31 + KPfvf.*GPf31-KPsvf.*GPs31;
v22 = KSvf.*(GS22-GS)-KEMvf.*(GEM22-GEM) + KPfvf.*GPf22-KPsvf.*GPs22;
v23 = KSvf.*GS32-KEMvf.*GEM32 + KPfvf.*GPf32-KPsvf.*GPs32;
v33 = KSvf.*(GS33-GS)-KEMvf.*(GEM33-GEM) + KPfvf.*GPf33-KPsvf.*GPs33;
% electric field
KSSef = -s.*chi.*zeta.*dems/Gfr;
KPsef =  s.^2.*rhoE.*L*C.*(s.^2*rhof./C-yps.^2).*dpspf./((H*M-C^2)*varsigma);
KPfef =  s.^2.*rhoE.*L*C.*(s.^2*rhof./C-ypf.^2).*dpspf./((H*M-C^2)*varsigma);
E11 = KSSef.*(GS11-GS-(GEM11-GEM)) + KPfef.*GPf11-KPsef.*GPs11;
E12 = KSSef.*(GS21-GEM31) + KPfef.*GPf21-KPsef.*GPs21;
E13 = KSSef.*(GS31-GEM31) + KPfef.*GPf31-KPsef.*GPs31;
E22 = KSSef.*(GS22-GS-(GEM22-GEM)) + KPfef.*GPf22-KPsef.*GPs22;
E23 = KSSef.*(GS32-GEM32) + KPfef.*GPf32-KPsef.*GPs32;
E33 = KSSef.*(GS33-GS-(GEM33-GEM)) + KPfef.*GPf33-KPsef.*GPs33;
% filtration velocity
KSwf =  s.*rhof.*(zeta.*varsigma-ys.^2).*dems./(rhoE*Gfr);
KEMwf =  s.*rhof.*(zeta.*varsigma-yem.^2).*dems./(rhoE*Gfr);
KPswf = -s*C.*(s.^2*rhof./C-yps.^2).*dpspf./(H*M-C^2);
KPfwf = -s*C.*(s.^2*rhof./C-ypf.^2).*dpspf./(H*M-C^2);
w11 = KSwf.*(GS11-GS)-KEMwf.*(GEM11-GEM) + KPfwf.*GPf11-KPswf.*GPs11;
w12 = KSwf.*GS21-KEMwf.*GEM31 + KPfwf.*GPf21-KPswf.*GPs21;
w13 = KSwf.*GS31-KEMwf.*GEM31 + KPfwf.*GPf31-KPswf.*GPs31;
w22 = KSwf.*(GS22-GS)-KEMwf.*(GEM22-GEM) + KPfwf.*GPf22-KPswf.*GPs22;
w23 = KSwf.*GS32-KEMwf.*GEM32 + KPfwf.*GPf32-KPswf.*GPs32;
w33 = KSwf.*(GS33-GS)-KEMwf.*(GEM33-GEM) + KPfwf.*GPf33-KPswf.*GPs33;
% acoustic pressure
KPPpf = -s.^2.*(rhoE.*etae*C./varsigma-rhof*M).*dpspf/(H*M-C^2);
p1 = KPPpf.*(GPf1-GPs1);
p2 = KPPpf.*(GPf2-GPs2);
p3 = KPPpf.*(GPf3-GPs3);
% magnetic field
KSSmf = -s.*chi.*dems/Gfr;
% H11 = 0;
H12 = -KSSmf.*(GS3-GEM3);
H13 = KSSmf.*(GS2-GEM2);
% H21 = -H12;
% H22 = 0;
H23 = -KSSmf.*(GS1-GEM1);
% H31 = -H13;
% H32 = -H23;
% H33 = 0;
% bulk stress field; tau_{ij}=tau_{ji}
KStf = -(zeta.*etae-ys.^2).*dems;
KEMtf = -(zeta.*etae-yem.^2).*dems;
KPstf = 2*Gfr*M*(s.^2.*rhoE.*etae./(M*varsigma)-yps.^2).*dpspf/(H*M-C^2);
KPftf = 2*Gfr*M*(s.^2.*rhoE.*etae./(M*varsigma)-ypf.^2).*dpspf/(H*M-C^2);
NPstf = (s.^2.*(rhoE.*H.*etae./varsigma-rhof*C)/(H*M-C^2)-yps.^2).*dpspf;
NPftf = (s.^2.*(rhoE.*H.*etae./varsigma-rhof*C)/(H*M-C^2)-ypf.^2).*dpspf;
tau111 = 2*(KStf.*(GS111-GS1)-KEMtf.*(GEM111-GEM1))+KPftf.*(GPf111-GPf1)...
    -KPstf.*(GPs111-GPs1) + NPftf.*GPf1-NPstf.*GPs1;
tau211 = KStf.*(2*GS211-GS2)-KEMtf.*(2*GEM211-GEM2)+KPftf.*GPf211-KPstf.*GPs211;
tau311 = KStf.*(2*GS311-GS3)-KEMtf.*(2*GEM311-GEM3)+KPftf.*GPf311-KPstf.*GPs311;
tau221 = 2*KStf.*GS221-2*KEMtf.*GEM221+KPftf.*(GPf221-GPf1)-KPstf.*(GPs221-GPs1) ...
      + NPftf.*GPf1-NPstf.*GPs1;
tau321 = 2*KStf.*GS321-2*KEMtf.*GEM321+KPftf.*GPf321-KPstf.*GPs321;
tau331 = 2*KStf.*GS331-2*KEMtf.*GEM331+KPftf.*(GPf331-GPf1)-KPstf.*(GPs331-GPs1) ...
      + NPftf.*GPf1-NPstf.*GPs1;
tau112 = 2*KStf.*GS211-2*KEMtf.*GEM211+KPftf.*(GPf211-GPf2)-KPstf.*(GPs211-GPs2)...
    + NPftf.*GPf2-NPstf.*GPs2;
tau212 = KStf.*(2*GS221-GS1)-KEMtf.*(2*GEM221-GEM1) + KPftf.*GPf221-KPstf.*GPs221;
% tau312 = tau321;
tau222 = 2*(KStf.*(GS222-GS2)-KEMtf.*(GEM222-GEM2))+KPftf.*(GPf222-GPf2)...
    -KPstf.*(GPs222-GPs2) + NPftf.*GPf2-NPstf.*GPs2;
tau322 = KStf.*(2*GS322-GS3)-KEMtf.*(2*GEM322-GEM3)+KPftf.*GPf322-KPstf.*GPs322;
tau332 = 2*KStf.*GS332-2*KEMtf.*GEM322+KPftf.*(GPf332-GPf2)-KPstf.*(GPs332-GPs2)...
    + NPftf.*GPf-NPstf.*GPs2;
tau113 = 2*KStf.*GS311-2*KEMtf.*GEM311+KPftf.*(GPf311-GPf3)-KPstf.*(GPs311-GPs3)...
    + NPftf.*GPf3-NPstf.*GPs3;
% tau213 = tau321;
tau313 = KStf.*(2*GS331-GS1)-KEMtf.*(2*GEM331-GEM1) + KPftf.*GPf331-KPstf.*GPs331;
tau223 = 2*KStf.*GS322-2*KEMtf.*GEM322+KPftf.*(GPf322-GPf3)-KPstf.*(GPs322-GPs3)...
    + NPftf.*GPf3-NPstf.*GPs3;
tau323 = KStf.*(2*GS332-GS2)-KEMtf.*(2*GEM332-GEM2) + KPftf.*GPf332-KPstf.*GPs332;
tau333 = 2*(KStf.*(GS333-GS3)-KEMtf.*(GEM333-GEM3))+KPftf.*(GPf333-GPf3)...
    -KPstf.*(GPs333-GPs3) + NPftf.*GPf3-NPstf.*GPs3;
% compute and plot some time-fields
tv11=2*real(ifft(v11.*wav,mf*nf,3))/dt;
tE11=2*real(ifft(E11.*wav,mf*nf,3))/dt;
tw11=2*real(ifft(w11.*wav,mf*nf,3))/dt;
tp1=2*real(ifft(p1.*wav,mf*nf,3))/dt;
tH12=2*real(ifft(H12.*wav,mf*nf,3))/dt;
ttau111=2*real(ifft(tau111.*wav,mf*nf,3))/dt;
for it=1:nf
    imagesc(xx1,xx2,squeeze(ttau111(:,:,it)));
    colormap(1-gray)
    colorbar;
    set(gca,'fontsize',18,'plotboxaspectratio',[1 1 1])
    xlabel('x distance from source [m]')
    ylabel('y distance from source [m]')
    pause(0.1)
end
