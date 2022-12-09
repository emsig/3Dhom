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

% force acting on the fluid

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
[GPs111,GPs211,GPs311,GPs221,GPs321,GPs331,GPs222,GPs322,GPs332,GPs333]=Greenijk(ypf,th1,th2,th3,R);
[GPf111,GPf211,GPf311,GPf221,GPf321,GPf331,GPf222,GPf322,GPf332,GPf333]=Greenijk(yps,th1,th2,th3,R);

% particle velocity
KSvff = s.*rhof.*(zeta.*varsigma-ys.^2).*dems./(rhoE*Gfr);
KEMvff = s.*rhof.*(zeta.*varsigma-yem.^2).*dems./(rhoE*Gfr);
KPsvff = -s*C.*(s.^2*rhof./C-yps.^2).*dpspf./(H*M-C^2);
KPfvff = -s*C.*(s.^2*rhof./C-ypf.^2).*dpspf./(H*M-C^2);
v11 = KSvff.*(GS11-GS)-KEMvff.*(GEM11-GEM) + KPfvff.*GPf11-KPsvff.*GPs11;
v12 = KSvff.*GS21-KEMvff.*GEM31 + KPfvff.*GPf21-KPsvff.*GPs21;
v13 = KSvff.*GS31-KEMvff.*GEM31 + KPfvff.*GPf31-KPsvff.*GPs31;
v22 = KSvff.*(GS22-GS)-KEMvff.*(GEM22-GEM) + KPfvff.*GPf22-KPsvff.*GPs22;
v32 = KSvff.*GS32-KEMvff.*GEM32 + KPfvff.*GPf32-KPsvff.*GPs32;
v33 = KSvff.*(GS33-GS)-KEMvff.*(GEM33-GEM) + KPfvff.*GPf33-KPsvff.*GPs33;
% electric field
KSeff = zeta.*L.*(s.^2.*rhoB/Gfr-ys.^2).*dems;
KEMeff = zeta.*L.*(s.^2.*rhoB/Gfr-yem.^2).*dems;
KPseff = -s.^2.*rhoE.*L.*H.*(s.^2.*rhoB/H-yps.^2).*dpspf./((H*M-C^2).*varsigma);
KPfeff = -s.^2.*rhoE.*L.*H.*(s.^2.*rhoB/H-ypf.^2).*dpspf./((H*M-C^2).*varsigma);
E11 = KSeff.*(GS11-GS)-KEMeff.*(GEM11-GEM) + KPfeff.*GPf11-KPseff.*GPs11;
E12 = KSeff.*GS21-KEMeff.*GEM31 + KPfeff.*GPf21-KPseff.*GPs21;
E13 = KSeff.*GS31-KEMeff.*GEM31 + KPfeff.*GPf31-KPseff.*GPs31;
E22 = KSeff.*(GS22-GS)-KEMeff.*(GEM22-GEM) + KPfeff.*GPf22-KPseff.*GPs22;
E32 = KSeff.*GS32-KEMeff.*GEM32 + KPfeff.*GPf32-KPseff.*GPs32;
E33 = KSeff.*(GS33-GS)-KEMeff.*(GEM33-GEM) + KPfeff.*GPf33-KPseff.*GPs33;
% filtration velocity
KSwff = -(s.^2.*rhoB/Gfr-ys.^2).*(zeta.*varsigma-ys.^2).*dems./(s.*rhoE);
KEMwff = -(s.^2.*rhoB/Gfr-yem.^2).*(zeta.*varsigma-yem.^2).*dems./(s.*rhoE);
KPswff = s.*H.*(s.^2.*rhoB/H-yps.^2).*dpspf./(H*M-C^2);
KPfwff = s.*H.*(s.^2.*rhoB/H-ypf.^2).*dpspf./(H*M-C^2);
w11 = KSwff.*(GS11-GS)-KEMwff.*(GEM11-GEM) + KPfwff.*GPf11-KPswff.*GPs11;
w12 = KSwff.*GS21-KEMwff.*GEM31 + KPfwff.*GPf21-KPswff.*GPs21;
w13 = KSwff.*GS31-KEMwff.*GEM31 + KPfwff.*GPf31-KPswff.*GPs31;
w22 = KSwff.*(GS22-GS)-KEMwff.*(GEM22-GEM) + KPfwff.*GPf22-KPswff.*GPs22;
w32 = KSwff.*GS32-KEMwff.*GEM32 + KPfwff.*GPf32-KPswff.*GPs32;
w33 = KSwff.*(GS33-GS)-KEMwff.*(GEM33-GEM) + KPfwff.*GPf33-KPswff.*GPs33;
% acoustic pressure
KPspff = -(s.^2*(rhoB.*M-rhof*C)./(H*M-C^2)-yps.^2).*dpspf;
KPfpff = -(s.^2*(rhoB.*M-rhof*C)./(H*M-C^2)-ypf.^2).*dpspf;
p1 = KPfpff.*GPf1-KPspff.*GPs1;
p2 = KPfpff.*GPf2-KPspff.*GPs2;
p3 = KPfpff.*GPf3-KPspff.*GPs3;
% magnetic field
KSmff = -L.*(s.^2.*rhoB/Gfr-ys.^2).*dems;
KEMmff = -L.*(s.^2.*rhoB/Gfr-yem.^2).*dems;
% H11 = 0;
H12 = -KSmff.*GS3+KEMmff.*GEM3;
H13 = KSmff.*GS2-KEMmff.*GEM2;
% H21 = -H12;
% H22 = 0;
H23 = -KSmff.*GS1+KEMmff.*GEM1;
% H31 = -H13;
% H32 = -H23;
% H33 = 0;
% bulk stress field
KStff = rhof.*(zeta.*varsigma-ys.^2).*dems./rhoE;
KEMtff = rhof.*(zeta.*varsigma-yem.^2).*dems./rhoE;
KPstff = -2*Gfr*C.*(s.^2*rhof/C-yps.^2).*dpspf./(H*M-C^2);
KPftff = -2*Gfr*C.*(s.^2*rhof/C-ypf.^2).*dpspf./(H*M-C^2);
NPftff = -s.^2*(rhof*H-rhoB*C).*dpspf./(H*M-C^2);
NPstff = NPftff;
tau111 = 2*(KStff.*(GS111-GS1)-KEMtff.*(GEM111-GEM1))+KPftff.*(GPf111-GPf1)...
    -KPstff.*(GPs111-GPs1) + NPftff.*GPf1-NPstff.*GPs1;
tau211 = KStff.*(2*GS211-GS2)-KEMtff.*(2*GEM211-GEM2)+KPftff.*GPf211-KPstff.*GPs211;
tau311 = KStff.*(2*GS311-GS3)-KEMtff.*(2*GEM311-GEM3)+KPftff.*GPf311-KPstff.*GPs311;
tau221 = 2*KStff.*GS221-2*KEMtff.*GEM221+KPftff.*(GPf221-GPf1)-KPstff.*(GPs221-GPs1) ...
      + NPftff.*GPf1-NPstff.*GPs1;
tau321 = 2*KStff.*GS321-2*KEMtff.*GEM321+KPftff.*GPf321-KPstff.*GPs321;
tau331 = 2*KStff.*GS331-2*KEMtff.*GEM331+KPftff.*(GPf331-GPf1)-KPstff.*(GPs331-GPs1) ...
      + NPftff.*GPf1-NPstff.*GPs1;
tau112 = 2*KStff.*GS211-2*KEMtff.*GEM211+KPftff.*(GPf211-GPf2)-KPstff.*(GPs211-GPs2)...
    + NPftff.*GPf2-NPstff.*GPs2;
tau212 = KStff.*(2*GS221-GS1)-KEMtff.*(2*GEM221-GEM1) + KPftff.*GPf221-KPstff.*GPs211;
% tau312 = tau321;
tau222 = 2*(KStff.*(GS222-GS2)-KEMtff.*(GEM222-GEM2))+KPftff.*(GPf222-GPf2)...
    -KPstff.*(GPs222-GPs2) + NPftff.*GPf2-NPstff.*GPs2;
tau322 = KStff.*(2*GS322-GS3)-KEMtff.*(2*GEM322-GEM3)+KPftff.*GPf322-KPstff.*GPs322;
tau332 = 2*KStff.*GS332-2*KEMtff.*GEM322+KPftff.*(GPf332-GPf2)-KPstff.*(GPs332-GPs2)...
    + NPftff.*GPf2-NPstff.*GPs2;
tau113 = 2*KStff.*GS311-2*KEMtff.*GEM311+KPftff.*(GPf311-GPf3)-KPstff.*(GPs311-GPs3)...
    + NPftff.*GPf3-NPstff.*GPs3;
% tau213 = tau321;
tau313 = KStff.*(2*GS331-GS1)-KEMtff.*(2*GEM331-GEM1) + KPftff.*GPf331-KPstff.*GPs331;
tau223 = 2*KStff.*GS322-2*KEMtff.*GEM322+KPftff.*(GPf322-GPf3)-KPstff.*(GPs322-GPs3)...
    + NPftff.*GPf3-NPstff.*GPs3;
tau323 = KStff.*(2*GS332-GS2)-KEMtff.*(2*GEM332-GEM2) + KPftff.*GPf332-KPstff.*GPs332;
tau333 = 2*(KStff.*(GS333-GS3)-KEMtff.*(GEM333-GEM3))+KPftff.*(GPf333-GPf3)...
    -KPstff.*(GPs333-GPs3) + NPftff.*GPf3-NPstff.*GPs3;
% compute and plot some time-fields
tv11=2*real(ifft(v11.*wav,mf*nf,3))/dt;
tE11=2*real(ifft(E11.*wav,mf*nf,3))/dt;
tw11=2*real(ifft(w11.*wav,mf*nf,3))/dt;
tp1=2*real(ifft(p1.*wav,mf*nf,3))/dt;
tH12=2*real(ifft(H12.*wav,mf*nf,3))/dt;
ttau111=2*real(ifft(tau111.*wav,mf*nf,3))/dt;
for it=1:nf
    imagesc(xx1,xx2,squeeze(tH12(:,:,it)));
    colormap(1-gray)
    colorbar;
    set(gca,'fontsize',18,'plotboxaspectratio',[1 1 1])
    xlabel('x distance from source [m]')
    ylabel('y distance from source [m]')
    pause(0.1)
end
