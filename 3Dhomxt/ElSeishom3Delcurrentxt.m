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
varsigma        =   etae - s.*rhof.*L.^2;
chi             =   s.*rhof*L;
% source signature
wav=-2*sqrt(1/pi)*(s/(2*pi*fc)).^2.*exp((s/(2*pi*fc)).^2)/fc;
% spherical wavenumbers of eqs B-9, B-10, B-12, B-13.
ypf     =   s./sqrt(2).*sqrt((rhoB.*M-2.*rhof.*C+rhof.*H.*etae./varsigma)./(H.*M-C.^2) - sqrt(((rhoB.*M-2.*rhof.*C+rhof.*H.*etae./varsigma)./(H.*M-C.^2)).^2 + 4.*(rhof.^2-rhoB.*rhof.*etae./varsigma)./(H.*M-C.^2)));
yps     =   s./sqrt(2).*sqrt((rhoB.*M-2.*rhof.*C+rhof.*H.*etae./varsigma)./(H.*M-C.^2) + sqrt(((rhoB.*M-2.*rhof.*C+rhof.*H.*etae./varsigma)./(H.*M-C.^2)).^2 + 4.*(rhof.^2-rhoB.*rhof.*etae./varsigma)./(H.*M-C.^2)));
ys      =   1/sqrt(2).*sqrt(s.^2.*rhoc./Gfr+zeta.*etae - sqrt((s.^2.*rhoc./Gfr-zeta.*etae).^2-4.*s.^3.*zeta.*(rhof.*L).^2./Gfr));
yem     =   1/sqrt(2).*sqrt(s.^2.*rhoc./Gfr+zeta.*etae + sqrt((s.^2.*rhoc./Gfr-zeta.*etae).^2-4.*s.^3.*zeta.*(rhof.*L).^2./Gfr));
% recurring coefficient factors
dems    =   1./(yem.^2-ys.^2);            
dpspf   =   1./(yps.^2-ypf.^2);     

% electric current source

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
% particle velocity
KSSve = -s.*chi.*zeta.*dems/Gfr;
KPsve = s.^2.*rhof.*L*C.*(s.^2*rhof./C-yps.^2).*dpspf./((H*M-C^2)*varsigma);
KPfve = s.^2.*rhof.*L*C.*(s.^2*rhof./C-ypf.^2).*dpspf./((H*M-C^2)*varsigma);
v11 = KSSve.*(GS11-GS-(GEM11-GEM)) + KPfve.*GPf11-KPsve.*GPs11;
v12 = KSSve.*(GS21-GEM31) + KPfve.*GPf21-KPsve.*GPs21;
v13 = KSSve.*(GS31-GEM31) + KPfve.*GPf31-KPsve.*GPs31;
v22 = KSSve.*(GS22-GS-(GEM22-GEM)) + KPfve.*GPf22-KPsve.*GPs22;
v32 = KSSve.*(GS32-GEM32) + KPfve.*GPf32-KPsve.*GPs32;
v33 = KSSve.*(GS33-GS-(GEM33-GEM)) + KPfve.*GPf33-KPsve.*GPs33;
% electric field
KSee = zeta.*(s.^2.*rhoc/Gfr-ys.^2).*dems;
KEMee = zeta.*(s.^2.*rhoc/Gfr-yem.^2).*dems;
KPsee = (s.*rhof.*L./varsigma).^2.*s.*H.*(s.^2.*rhoB./H-yps.^2).*dpspf./(H*M-C^2);
KPfee = (s.*rhof.*L./varsigma).^2.*s.*H.*(s.^2.*rhoB./H-ypf.^2).*dpspf./(H*M-C^2);
E11 = KSee.*(GS11-GS)-KEMee.*(GEM11-GEM) + KPfee.*GPf11-KPsee.*GPs11;
E12 = KSee.*GS21-KEMee.*GEM31 + KPfee.*GPf21-KPsee.*GPs21;
E13 = KSee.*GS31-KEMee.*GEM31 + KPfee.*GPf31-KPsee.*GPs31;
E22 = KSee.*(GS22-GS)-KEMee.*(GEM22-GEM) + KPfee.*GPf22-KPsee.*GPs22;
E32 = KSee.*GS32-KEMee.*GEM32 + KPfee.*GPf32-KPsee.*GPs32;
E33 = KSee.*(GS33-GS)-KEMee.*(GEM33-GEM) + KPfee.*GPf33-KPsee.*GPs33;
% filtration velocity
KSwe = zeta.*L.*(s.^2.*rhoB/Gfr-ys.^2).*dems;
KEMwe = zeta.*L.*(s.^2.*rhoB/Gfr-yem.^2).*dems;
KPswe = -s.^2.*rhof.*L.*H.*(s.^2.*rhoB/H-yps.^2).*dpspf./((H*M-C^2).*varsigma);
KPfwe = -s.^2.*rhof.*L.*H.*(s.^2.*rhoB/H-ypf.^2).*dpspf./((H*M-C^2).*varsigma);
w11 = KSwe.*(GS11-GS)-KEMwe.*(GEM11-GEM) + KPfwe.*GPf11-KPswe.*GPs11;
w12 = KSwe.*GS21-KEMwe.*GEM31 + KPfwe.*GPf21-KPswe.*GPs21;
w13 = KSwe.*GS31-KEMwe.*GEM31 + KPfwe.*GPf31-KPswe.*GPs31;
w22 = KSwe.*(GS22-GS)-KEMwe.*(GEM22-GEM) + KPfwe.*GPf22-KPswe.*GPs22;
w32 = KSwe.*GS32-KEMwe.*GEM32 + KPfwe.*GPf32-KPswe.*GPs32;
w33 = KSwe.*(GS33-GS)-KEMwe.*(GEM33-GEM) + KPfwe.*GPf33-KPswe.*GPs33;
% acoustic pressure
KPspe = s.*rhof.*L.*(s.^2.*(rhoB*M-rhof*C)/(H*M-C^2)-yps.^2).*dpspf./varsigma;
KPfpe = s.*rhof.*L.*(s.^2.*(rhoB*M-rhof*C)/(H*M-C^2)-ypf.^2).*dpspf./varsigma;
p1 = KPfpe.*GPf1-KPspe.*GPs1;
p2 = KPfpe.*GPf2-KPspe.*GPs2;
p3 = KPfpe.*GPf3-KPspe.*GPs3;
% magnetic field
KSme = (s.^2.*rhoc./Gfr-ys.^2).*dems;
KEMme = (s.^2.*rhoc./Gfr-yem.^2).*dems;
% H11 = 0;
H12 = -KSme.*GS3+KEMme.*GEM3;
H13 = KSme.*GS2-KEMme.*GEM2;
% H21 = -H12;
% H22 = 0;
H23 = -KSme.*GS1+KEMme.*GEM1;
% H31 = -H13;
% H32 = -H23;
% H33 = 0;
% bulk stress field
KSte = -chi.*zeta.*dems;
KEMte = KSte;
KPste = 2*s.*rhof.*L.*C.*Gfr.*(s.^2.*rhof/C-yps.^2).*dpspf./((H*M-C^2).*varsigma);
KPfte = 2*s.*rhof.*L.*C.*Gfr.*(s.^2.*rhof/C-ypf.^2).*dpspf./((H*M-C^2).*varsigma);
NPste = s.^3.*rhof.*L.*(rhof*H-rhoB*C).*dpspf./((H*M-C^2).*varsigma);
NPfte = NPste;
tau111 = 2*(KSte.*(GS111-GS1)-KEMte.*(GEM111-GEM1))+KPfte.*(GPf111-GPf1)...
    -KPste.*(GPs111-GPs1) + NPfte.*GPf1-NPste.*GPs1;
tau211 = KSte.*(2*GS211-GS2)-KEMte.*(2*GEM211-GEM2)+KPfte.*GPf211-KPste.*GPs211;
tau311 = KSte.*(2*GS311-GS3)-KEMte.*(2*GEM311-GEM3)+KPfte.*GPf311-KPste.*GPs311;
tau221 = 2*KSte.*GS221-2*KEMte.*GEM221+KPfte.*(GPf221-GPf1)-KPste.*(GPs221-GPs1) ...
      + NPfte.*GPf1-NPste.*GPs1;
tau321 = 2*KSte.*GS321-2*KEMte.*GEM321+KPfte.*GPf321-KPste.*GPs321;
tau331 = 2*KSte.*GS331-2*KEMte.*GEM331+KPfte.*(GPf331-GPf1)-KPste.*(GPs331-GPs1) ...
      + NPfte.*GPf1-NPste.*GPs1;
tau112 = 2*KSte.*GS211-2*KEMte.*GEM211+KPfte.*(GPf211-GPf2)-KPste.*(GPs211-GPs2)...
    + NPfte.*GPf2-NPste.*GPs2;
tau212 = KSte.*(2*GS221-GS1)-KEMte.*(2*GEM221-GEM1) + KPfte.*GPf221-KPste.*GPs211;
% tau312 = tau321;
tau222 = 2*(KSte.*(GS222-GS2)-KEMte.*(GEM222-GEM2))+KPfte.*(GPf222-GPf2)...
    -KPste.*(GPs222-GPs2) + NPfte.*GPf2-NPste.*GPs2;
tau322 = KSte.*(2*GS322-GS3)-KEMte.*(2*GEM322-GEM3)+KPfte.*GPf322-KPste.*GPs322;
tau332 = 2*KSte.*GS332-2*KEMte.*GEM322+KPfte.*(GPf332-GPf2)-KPste.*(GPs332-GPs2)...
    + NPfte.*GPf2-NPste.*GPs2;
tau113 = 2*KSte.*GS311-2*KEMte.*GEM311+KPfte.*(GPf311-GPf3)-KPste.*(GPs311-GPs3)...
    + NPfte.*GPf3-NPste.*GPs3;
% tau213 = tau321;
tau313 = KSte.*(2*GS331-GS1)-KEMte.*(2*GEM331-GEM1) + KPfte.*GPf331-KPste.*GPs331;
tau223 = 2*KSte.*GS322-2*KEMte.*GEM322+KPfte.*(GPf322-GPf3)-KPste.*(GPs322-GPs3)...
    + NPfte.*GPf3-NPste.*GPs3;
tau323 = KSte.*(2*GS332-GS2)-KEMte.*(2*GEM332-GEM2) + KPfte.*GPf332-KPste.*GPs332;
tau333 = 2*(KSte.*(GS333-GS3)-KEMte.*(GEM333-GEM3))+KPfte.*(GPf333-GPf3)...
    -KPste.*(GPs333-GPs3) + NPfte.*GPf3-NPste.*GPs3;
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
