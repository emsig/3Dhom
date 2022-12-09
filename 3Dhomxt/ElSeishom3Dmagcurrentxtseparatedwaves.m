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
% spherical wavenumbers of eqs B-9 and B-10.
ys      =   1/sqrt(2).*sqrt(s.^2.*rhoc./Gfr+zeta.*etae - sqrt((s.^2.*rhoc./Gfr-zeta.*etae).^2-4.*s.^3.*zeta.*(rhof.*L).^2./Gfr));
yem     =   1/sqrt(2).*sqrt(s.^2.*rhoc./Gfr+zeta.*etae + sqrt((s.^2.*rhoc./Gfr-zeta.*etae).^2-4.*s.^3.*zeta.*(rhof.*L).^2./Gfr));
% recurring coefficient factors
dems    =   1./(yem.^2-ys.^2);            

% magnetic current, J_1^m only

% Green's functions of eq C-10
GS = exp(-ys.*R)./(4*pi*R);
GEM = exp(-yem.*R)./(4*pi*R);
% grad (equation C-15)
[GS1,GS2,GS3]=Greeni(ys,th1,th2,th3,R);
[GEM1,GEM2,GEM3]=Greeni(yem,th1,th2,th3,R);
% grad-grad divided by gamma^2 (equation C-16)
[GS11,GS21,GS31]=Greenij(ys,th1,th2,th3,R);
[GEM11,GEM21,GEM31]=Greenij(yem,th1,th2,th3,R);
% for stress field we need GS22, GS32, and GS33
GS22 = ((3*th2.^2-1).*(1./(ys.*R).^2+1./(ys.*R))+th2.^2).*GS;
GEM22 = ((3*th2.^2-1).*(1./(yem.*R).^2+1./(yem.*R))+th2.^2).*GEM;
GS32 = (3./(ys.*R).^2+3./(ys.*R)+1).*th2.*th3.*GS;
GEM32 = (3./(yem.*R).^2+3./(yem.*R)+1).*th2.*th3.*GEM;
GS33 = ((3*th3.^2-1).*(1./(ys.*R).^2+1./(ys.*R))+th3.^2).*GS;
GEM33 = ((3*th3.^2-1).*(1./(yem.*R).^2+1./(yem.*R))+th3.^2).*GEM;
% grad-grad - delta_{ij} is coded on the fly

% particle velocity
KSSvm = -s.*chi.*dems/Gfr;
v11 = 0;
v21S = -KSSvm.*GS3;
v21EM = KSSvm.*GEM3;
v31 = KSSvm.*(GS2-GEM2);
% electric field
KSem = (s.^2.*rhoc./Gfr-ys.^2).*dems;
KEMem = (s.^2.*rhoc./Gfr-yem.^2).*dems;
E11 = 0;
E21 = -KSem.*GS3-KEMem.*GEM3;
E31 = KSem.*GS2-KEMem.*GEM2;
% filtration velocity
KSwm = L.*(s.^2.*rhoB/Gfr-ys.^2).*dems;
KEMwm = L.*(s.^2.*rhoB/Gfr-yem.^2).*dems;
w11 = 0;
w21S = -KSwm.*GS3;
w21EM = KEMwm.*GEM3;
w31 = KSwm.*GS2-KEMwm.*GEM2;
% acoustic pressure is zero
% magnetic field
KSmm = (s.^2.*rhoc/Gfr-ys.^2).*dems;
KEMmm = (s.^2.*rhoc/Gfr-yem.^2).*dems;
NSSmm = -chi.^2.*s.*dems/Gfr;
H11 = KSmm.*(ys.^2.*GS11./zeta-etae.*GS) - ...
      KEMmm.*(yem.^2.*GEM11./zeta-etae.*GEM) + NSSmm.*(GS - GEM);
H21 = (KSmm.*ys.^2.*GS21 - KEMmm.*yem.^2.*GEM21)./zeta;
H31 = (KSmm.*ys.^2.*GS31 - KEMmm.*yem.^2.*GEM31)./zeta;
% bulk stress field
KSStm = chi.*dems;
tau11 = 0;
tau21 = KSStm.*(ys.^2.*GS31-yem.^2.*GEM31);
tau31 = -KSStm.*(ys.^2.*GS21-yem.^2.*GEM21);
tau22 = 2*KSStm.*(ys.^2.*GS32-yem.^2.*GEM32);
tau32 = -KSStm.*(ys.^2.*(GS22-GS33)-yem.^2.*(GEM22-GEM33));
tau33 = -2*KSStm.*(ys.^2.*GS32-yem.^2.*GEM32);
% compute and plot some time-fields
tv21S=2*real(ifft(v21S.*wav,mf*nf,3))/dt;
tv21EM=2*real(ifft(v21EM.*wav,mf*nf,3))/dt;
tE21=2*real(ifft(E21.*wav,mf*nf,3))/dt;
tw21S=2*real(ifft(w21S.*wav,mf*nf,3))/dt;
tw21EM=2*real(ifft(w21EM.*wav,mf*nf,3))/dt;
tH11=2*real(ifft(H11.*wav,mf*nf,3))/dt;
ttau21=2*real(ifft(tau21.*wav,mf*nf,3))/dt;
it=nf/2;
figure(1)
imagesc(xx1,xx2,squeeze(tv21S(:,:,it)));colorbar;
colormap(1-gray)
set(gca,'fontsize',18,'plotboxaspectratio',[1 1 1])
xlabel('x distance from source [m]')
ylabel('y distance from source [m]')
title(['S-wave of G_{11}^{ve}(x,y,z=50m,t=98ms)'])
print -deps2 Svm.eps
figure(2)
imagesc(xx1,xx2,squeeze(tv21EM(:,:,it)));colorbar;
colormap(1-gray)
set(gca,'fontsize',18,'plotboxaspectratio',[1 1 1])
xlabel('x distance from source [m]')
ylabel('y distance from source [m]')
title(['EM field of G_{11}^{ve}(x,y,z=50m,t=98 ms)'])
print -deps2 EMvm.eps
figure(3)
imagesc(xx1,xx2,squeeze(tw21S(:,:,it)));colorbar;
colormap(1-gray)
set(gca,'fontsize',18,'plotboxaspectratio',[1 1 1])
xlabel('x distance from source [m]')
ylabel('y distance from source [m]')
title(['S-wave of G_{11}^{we}(x,y,z=50m,t=98ms)'])
print -deps2 Swm.eps
figure(4)
imagesc(xx1,xx2,squeeze(tw21EM(:,:,it)));colorbar;
colormap(1-gray)
set(gca,'fontsize',18,'plotboxaspectratio',[1 1 1])
xlabel('x distance from source [m]')
ylabel('y distance from source [m]')
title(['EM field of G_{11}^{we}(x,y,z=50m,t=98ms)'])
print -deps2 EMwm.eps

