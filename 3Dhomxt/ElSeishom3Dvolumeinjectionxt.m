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
% spherical wavenumbers of eqs B-12 and B-13.
ypf     =   s./sqrt(2).*sqrt((rhoB.*M-2.*rhof.*C+rhoE.*H.*etae./varsigma)./(H.*M-C.^2) - sqrt(((rhoB.*M-2.*rhof.*C+rhoE.*H.*etae./varsigma)./(H.*M-C.^2)).^2 + 4.*(rhof.^2-rhoB.*rhoE.*etae./varsigma)./(H.*M-C.^2)));
yps     =   s./sqrt(2).*sqrt((rhoB.*M-2.*rhof.*C+rhoE.*H.*etae./varsigma)./(H.*M-C.^2) + sqrt(((rhoB.*M-2.*rhof.*C+rhoE.*H.*etae./varsigma)./(H.*M-C.^2)).^2 + 4.*(rhof.^2-rhoB.*rhoE.*etae./varsigma)./(H.*M-C.^2)));
% recurring coefficient factors
dpspf   =   1./(yps.^2-ypf.^2);     

% volume injection rate source

% Green's functions of eq C-10
GPs = exp(-yps.*R)./(4*pi*R);
GPf = exp(-ypf.*R)./(4*pi*R);
% grad (equation C-15)
[GPf1,GPf2,GPf3]=Greeni(ypf,th1,th2,th3,R);
[GPs1,GPs2,GPs3]=Greeni(yps,th1,th2,th3,R);
% grad-grad divided by gamma^2 (equation C-16)
[GPf11,GPf21,GPf31,GPf22,GPf32,GPf33]=Greenij(ypf,th1,th2,th3,R);
[GPs11,GPs21,GPs31,GPs22,GPs32,GPs33]=Greenij(yps,th1,th2,th3,R);
% grad-grad - delta_{ij} is coded on the fly

% particle velocity
KPPvq = -s.^2.*(rhoE.*etae*C./varsigma-rhof*M).*dpspf/(H*M-C^2);
v1 = KPPvq.*(GPf1-GPs1);
v2 = KPPvq.*(GPf2-GPs2);
v3 = KPPvq.*(GPf3-GPs3);
% electric field
KPseq = s.*rhoE.*L.*(s.^2.*(rhoB*M-rhof*C)/(H*M-C^2)-yps.^2).*dpspf./varsigma;
KPfeq = s.*rhoE.*L.*(s.^2.*(rhoB*M-rhof*C)/(H*M-C^2)-ypf.^2).*dpspf./varsigma;
E1 = KPfeq.*GPf1-KPseq.*GPs1;
E2 = KPfeq.*GPf2-KPseq.*GPs2;
E3 = KPfeq.*GPf3-KPseq.*GPs3;
% filtration velocity
KPswq = -(s.^2*(rhoB.*M-rhof*C)./(H*M-C^2)-yps.^2).*dpspf;
KPfwq = -(s.^2*(rhoB.*M-rhof*C)./(H*M-C^2)-ypf.^2).*dpspf;
w1 = KPfwq.*GPf1-KPswq.*GPs1;
w2 = KPfwq.*GPf2-KPswq.*GPs2;
w3 = KPfwq.*GPf3-KPswq.*GPs3;
% acoustic pressure
KPspq = (s.^2.*M.*(rhoB-rhof^2.*varsigma./(rhoE.*etae))/(H*M-C^2) ...
        -yps.^2).*s.*rhoE.*etae.*dpspf./varsigma;
KPfpq = (s.^2.*M.*(rhoB-rhof^2.*varsigma./(rhoE.*etae))/(H*M-C^2) ...
        -ypf.^2).*s.*rhoE.*etae.*dpspf./varsigma;
p = KPfpq.*GPf-KPspq.*GPs;
% magnetic field is zero
% bulk stress field
KPPtq = 2*Gfr*s.*(rhoE.*etae*C./varsigma-rhof*M).*dpspf/(H*M-C^2);
NPstq = -s.*(s.^2*C.*(rhoB*rhoE.*etae./varsigma-rhof^2)...
        ./(H*M-C^2)-rhof.*yps.^2).*dpspf;
NPftq = -s.*(s.^2*C.*(rhoB*rhoE.*etae./varsigma-rhof^2)...
        ./(H*M-C^2)-rhof.*ypf.^2).*dpspf;
tau11 = -KPPtq.*(ypf.^2.*(GPf11-GPf)-yps.^2.*(GPs11-GPs)) - (NPftq.*GPf-NPstq.*GPs);
tau21 = -KPPtq.*(ypf.^2.*GPf21-yps.^2.*GPs21);
tau31 = -KPPtq.*(ypf.^2.*GPf31-yps.^2.*GPs31);
tau22 = -KPPtq.*(ypf.^2.*(GPf22-GPf)-yps.^2.*(GPs22-GPs)) - (NPftq.*GPf-NPstq.*GPs);
tau32 = -KPPtq.*(ypf.^2.*GPf32-yps.^2.*GPs32);
tau33 = -KPPtq.*(ypf.^2.*(GPf33-GPf)-yps.^2.*(GPs33-GPs)) - (NPftq.*GPf-NPstq.*GPs);
% compute and plot some time-fields
tv1=2*real(ifft(v1.*wav,mf*nf,3))/dt;
tE1=2*real(ifft(E1.*wav,mf*nf,3))/dt;
tw1=2*real(ifft(w1.*wav,mf*nf,3))/dt;
tp=2*real(ifft(p.*wav,mf*nf,3))/dt;
ttau11=2*real(ifft(tau11.*wav,mf*nf,3))/dt;
for it=1:nf
    imagesc(xx1,xx2,squeeze(tE1(:,:,it)));
    colormap(1-gray)
    colorbar;
    set(gca,'fontsize',18,'plotboxaspectratio',[1 1 1])
    xlabel('x distance from source [m]')
    ylabel('y distance from source [m]')
    pause(0.1)
end
