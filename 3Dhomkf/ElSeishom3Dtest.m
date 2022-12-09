% SPDX-License-Identifier: CC0-1.0

% defining all the variables
clear all;
clc;
close all;
% predetermine fontsize
fs=14;
% Setting the parameters for a source wave
nf              =   64;
mf              =   4; % multiplication factor to have sinc-interpolation in time-domain
% Define the center frequency of the source waveform
fc              =   40;
% Define frequency axis
freq            =   linspace(0,4*fc,nf);
freq(1)         =   1e-6; % To prevent unnecessary problems upon division by zero. Can't be too close to zero, though, 
                         % since some variables either blow up or vanish in comparison to the rest of the entries in the frequency vector.
s               =   2i*pi*freq(nf/4)+2*pi*freq(2);
% compute the frequency step size
df              =   freq(3)-freq(2);
dt              =   1/(mf*nf*df);
t               =   linspace(-(mf-1)*nf-1,(mf-1)*nf,mf*nf)*dt;


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
L               =   0; %because omega << omega_c (critical frequency)
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
k               =   k0./(sqrt(1+4.*s./(similaritypar.*omegac))+s./omegac);% frequency-dependent dynamic permeability
rhoE            =   eta./(s.*k); % effective density
rhoc            =   rhoB-((rhof.*rhof)./rhoE); % Complex Density rhoC
zeta            =   sigmaM+s.*mu0; 
etae            =   sigmaE+s.*epsilon;
varsigma        =   etae - s.*rhoE.*L.^2;
chi             =   s.*rhof*L;

ypf     =   s./sqrt(2).*sqrt((rhoB.*M-2.*rhof.*C+rhoE.*H.*etae./varsigma)./(H.*M-C.^2) - sqrt(((rhoB.*M-2.*rhof.*C+rhoE.*H.*etae./varsigma)./(H.*M-C.^2)).^2 + 4.*(rhof.^2-rhoB.*rhoE.*etae./varsigma)./(H.*M-C.^2)));
yps     =   s./sqrt(2).*sqrt((rhoB.*M-2.*rhof.*C+rhoE.*H.*etae./varsigma)./(H.*M-C.^2) + sqrt(((rhoB.*M-2.*rhof.*C+rhoE.*H.*etae./varsigma)./(H.*M-C.^2)).^2 + 4.*(rhof.^2-rhoB.*rhoE.*etae./varsigma)./(H.*M-C.^2)));
ys      =   1/sqrt(2).*sqrt(s.^2.*rhoc./Gfr+zeta.*etae - sqrt((s.^2.*rhoc./Gfr-zeta.*etae).^2-4.*s.^3.*zeta.*(rhof.*L).^2./Gfr));
yem     =   1/sqrt(2).*sqrt(s.^2.*rhoc./Gfr+zeta.*etae + sqrt((s.^2.*rhoc./Gfr-zeta.*etae).^2-4.*s.^3.*zeta.*(rhof.*L).^2./Gfr));
dems    =   1./(yem.^2-ys.^2);            
dpspf   =   1./(yps.^2-ypf.^2);     

nx=512;
x1=linspace(0,500,nx/2+1);
dx=x1(2);
dk=2*pi/(dx*nx);
kk=(0:nx/2+1)*dk;
kk2=(0:nx/2)*dk;
[k1,k2]=meshgrid(kk,kk2);
nn=64;
for ik = 1:nx/nn
    k3 = (ik-1)*dk*nn +1e-3;
    kap = sqrt(k1.^2 + k2.^2 + k3.^2);
    kmat11 = (k1./kap).^2;
    kmat12 = k1.*k2./kap.^2;
    kmat13 = k1.*k3./kap.^2;
    kmat21 = k2.*k1./kap.^2;
    kmat22 = (k2./kap).^2;
    kmat23 = k2.*k3./kap.^2;
    kmat31 = k3.*k1./kap.^2;
    kmat32 = k3.*k2./kap.^2;
    kmat33 = (k3./kap).^2;
    GSS = (yem.^2-ys.^2)./((kap.^2+yem.^2).*(kap.^2+ys.^2));
    GPP = (yps.^2-ypf.^2)./((kap.^2+ypf.^2).*(kap.^2+yps.^2));
% force acting on the bulk
    KSSvf = -s.*(zeta.*etae+kap.^2).*dems/Gfr;
    KPPvf = s.*M.*(s.^2.*rhoE.*etae./(M.*varsigma)+kap.^2).*dpspf/(H*M-C^2);
    KSSef = -s.*chi.*zeta.*dems/Gfr;
    KPPef = s.^2.*rhoE.*L*C.*(s.^2*rhof./C+kap.^2).*dpspf./((H*M-C^2)*varsigma);
    KSSwf = s.*rhof.*(zeta.*varsigma+kap.^2).*dems./(rhoE*Gfr);
    KPPwf = -s*C.*(s.^2*rhof./C+kap.^2).*dpspf./(H*M-C^2);
    KPPpf = -s.^2.*(rhoE.*etae*C./varsigma-rhof*M).*dpspf/(H*M-C^2);
    KSSmf = s.*chi.*dems/Gfr;
    KSStf = -(zeta.*etae+kap.^2).*dems;
    KPPtf = 2*Gfr*M*(s.^2.*rhoE.*etae./(M*varsigma)+kap.^2).*dpspf/(H*M-C^2);
    NPPtf = (s.^2.*(rhoE.*H.*etae./varsigma-rhof*C)/(H*M-C^2)+kap.^2).*dpspf;
    v11 = KSSvf.*GSS.*(kmat11-1)+KPPvf.*GPP.*kmat11;
    v21 = KSSvf.*GSS.*kmat21+KPPvf.*GPP.*kmat21;
    v31 = KSSvf.*GSS.*kmat31+KPPvf.*GPP.*kmat31;
    E11 = KSSef.*GSS.*(kmat11-1)+KPPef.*GPP.*kmat11;
    E21 = KSSef.*GSS.*kmat21+KPPef.*GPP.*kmat21;
    E31 = KSSef.*GSS.*kmat31+KPPef.*GPP.*kmat31;
    w11 = KSSwf.*GSS.*(kmat11-1)+KPPwf.*GPP.*kmat11;
    w21 = KSSwf.*GSS.*kmat21+KPPwf.*GPP.*kmat21;
    w31 = KSSwf.*GSS.*kmat31+KPPwf.*GPP.*kmat31;
    H11 = 0;
    H21 = KSSmf.*1i.*k3.*GSS;
    H31 = -KSSmf.*1i.*k2.*GSS;
    p1 = -1i*k1.*KPPpf.*GPP;
    tau11 = -2*KSStf.*GSS.*(kmat11-1).*1i.*k1-KPPtf.*GPP.*(kmat11-1).*1i.*k1...
        - NPPtf.*GPP.*1i.*k1;
    tau21 = -KSStf.*GSS.*(2*kmat11-1).*1i.*k2-KPPtf.*GPP.*kmat11.*1i.*k2;
    tau22 = -2*KSStf.*GSS.*kmat22.*1i.*k1-KPPtf.*GPP.*(kmat22-1).*1i.*k1...
        - NPPtf.*GPP.*1i.*k1;
    tau31 = -KSStf.*GSS.*(2*kmat11-1).*1i.*k3-KPPtf.*GPP.*kmat11.*1i.*k3;
    tau32 = -2*KSStf.*GSS.*kmat12.*1i.*k3-KPPtf.*GPP.*kmat12.*1i.*k3;
    tau33 = -2*KSStf.*GSS.*kmat33.*1i.*k1-KPPtf.*GPP.*(kmat33-1).*1i.*k1...
        - NPPtf.*GPP.*1i.*k1;
% electric current source
    KSSve = KSSef;
    KPPve = KPPef;
    KSSee = zeta.*(s.^2.*rhoc/Gfr+kap.^2).*dems;
    KPPee = (s.*rhoE.*L./varsigma).^2.*s.*H.*(s.^2.*rhoB./H+kap.^2).*dpspf./(H*M-C^2);
    KSSwe = zeta.*L.*(s.^2.*rhoB/Gfr+kap.^2).*dems;
    KPPwe = -s.^2.*rhoE.*L.*H.*(s.^2.*rhoB/H+kap.^2).*dpspf./((H*M-C^2).*varsigma);
    KSSme = -(s.^2.*rhoc./Gfr+kap.^2).*dems;
    KPPpe = s.*rhoE.*L.*(s.^2.*(rhoB*M-rhof*C)/(H*M-C^2)+kap.^2).*dpspf./varsigma;
    KSSte = -chi.*zeta.*dems;
    KPPte = 2*s.*rhoE.*L.*C.*Gfr.*(s.^2.*rhof/C+kap.^2).*dpspf./((H*M-C^2).*varsigma);
    NPPte = s.^3.*rhoE.*L.*(rhof*H-rhoB*C).*dpspf./((H*M-C^2).*varsigma);
    v11 = KSSve.*GSS.*(kmat11-1)+KPPve.*GPP.*kmat11;
    v21 = KSSve.*GSS.*kmat21+KPPve.*GPP.*kmat21;
    v31 = KSSve.*GSS.*kmat31+KPPve.*GPP.*kmat31;
    E11 = KSSee.*GSS.*(kmat11-1)+(KPPee.*GPP-1/varsigma).*kmat11;
    E21 = KSSee.*GSS.*kmat21+(KPPee.*GPP-1/varsigma).*kmat21;
    E31 = KSSee.*GSS.*kmat31+(KPPee.*GPP-1/varsigma).*kmat31;
    w11 = KSSwe.*GSS.*(kmat11-1)+KPPwe.*GPP.*kmat11;
    w21 = KSSwe.*GSS.*kmat21+KPPwe.*GPP.*kmat21;
    w31 = KSSwe.*GSS.*kmat31+KPPwe.*GPP.*kmat31;
    H11 = 0;
    H21 = KSSme.*1i.*k3.*GSS;
    H31 = -KSSme.*1i.*k2.*GSS;
    p1 = -KPPpe.*GPP.*1i.*k1;
    tau11 = -2*KSSte.*GSS.*(kmat11-1).*1i.*k1-KPPte.*GPP.*(kmat11-1).*1i.*k1...
        - NPPte.*GPP.*1i.*k1;
    tau21 = -KSSte.*GSS.*(2*kmat11-1).*1i.*k2-KPPte.*GPP.*kmat11.*1i.*k2;
    tau22 = -2*KSSte.*GSS.*kmat22.*1i.*k1-KPPte.*GPP.*(kmat22-1).*1i.*k1...
        - NPPte.*GPP.*1i.*k1;
    tau31 = -KSSte.*GSS.*(2*kmat11-1).*1i.*k3-KPPte.*GPP.*kmat11.*1i.*k3;
    tau32 = -2*KSSte.*GSS.*kmat12.*1i.*k3-KPPte.*GPP.*kmat12.*1i.*k3;
    tau33 = -2*KSSte.*GSS.*kmat33.*1i.*k1-KPPte.*GPP.*(kmat33-1).*1i.*k1...
        - NPPte.*GPP.*1i.*k1;
% 
    Maxw111 = 1i*(k2.*H31-k3.*H21)+varsigma.*E11+s.*rhoE.*L.*w11+1;
    Maxw121 = 1i*(k3.*H11-k1.*H31)+varsigma.*E21+s.*rhoE.*L.*w21;
    Maxw131 = 1i*(k1.*H21-k2.*H11)+varsigma.*E31+s.*rhoE.*L.*w31;
    Maxw211 = -1i*(k2.*E31-k3.*E21)+zeta.*H11;
    Maxw221 = -1i*(k3.*E11-k1.*E31)+zeta.*H21;
    Maxw231 = -1i*(k1.*E21-k2.*E11)+zeta.*H31;
    Newt1=s.*(rhoB.*v11+rhof.*w11)+1i.*(k1.*tau11+k2.*tau21+k3.*tau31);
    Newt2=s.*(rhoB.*v21+rhof.*w21)+1i.*(k1.*tau21+k2.*tau22+k3.*tau32);
    Newt3=s.*(rhoB.*v31+rhof.*w31)+1i.*(k1.*tau31+k2.*tau32+k3.*tau33);
    Newf1=s.*(rhof.*v11+rhoE.*(w11-L*E11))-1i.*k1.*p1;
    Newf2=s.*(rhof.*v21+rhoE.*(w21-L*E21))-1i.*k2.*p1;
    Newf3=s.*(rhof.*v31+rhoE.*(w31-L*E31))-1i.*k3.*p1;
    Hookf=s.*p1-C*1i*(k1.*v11+k2.*v21+k3.*v31)-M*1i*(k1.*w11+k2.*w21+k3.*w31);
    Hke11=s*tau11+(H-2*Gfr)*1i*(k1.*v11+k2.*v21+k3.*v31)+Gfr*2i*(k1.*v11) ...
        + C*1i*(k1.*w11+k2.*w21+k3.*w31);
    Hke21=s*tau21+Gfr*1i*(k1.*v21+k2.*v11);
    Hke22=s*tau22+(H-2*Gfr)*1i*(k1.*v11+k2.*v21+k3.*v31)+Gfr*2i*(k2.*v21) ...
        + C*1i*(k1.*w11+k2.*w21+k3.*w31);
    Hke31=s*tau31+Gfr*1i*(k1.*v31+k3.*v11);
    Hke32=s*tau32+Gfr*1i*(k2.*v31+k3.*v21);
    Hke33=s*tau33+(H-2*Gfr)*1i*(k1.*v11+k2.*v21+k3.*v31)+Gfr*2i*(k3.*v31) ...
        + C*1i*(k1.*w11+k2.*w21+k3.*w31);
    Hkf11=s.*p1-1i*(k1.*v11+k2.*v21+k3.*v31)*C-M*1i*(k1.*w11+k2.*w21+k3.*w31);
    if ik > 1
        close(1)
    end
    figure(1)
    subplot(3,3,1);imagesc(log10(abs(Hke11)));colorbar;title('Hke_{21}');
    subplot(3,3,2);imagesc(log10(abs(Hke21)));colorbar;title('Hke_{22}');
    subplot(3,3,3);imagesc(log10(abs(Hke31)));colorbar;title('Hke_{33}');
    subplot(3,3,4);imagesc(log10(abs(Newt1)));colorbar;title('NWt_{11}');
    subplot(3,3,5);imagesc(log10(abs(Newt2)));colorbar;title('NWt_{21}');
    subplot(3,3,6);imagesc(log10(abs(Newt3)));colorbar;title('NWt_{31}');
    subplot(3,3,7);imagesc(log10(abs(Newf1)));colorbar;title('NWf_{11}');
    subplot(3,3,8);imagesc(log10(abs(Newf2)));colorbar;title('NWf_{21}');
    subplot(3,3,9);imagesc(log10(abs(Newf3)));colorbar;title('NWf_{31}');
    pause(0.1)
    figure(2)
    subplot(2,3,1);imagesc(log(abs(Maxw111)));colorbar;title('ME1_{11}');
    subplot(2,3,2);imagesc(log(abs(Maxw121)));colorbar;title('ME1_{21}');
    subplot(2,3,3);imagesc(log(abs(Maxw131)));colorbar;title('ME1_{31}');
    subplot(2,3,4);imagesc(log(abs(Maxw211)));colorbar;title('ME2_{11}');
    subplot(2,3,5);imagesc(log(abs(Maxw221)));colorbar;title('ME2_{21}');
    subplot(2,3,6);imagesc(log(abs(Maxw231)));colorbar;title('ME2_{31}');
    pause(0.1)
% force acting on the fluid
    KSSvff = s.*rhof.*(zeta.*varsigma+kap.^2).*dems./(Gfr.*rhoE);
    KPPvff = -s.*C.*(s.^2.*rhof/C+kap.^2).*dpspf./(H*M-C^2);
    KSSeff = KSSwe;
    KPPeff = KPPwe;
    KSSwff = -(s.^2.*rhoB/Gfr+kap.^2).*(zeta.*varsigma+kap.^2).*dems./(s.*rhoE);
    KPPwff = s.*H.*(s.^2.*rhoB/H+kap.^2).*dpspf./(H*M-C^2);
    KSSmff = -L.*(s.^2.*rhoB/Gfr+kap.^2).*dems;
    KPPpff = -(s.^2*(rhoB.*M-rhof*C)./(H*M-C^2)+kap.^2).*dpspf;
    KSStff = rhof.*(zeta.*varsigma+kap.^2).*dems./rhoE;
    KPPtff = -2*Gfr*C.*(s.^2*rhof/C+kap.^2).*dpspf./(H*M-C^2);
    NPPtff = -s.^2*(rhof*H-rhoB*C).*dpspf./(H*M-C^2);
    v11 = KSSvff.*GSS.*(kmat11-1)+KPPvff.*GPP.*kmat11;
    v21 = KSSvff.*GSS.*kmat21+KPPvff.*GPP.*kmat21;
    v31 = KSSvff.*GSS.*kmat31+KPPvff.*GPP.*kmat31;
    E11 = KSSeff.*GSS.*(kmat11-1)+KPPeff.*GPP.*kmat11;
    E21 = KSSeff.*GSS.*kmat21+KPPeff.*GPP.*kmat21;
    E31 = KSSeff.*GSS.*kmat31+KPPeff.*GPP.*kmat31;
    w11 = KSSwff.*GSS.*(kmat11-1)+KPPwff.*GPP.*kmat11;
    w21 = KSSwff.*GSS.*kmat21+KPPwff.*GPP.*kmat21;
    w31 = KSSwff.*GSS.*kmat31+KPPwff.*GPP.*kmat31;
    H11 = 0;
    H21 = KSSmff.*1i.*k3.*GSS;
    H31 = -KSSmff.*1i.*k2.*GSS;
    p1 = -KPPpff.*GPP.*1i.*k1;
    tau11 = -2*KSStff.*GSS.*(kmat11-1).*1i.*k1-KPPtff.*GPP.*(kmat11-1).*1i.*k1...
        - NPPtff.*GPP.*1i.*k1;
    tau21 = -KSStff.*GSS.*(2*kmat11-1).*1i.*k2-KPPtff.*GPP.*kmat11.*1i.*k2;
    tau22 = -2*KSStff.*GSS.*kmat22.*1i.*k1-KPPtff.*GPP.*(kmat22-1).*1i.*k1...
        - NPPtff.*GPP.*1i.*k1;
    tau31 = -KSSte.*GSS.*(2*kmat11-1).*1i.*k3-KPPtff.*GPP.*kmat11.*1i.*k3;
    tau32 = -2*KSStff.*GSS.*kmat12.*1i.*k3-KPPtff.*GPP.*kmat12.*1i.*k3;
    tau33 = -2*KSStff.*GSS.*kmat33.*1i.*k1-KPPtff.*GPP.*(kmat33-1).*1i.*k1...
        - NPPtff.*GPP.*1i.*k1;
% volume injection rate source
    KPPvq = -s.^2.*(C*rhoE.*etae./varsigma-M*rhof).*dpspf/(H*M-C^2);
    KPPeq = s.*rhoE.*L.*(s.^2*(rhoB.*M-rhof*C)./(H*M-C^2)+kap.^2).*dpspf./varsigma;
    KPPwq = KPPpff;
    KPPmq = 0;
    KPPpq = s.*rhoB.*etae.*(s.^2.*(rhoc+chi.*rhof.*L/etae)*M/(H*M-C^2)+...
        kap.^2).*dpspf./varsigma;
    KPPtq = 2*Gfr*s.*(rhoE.*etae*C./varsigma-rhof*M).*kap.^2.*dpspf;
    NPPtq = -s.*rhof.*(s.^2*C.*(rhoB*rhoE*etae./(rhof*varsigma)-rhof)...
        ./(H*M-C^2)+kap.^2).*dpspf;
    v11 = -1i*k1.*KPPvq.*GPP;
    v21 = -1i*k2.*KPPvq.*GPP;
    v31 = -1i*k3.*KPPvq.*GPP;
    E11 = -1i*k1.*KPPeq.*GPP;
    E21 = -1i*k2.*KPPeq.*GPP;
    E31 = -1i*k3.*KPPeq.*GPP;
    w11 = -1i*k1.*KPPwq.*GPP;
    w21 = -1i*k2.*KPPwq.*GPP;
    w31 = -1i*k3.*KPPwq.*GPP;
    H11 = 0;
    H21 = 0;
    H31 = 0;
    p = KPPpq.*GPP;
    tau11 = KPPtq.*GPP.*((k1./kap).^2-1)+NPPtq;
    tau21 = KPPtq.*GPP.*k2.*k1./kap.^2;
    tau31 = KPPtq.*GPP.*k3.*k1./kap.^2;
%     Maxw111 = 1i*(k2.*H31-k3.*H21)+varsigma.*E11+s.*rhoE.*L.*w11;
%     Maxw121 = 1i*(k3.*H11-k1.*H31)+varsigma.*E21+s.*rhoE.*L.*w21;
%     Maxw131 = 1i*(k1.*H21-k2.*H11)+varsigma.*E31+s.*rhoE.*L.*w31;
%     Maxw211 = -1i*(k2.*E31-k3.*E21)+zeta.*H11;
%     Maxw221 = -1i*(k3.*E11-k1.*E31)+zeta.*H21;
%     Maxw231 = -1i*(k1.*E21-k2.*E11)+zeta.*H31;
% magnetic current source
    KSSvm = -KSSmf;
    KSSem = -KSSme;
    KSSwm = -KSSmff;
    KSSmm = -(s.^2*rhoc/Gfr+kap.^2).*dems;
    NSSmm = -chi.^2.*s.*dems/Gfr;
    KPPpm = 0;
    KSStm = chi.*dems;
    v11 = 0;
    v21 = 1i*k3.*KSSvm.*GSS;
    v31 = -1i*k2.*KSSvm.*GSS;
    E11 = 0;
    E21 = 1i*k3.*KSSem.*GSS;
    E31 = -1i*k2.*KSSem.*GSS;
    w11 = 0;
    w21 = 1i*k3.*KSSwm.*GSS;
    w31 = -1i*k2.*KSSwm.*GSS;
    H11 = KSSmm.*GSS.*(etae+k1.^2./zeta)+NSSmm.*GSS;
    H21 = KSSmm.*GSS.*k1.*k2./zeta;
    H31 = KSSmm.*GSS.*k1.*k3./zeta;
    p1 = 0;
    tau11 = 0;
    tau21 = -k1.*k3.*KSStm.*GSS;
    tau31 = k1.*k2.*KSStm.*GSS;
% deformation rate source h_{11}
    KSSvh = -KSStf;
    KPPvh = -KPPtf;
    NPPvh = -NPPtf;
    KSSeh = -KSSte;
    KPPeh = -KPPte;
    NPPeh = -NPPte;
    KSSwh = -KSStff;
    KPPwh = -KPPtff;
    NPPwh = -NPPtff;
    KSSmh = KSStm;
    KPPph = -KPPtq;
    NPPph = -NPPtq;
    KSSth = Gfr.*(zeta.*etae+kap.^2).*dems./s;
    KPPth = 2*Gfr.*(kap.^4+kap.^2.*s.^2.*(rhoE*H*etae./varsigma-rhof*C)...
        ./(H*M-C^2)).*dpspf./s;
    NPPth = H*(kap.^4+kap.^2.*s.^2.*(rhoE*H*etae./varsigma-2*rhof*C+rhoB*C^2/H)...
        ./(H*M-C^2)).*dpspf./s;
    QPPth = 4*Gfr*(M*kap.^4+kap.^2.*s.^2.*rhoE.*etae./varsigma)./(s.*(H*M-C^2));
    v11 = 2*KSSvh.*GSS.*(1-kmat11).*1i.*k1+KPPvh.*GPP.*(1-kmat11).*1i.*k1...
        - NPPvh.*GPP.*1i.*k1;
    v21 = -2*KSSvh.*GSS.*kmat11.*1i.*k2+KPPvh.*GPP.*(1-kmat11).*1i.*k2...
        - NPPvh.*GPP.*1i.*k2;
    v31 = -2*KSSvh.*GSS.*kmat11.*1i.*k3+KPPvh.*GPP.*(1-kmat11).*1i.*k3...
        - NPPvh.*GPP.*1i.*k3;
    E11 = 2*KSSeh.*GSS.*(1-kmat11).*1i.*k1+KPPeh.*GPP.*(1-(k1./kap).^2).*1i.*k1...
        - NPPeh.*GPP.*1i.*k1;
    E21 = -2*KSSeh.*GSS.*(k1./kap).^2.*1i.*k2+KPPeh.*GPP.*(1-(k1./kap).^2).*1i.*k2...
        - NPPeh.*GPP.*1i.*k2;
    E31 = -2*KSSeh.*GSS.*(k1./kap).^2.*1i.*k3+KPPeh.*GPP.*(1-(k1./kap).^2).*1i.*k3...
        - NPPeh.*GPP.*1i.*k3;
    w11 = 2*KSSwh.*GSS.*(1-(k1./kap).^2).*1i.*k1+KPPwh.*GPP.*(1-(k1./kap).^2).*1i.*k1...
        - NPPwh.*GPP.*1i.*k1;
    w21 = -2*KSSwh.*GSS.*(k1./kap).^2.*1i.*k2+KPPwh.*GPP.*(1-(k1./kap).^2).*1i.*k2...
        - NPPwh.*GPP.*1i.*k2;
    w31 = -2*KSSwh.*GSS.*(k1./kap).^2.*1i.*k3+KPPwh.*GPP.*(1-(k1./kap).^2).*1i.*k3...
        - NPPwh.*GPP.*1i.*k3;
    H11 = 0;
    H21 = -2*k1.*k3.*KSSmh.*GSS;
    H31 = 2*k1.*k2.*KSSmh.*GSS;
    p1 = KPPph.*GPP.*((k1./kap).^2-1)+NPPph.*GPP;
    tau11 = 2*KPPth.*GPP.*((k1./kap).^2-1)+NPPth.*GPP.*((k1./kap).^2-1).^2 ...
       + QPPth.*GPP+4*KSSth.*GSS.*(1-(k1./kap).^2).*k1.^2+H./s;
    tau21 = KPPth.*GPP.*k2.*k1./kap.^2-NPPth.*k1.*k2.*(1-(k1./kap).^2)./kap.^2 ...
        -KSSth.*GSS.*(2*k1.^3.*k2./kap.^2+2*(1-(k1./kap).^2).*k1.*k2);
    tau31 = KPPth.*GPP.*k3.*k1./kap.^2-NPPth.*k1.*k3.*(1-(k1./kap).^2)./kap.^2 ...
        -KSSth.*GSS.*(2*k1.^3.*k3./kap.^2+2*(1-(k1./kap).^2).*k1.*k3);
%     Maxw111 = 1i*(k2.*H31-k3.*H21)+varsigma.*E11+s.*rhoE.*L.*w11;
%     Maxw121 = 1i*(k3.*H11-k1.*H31)+varsigma.*E21+s.*rhoE.*L.*w21;
%     Maxw131 = 1i*(k1.*H21-k2.*H11)+varsigma.*E31+s.*rhoE.*L.*w31;
%     Maxw211 = -1i*(k2.*E31-k3.*E21)+zeta.*H11;
%     Maxw221 = -1i*(k3.*E11-k1.*E31)+zeta.*H21;
%     Maxw231 = -1i*(k1.*E21-k2.*E11)+zeta.*H31;
%     n1=abs(k1.*tau11+k2.*tau21+k3.*tau31);
%     n2=abs(k1.*tau21+k2.*tau22+k3.*tau32)+1e-10;
%     n3=abs(k1.*tau31+k2.*tau32+k3.*tau33)+1e-10;
%     Newt1=(s.*(rhoB.*v11+rhof.*w11)+1i.*(k1.*tau11+k2.*tau21+k3.*tau31));
%     Newt2=(s.*(rhoB.*v21+rhof.*w21)+1i.*(k1.*tau21+k2.*tau22+k3.*tau32));
%     Newt3=(s.*(rhoB.*v31+rhof.*w31)+1i.*(k1.*tau31+k2.*tau32+k3.*tau33));
%     Newf1=(s.*rhof.*v11+s.*rhoE.*(w11-L*E11)-1i.*k1.*p1)./abs(s.*rhof.*v11);
%     Newf2=(s.*rhof.*v21+s.*rhoE.*(w21-L*E21)-1i.*k2.*p1)./abs(s.*rhof.*v21);
%     Newf3=(s.*rhof.*v31+s.*rhoE.*(w31-L*E31)-1i.*k3.*p1)./abs(s.*rhof.*v31);
%     Hookf=s.*p1-C*1i*(k1.*v11+k2.*v21+k3.*v31)-M*1i*(k1.*w11+k2.*w21+k3.*w31);
%     Hke11=(s*tau11+(H-2*Gfr)*1i*(k1.*v11+k2.*v21+k3.*v31)+Gfr*2i*(k1.*v11) ...
%         + C*1i*(k1.*w11+k2.*w21+k3.*w31))./abs(s*tau11);
%     Hke21=(s*tau21+Gfr*1i*(k1.*v21+k2.*v11))./abs(s*tau21);
%     Hke22=(s*tau22+(H-2*Gfr)*1i*(k1.*v11+k2.*v21+k3.*v31)+Gfr*2i*(k2.*v21) ...
%         + C*1i*(k1.*w11+k2.*w21+k3.*w31))./abs(s*tau22);
%     Hke31=(s*tau31+Gfr*1i*(k1.*v31+k3.*v11))./abs(s*tau31);
%     Hke32=(s*tau32+Gfr*1i*(k2.*v31+k3.*v21))./abs(s*tau32);
%     Hke33=(s*tau33+(H-2*Gfr)*1i*(k1.*v11+k2.*v21+k3.*v31)+Gfr*2i*(k3.*v31) ...
%         + C*1i*(k1.*w11+k2.*w21+k3.*w31))./abs(s*tau33);
%     Hkf11=(s.*p1-1i*(k1.*v11+k2.*v21+k3.*v31)*C-M*1i*(k1.*w11+k2.*w21+k3.*w31)) ...
%         ./abs(M*(k1.*w11+k2.*w21+k3.*w31));
    

%     figure(1)
%     subplot(2,3,1);imagesc(log(abs(Maxw111)));colorbar;title('ME1_{11}');
%     subplot(2,3,2);imagesc(log(abs(Maxw121)));colorbar;title('ME1_{21}');
%     subplot(2,3,3);imagesc(log(abs(Maxw131)));colorbar;title('ME1_{31}');
%     subplot(2,3,4);imagesc(log(abs(Maxw211)));colorbar;title('ME2_{11}');
%     subplot(2,3,5);imagesc(log(abs(Maxw221)));colorbar;title('ME2_{21}');
%     subplot(2,3,6);imagesc(log(abs(Maxw231)));colorbar;title('ME2_{31}');
%     pause(0.1)
%     figure(1)
%     subplot(3,3,1);imagesc(log10(abs(Hke11)));colorbar;title('Hke_{21}');
%     subplot(3,3,2);imagesc(log10(abs(Hke22)));colorbar;title('Hke_{22}');
%     subplot(3,3,3);imagesc(log10(abs(Hke33)));colorbar;title('Hke_{33}');
%     subplot(3,3,4);imagesc(log10(abs(Newt1)));colorbar;title('NWt_{11}');
%     subplot(3,3,5);imagesc(log10(abs(Newt2)));colorbar;title('NWt_{21}');
%     subplot(3,3,6);imagesc(log10(abs(Newt3)));colorbar;title('NWt_{31}');
%     subplot(3,3,7);imagesc(log10(abs(Newf1)));colorbar;title('NWf_{11}');
%     subplot(3,3,8);imagesc(log10(abs(Newf2)));colorbar;title('NWf_{21}');
%     subplot(3,3,9);imagesc(log10(abs(Newf3)));colorbar;title('NWf_{31}');
%     pause(0.1)
end
