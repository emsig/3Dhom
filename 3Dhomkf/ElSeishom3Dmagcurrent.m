% SPDX-License-Identifier: CC0-1.0
%
% defining all the variables
clear all;
clc;
close all;
% Define the center frequency of the source waveform
fc              =   40;
% for this test of the Green's fucntions we use the basic equations and a
% single frrequency value
s               =   2i*pi*fc;

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
L               =   10000; %because omega << omega_c (critical frequency)
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

% define the spherical wavenumebrs of eqs B-9 and B-10.
ys      =   1/sqrt(2).*sqrt(s.^2.*rhoc./Gfr+zeta.*etae - sqrt((s.^2.*rhoc./Gfr-zeta.*etae).^2-4.*s.^3.*zeta.*(rhof.*L).^2./Gfr));
yem     =   1/sqrt(2).*sqrt(s.^2.*rhoc./Gfr+zeta.*etae + sqrt((s.^2.*rhoc./Gfr-zeta.*etae).^2-4.*s.^3.*zeta.*(rhof.*L).^2./Gfr));
% define the transversal factor that occurs in all coefficients
dems    =   1./(yem.^2-ys.^2);            

% define the grid in k-space, we make a square grid in (kx,ky)-plane and a 
% small number of vertical wavenumber values to test the accuracy of
% Green's functions for different values of kz
nx=512;
x1=linspace(0,500,nx/2+1);
dx=x1(2);
dk=2*pi/(dx*nx);
kk=(0:nx/2+1)*dk;
[k1,k2]=meshgrid(kk,kk);
nn=64;
for ik = 1:nx/nn
% compute the vertical wavenumber
    k3 = (ik-1)*dk*nn/2 +1e-3;
% compute the spherical wavenumber
    kap = sqrt(k1.^2 + k2.^2 + k3.^2);
% compute the scalar transversal Green's functions
    GSS = 1./(kap.^2+ys.^2) - 1./(kap.^2+yem.^2);
% force acting on the bulk eq 39 and reciprocity
    KSSmf = s.*chi.*dems/Gfr;
    KSSvm = -KSSmf;
% electric current source eq 48 and reciprocity
    KSSme = -(s.^2.*rhoc./Gfr+kap.^2).*dems;
    KSSem = -KSSme;
% force acting on the fluid eq 55 and reciprocity
    KSSmff = -L.*(s.^2.*rhoB/Gfr+kap.^2).*dems;
    KSSwm = -KSSmff;
% magnetic current source
    KSSmm = (s.^2*rhoc/Gfr+kap.^2).*dems;
    NSSmm = -chi.^2.*s.*dems/Gfr;
    KSStm = chi.*dems;
% particle velocity as given in reciprocal form in eq 15
    v11 = 0;
    v21 = 1i*k3.*KSSvm.*GSS;
    v31 = -1i*k2.*KSSvm.*GSS;
% electric field as given in reciprocal form in eq 15
    E11 = 0;
    E21 = 1i*k3.*KSSem.*GSS;
    E31 = -1i*k2.*KSSem.*GSS;
% filtration velocity as given in reciprocal form in eq 15
    w11 = 0;
    w21 = 1i*k3.*KSSwm.*GSS;
    w31 = -1i*k2.*KSSwm.*GSS;
% acoustic pressure as given in eq 18
    p1 = 0;
% magnetic field in eq 20
    H11 = -KSSmm.*GSS.*(etae+k1.^2./zeta)+NSSmm.*GSS;
    H21 = -KSSmm.*GSS.*k1.*k2./zeta;
    H31 = -KSSmm.*GSS.*k1.*k3./zeta;
% buk stress as given in eq 21
    tau11 = 0;
    tau21 = -k1.*k3.*KSStm.*GSS;
    tau31 = k1.*k2.*KSStm.*GSS;
    tau22 = -2*k2.*k3.*KSStm.*GSS;
    tau32 = (k2.^2-k3^2).*KSStm.*GSS;
    tau33 = 2*k2.*k3.*KSStm.*GSS;
% test eq 3
    Maxw111 = 1i*(k2.*H31-k3.*H21)+varsigma.*E11+s.*rhoE.*L.*w11;
    Maxw121 = 1i*(k3.*H11-k1.*H31)+varsigma.*E21+s.*rhoE.*L.*w21;
    Maxw131 = 1i*(k1.*H21-k2.*H11)+varsigma.*E31+s.*rhoE.*L.*w31;
% test eq 4
    Maxw211 = -1i*(k2.*E31-k3.*E21)+zeta.*H11+1;
    Maxw221 = -1i*(k3.*E11-k1.*E31)+zeta.*H21;
    Maxw231 = -1i*(k1.*E21-k2.*E11)+zeta.*H31;
% test eq 7
    Newt1=s.*(rhoB.*v11+rhof.*w11)+1i.*(k1.*tau11+k2.*tau21+k3.*tau31);
    Newt2=s.*(rhoB.*v21+rhof.*w21)+1i.*(k1.*tau21+k2.*tau22+k3.*tau32);
    Newt3=s.*(rhoB.*v31+rhof.*w31)+1i.*(k1.*tau31+k2.*tau32+k3.*tau33);
% test eq 8
    Newf1=s.*(rhof.*v11+rhoE.*(w11-L*E11))-1i.*k1.*p1;
    Newf2=s.*(rhof.*v21+rhoE.*(w21-L*E21))-1i.*k2.*p1;
    Newf3=s.*(rhof.*v31+rhoE.*(w31-L*E31))-1i.*k3.*p1;
% test eq 9
    Hke11=(s*tau11+(H-2*Gfr)*1i*(k1.*v11+k2.*v21+k3.*v31)+Gfr*2i*(k1.*v11) ...
        + C*1i*(k1.*w11+k2.*w21+k3.*w31))/C;
    Hke21=s*tau21+Gfr*1i*(k1.*v21+k2.*v11);
    Hke22=(s*tau22+(H-2*Gfr)*1i*(k1.*v11+k2.*v21+k3.*v31)+Gfr*2i*(k2.*v21) ...
        + C*1i*(k1.*w11+k2.*w21+k3.*w31))/C;
    Hke31=s*tau31+Gfr*1i*(k1.*v31+k3.*v11);
    Hke32=s*tau32+Gfr*1i*(k2.*v31+k3.*v21);
    Hke33=(s*tau33+(H-2*Gfr)*1i*(k1.*v11+k2.*v21+k3.*v31)+Gfr*2i*(k3.*v31) ...
        + C*1i*(k1.*w11+k2.*w21+k3.*w31))/C;
% test eq 10
    Hkf11=s.*p1-1i*(k1.*v11+k2.*v21+k3.*v31)*C-M*1i*(k1.*w11+k2.*w21+k3.*w31);
% ploit test results
    if ik > 1
        close(1)
        close(2)
    end
    figure(1)
    subplot(3,3,1);imagesc(log10(abs(Hke11)));colorbar;title('Hke_{11}');
    subplot(3,3,2);imagesc(log10(abs(Hke21)));colorbar;title('Hke_{21}');
    subplot(3,3,3);imagesc(log10(abs(Hke31)));colorbar;title('Hke_{31}');
    subplot(3,3,4);imagesc(log10(abs(Newt1)));colorbar;title('NWt_{11}');
    subplot(3,3,5);imagesc(log10(abs(Newt2)));colorbar;title('NWt_{21}');
    subplot(3,3,6);imagesc(log10(abs(Newt3)));colorbar;title('NWt_{31}');
    subplot(3,3,7);imagesc(log10(abs(Newf1)));colorbar;title('NWf_{11}');
    subplot(3,3,8);imagesc(log10(abs(Newf2)));colorbar;title('NWf_{21}');
    subplot(3,3,9);imagesc(log10(abs(Newf3)));colorbar;title('NWf_{31}');
    pause(0.1)
    figure(2)
    subplot(3,3,1);imagesc(log10(abs(Hke22)));colorbar;title('Hke_{22}');
    subplot(3,3,2);imagesc(log10(abs(Hke33)));colorbar;title('Hke_{33}');
    subplot(3,3,3);imagesc(log10(abs(Hke32)));colorbar;title('Hke_{32}');
    subplot(3,3,4);imagesc(log10(abs(Maxw111)));colorbar;title('ME1_{11}');
    subplot(3,3,5);imagesc(log10(abs(Maxw121)));colorbar;title('ME1_{21}');
    subplot(3,3,6);imagesc(log10(abs(Maxw131)));colorbar;title('ME1_{31}');
    subplot(3,3,7);imagesc(log10(abs(Maxw211)));colorbar;title('ME2_{11}');
    subplot(3,3,8);imagesc(log10(abs(Maxw221)));colorbar;title('ME2_{21}');
    subplot(3,3,9);imagesc(log10(abs(Maxw231)));colorbar;title('ME2_{31}');
    pause(0.1)
end
