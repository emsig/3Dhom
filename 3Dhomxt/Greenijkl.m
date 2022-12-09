% SPDX-License-Identifier: CC0-1.0

function [G1111,G2111,G3111,G2211,G3211,G3311,G2221,G3221,G3321,...
    G3331,G2222,G3222,G3322,G3332,G3333]=Greenijkl(gam,th1,th2,th3,R)
% eq C-14 divided by gamma^4 to compute the first expression in the
% right-hand side of eq C-20; eq C-19 can be computed calling this
% function and the function greenij.m and multiplying the results by gamma^2
% notice that all subscripts can be interchanged and hence only 15
% different combinations exist
G = exp(-gam.*R)./(4*pi*R);
fac1 = 3./(gam.*R).^4+3./(gam.*R).^3+1./(gam.*R).^2;
fac2 = 15./(gam.*R).^4+15./(gam.*R).^3+6./(gam.*R).^2+1./(gam.*R);
fac3 = 105./(gam.*R).^4+105./(gam.*R).^3+45./(gam.*R).^2+10./(gam.*R)+1;
G1111 = (3*fac1 - 6*th1.^2.*fac2 + th1.^4.*fac3).*G;
G2111 = (-3*th1.*th2.*fac2 + th2.*th1.^3.*fac3).*G;
G3111 = (-3*th1.*th3.*fac2 + th3.*th1.^3.*fac3).*G;
G2211 = (fac1 - (th2.^2+th1.^2).*fac2 + th2.^2.*th1.^2.*fac3).*G;
G3211 = (-th2.*th3.*fac2 + th3.*th2.*th1.^2.*fac3).*G;
G3311 = (fac1 - (th3.^2+th1.^2).*fac2 + th3.^2.*th1.^2.*fac3).*G;
G2221 = (-3*th1.*th2.*fac2 + th1.*th2.^3.*fac3).*G;
G3221 = (-th1.*th3.*fac2 + th1.*th2.^2.*th3.*fac3).*G;
G3321 = (-th1.*th2.*fac2 + th1.*th2.*th3.^2.*fac3).*G;
G3331 = (-3*th1.*th3.*fac2 + th1.*th3.^3.*fac3).*G;
G2222 = (3*fac1 - 6*th2.^2.*fac2 + th2.^4.*fac3).*G;
G3222 = (-3*th2.*th3.*fac2 + th3.*th2.^3.*fac3).*G;
G3322 = (fac1 - (th2.^2+th3.^2).*fac2 + th2.^2.*th3.^2.*fac3).*G;
G3332 = (-3*th2.*th3.*fac2 + th2.*th3.^3.*fac3).*G;
G3333 = (3*fac1 - 6*th3.^2.*fac2 + th3.^4.*fac3).*G;
