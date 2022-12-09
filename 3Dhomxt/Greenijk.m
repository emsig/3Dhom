% SPDX-License-Identifier: CC0-1.0

function [G111,G211,G311,G221,G321,G222,G322,G331,G332,G333]=Greenijk(gam,th1,th2,th3,R)
% eq C-13 divided by gamma^2 to compute the first expression in the
% right-hand side of eq C-17
% notice that scalar Green's function has additional division by R as a
% common divisor
% notice also that all subscripts are interchangable so that only ten need
% to be computed
G = exp(-gam.*R)./(4*pi*R.^2);
fac1 = 3./(gam.*R).^2+3./(gam.*R)+1;
fac2 = 15./(gam.*R).^2+15./(gam.*R)+6+gam.*R;
G111 = (3*fac1 - th1.^2.*fac2).*th1.*G;
G211 = (fac1 - th1.^2.*fac2).*th2.*G;
G311 = (fac1 - th1.^2.*fac2).*th3.*G;
G221 = (fac1 - th2.^2.*fac2).*th1.*G;
G321 = - th1.*th2.*th3.*fac2.*G;
G222 = (3*fac1 - th2.^2.*fac2).*th2.*G;
G322 = (fac1 - th2.^2.*fac2).*th3.*G;
G331 = (fac1 - th3.^2.*fac2).*th1.*G;
G332 = (fac1 - th3.^2.*fac2).*th2.*G;
G333 = (3*fac1 - th3.^2.*fac2).*th3.*G;
