% SPDX-License-Identifier: CC0-1.0

function [G1,G2,G3]=Greeni(gam,th1,th2,th3,R)
% eq C-15 from C-11
G = exp(-gam.*R)./(4*pi*R);
G1 = -(1./R+gam).*th1.*G;
G2 = -(1./R+gam).*th2.*G;
G3 = -(1./R+gam).*th3.*G;
