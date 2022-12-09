% SPDX-License-Identifier: CC0-1.0

function [G11,G21,G31,G22,G32,G33]=Greenij(gam,th1,th2,th3,R)
% eq C-16 from C-12
G = exp(-gam.*R)./(4*pi*R);
G11 = ((3*th1.^2-1).*(1./(gam.*R).^2+1./(gam.*R))+th1.^2).*G;
G21 = (3./(gam.*R).^2+3./(gam.*R)+1).*th1.*th2.*G;
G31 = (3./(gam.*R).^2+3./(gam.*R)+1).*th1.*th3.*G;
G22 = ((3*th2.^2-1).*(1./(gam.*R).^2+1./(gam.*R))+th2.^2).*G;
G32 = (3./(gam.*R).^2+3./(gam.*R)+1).*th2.*th3.*G;
G33 = ((3*th3.^2-1).*(1./(gam.*R).^2+1./(gam.*R))+th3.^2).*G;
