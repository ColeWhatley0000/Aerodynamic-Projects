clc
clear

%% initial values
% number of panels
n = 128
r = 1
% alpha = 30/360*2*pi;
 alpha = pi*(7/6) ;
 U = 4
%% Geometry 
[x, y, X, Y, ptheta] = cylPanelGeometry(n, r)
[At, An, Cn1, Cn2, Ct1, Ct2, Cp, Ut, S] = cylPanelCoef(x, y, X, Y, ptheta, n, alpha, U)

thetacp = atan2(y,x)

figure(1)
plot(x,y)
hold on
quiver(x,y, (cos(thetacp)).*abs(Cp), (sin(thetacp)).*abs(Cp),'b')
hold off
axis equal

figure(2)
plot(thetacp,Cp)
title("coefficient of pressure and \theta")
figure(3)
plot(thetacp, Ut)
title("Tangential Velocity and \theta")
figure(4)
plot(x,Cp)
figure(5)
plot(x,Ut)
