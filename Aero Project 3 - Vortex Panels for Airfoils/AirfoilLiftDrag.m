function [CL, ratio] = AirfoilLiftDrag(x,y,th,S,alpha,Cp)

n_xy = [cos(th+pi/2); sin(th+ pi/2)]
f = -n_xy.*S.*Cp
F = [sum(f(1,:)), sum(f(2,:))]
rot = [cos(alpha), sin(alpha);...
    -sin(alpha), cos(alpha)]

FDL = rot*F'
CL = FDL(2)
CD = FDL(1)
ratio = CL./(CD)
end