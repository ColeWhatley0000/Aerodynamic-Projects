clc
clear
close all

%% initial geometry
M = flip(readmatrix('airfoil_data.xlsx'),1)
% airfoil for NACA 0012
n = max(size(M)) - 1
X1 = M(:,1);
    X1 = X1'
Y1 = M(:,2);
    Y1 = Y1'
% airfoil for NACA 4412
X2 = M(:,3);
    X2 = X2'
Y2 = M(:,4);
    Y2 = Y2'

    % plotting initial geometry
    f = figure(15)
    f.Position = [100 100 600 400]
    plot(X1,Y1, 'k')
    hold on
    plot(X1,Y1*0, 'r')
    hold off
    axis equal
    grid on

    h = figure(16)
    h.Position = [100 100 600 400]
    plot(X2,Y2, 'k')
    hold on
    plot(X2,Y2*0, 'r')
    hold off
    axis equal
    grid on


%% Initial Conditions
alpha = deg2rad([-16 -8 0 8 16])
n = numel(X1)-1
%% First Figure - del Cp curves
% NACA-0012
figure(1)
subplot(1, 2, 1)
for i = 1:numel(alpha)
    [x,y,S,th, Ut, Cp]  = panelCoef(X1,Y1, n, alpha(i))
    dCp = Cp(1:n/2)-Cp(n:-1:n/2+1)
    hold on
    plot(x(1:n/2), dCp)
    hold off
end

grid on
legend('\alpha = -16 deg.','\alpha = -8 deg.', '\alpha = 0 deg.', '\alpha = 8 deg.', '\alpha = 16 deg.')
xlabel('cord position x/c')
ylabel('\Delta C_p')
% NACA-4412
subplot(1,2,2)
for i = 1:numel(alpha)
    [x,y,S,th, Ut, Cp]  = panelCoef(X2,Y2, n, alpha(i))
    dCp = Cp(1:n/2)-Cp(n:-1:n/2+1)
    hold on
    plot(x(1:n/2), dCp)
    hold off
end

grid on
legend('\alpha = -16 deg.','\alpha = -8 deg.', '\alpha = 0 deg.', '\alpha = 8 deg.', '\alpha = 16 deg.')
xlabel('cord position x/c')
ylabel('\Delta C_p')


%% Second Figure - Lift Shenanigans

alphad = -16:0.5:16
alpha = deg2rad(alphad)
%NACA-0012

for i = 1:numel(alpha)
    [x,y,S,th, Ut, Cp]  = panelCoef(X1,Y1, n, alpha(i))
    [CL(i), ratio(i)] = AirfoilLiftDrag(x,y,th,S,alpha(i),Cp)
    [pu0(i,:), pl0(i,:)] = thwaitesAirfoil (x,y,Ut,S, n)
end

% determining zero lift angle
for i = 1:numel(alpha)-1
    if sign(CL(i+1)) ~= sign(CL(i))
        al_0 = -CL(i)/(CL(i+1)-CL(i))*(alpha(i+1)-alpha(i)) + alpha(i);

        break
    end
end



% reading xfoil data
Pol_a = readmatrix('NACA-0012_polars.xlsx')

% plotting Coefficient of Lift vs alpha
figure(2)
subplot(1,2,1)
plot(alphad, CL, 'b')
hold on
plot(Pol_a(:,1), Pol_a(:,2), 'r')
plot(alphad, 2*pi*(alpha-al_0))
plot(rad2deg(al_0), 0, 'kx')
hold off
% a e s t h e t i c
grid on
legend('Numerical', 'Xfoil Prediction (Re = 10^6)', 'Analytical', 'zero lift angle (numerical)')
ylabel('Lift Coefficient')
xlabel('airfoil angle of attack \alpha (degrees)')

% plotting efficiency vs alpha
eff = ratio
figure(3)
subplot(1,2,1)
plot(alphad, eff)
hold on
plot(Pol_a(:,1), Pol_a(:,2)./Pol_a(:,4))
hold off

grid on
xlabel('Angle of Attack, \alpha (degrees)')
ylabel('Aerodynamic Efficiency (C_L/C_D)')

% NACA-4412
for i = 1:numel(alpha)
    [x,y,S,th, Ut, Cp]  = panelCoef(X2,Y2, n, alpha(i))
    [CL(i), ratio(i)] = AirfoilLiftDrag(x,y,th,S,alpha(i),Cp)
    [pu4(i,:), pl4(i,:)] = thwaitesAirfoil (x,y,Ut,S, n)
end

% determining zero lift angle
for i = 1:numel(alpha)-1
    if sign(CL(i+1)) ~= sign(CL(i))
        al_0 = -CL(i)/(CL(i+1)-CL(i))*(alpha(i+1)-alpha(i)) + alpha(i);
        break
    end
end

Pol_b = readmatrix('NACA-4412_polars.xlsx')

% plotting
figure(2)
subplot(1,2,2)
plot(alphad, CL, 'b')
hold on
plot(Pol_b(:,1), Pol_b(:,2), 'r')
plot(alphad, 2*pi*(alpha-al_0))
plot(rad2deg(al_0), 0, 'kx')
hold off
% a e s t h e t i c
grid on
legend('Numerical', 'Xfoil Prediction (Re = 10^6)', 'Analytical', 'zero lift angle (numerical)')
ylabel('Lift Coefficient')
xlabel('airfoil angle of attack \alpha (degrees)')

% plotting efficiency vs alpha
eff = ratio
figure(3)
subplot(1,2,2)
plot(alphad, eff)
hold on
plot(Pol_b(:,1), Pol_b(:,2)./Pol_b(:,4))
hold off

grid on
xlabel('Angle of Attack, \alpha (degrees)')
ylabel('Aerodynamic Efficiency (C_L/C_D)')

%% plot 3 - separation points
figure(4)
plot(alphad, pu0(:,1))
hold on
plot(alphad, pu4(:,1))
plot(alphad, ones(numel(alphad),1)*.2)
hold off

grid on
legend('NACA-0012', 'NACA-4412')
xlabel('angle of attack \alpha (degrees)')
ylabel('upper separation point')



