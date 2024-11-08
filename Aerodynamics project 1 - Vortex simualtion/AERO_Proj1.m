clc
clear
close all

%% Parameters
b = 2;              %parameter for x & y values
S = 3;              % vortex strength
dt = .001;           % time step
N = 5000;            % number of time steps
m = 0;
G = [0.2, 0.4, 0.6]
for M = 1:3;
    g = G(M)
 %% Test Case Inputs
             
 P = [10 20 50]                                 % number of vortices
 n = P(M)
 z = n/2
 pt_dist = linspace(-1,1,n);
 xi_n = pt_dist*b/2;
 %xi_n = [-b/2, -b/4, b/4, b/2, -b/2, -b/4, b/4, b/2];                  % initial x values
 yi_n = 0*pt_dist
% yi_n = [b/2, b/2, b/2, b/2, -b/2, -b/2, -b/2, -b/2];                            % inital y values
 S_n = pt_dist
 for i= 1:n
    j = i-(z+1)
    S_n(i) = ((2*j+1)/z)/sqrt(1-((j+0.5)^2)/z^2)/10
end
%  S_n = [-S, -g*S, g*S, S, S, g*S, -g*S, -S ];                        % Circulation Strength Values
%  n = numel(S_n)

dx = @(xi, xj, yi, yj, Sj) (Sj.*(yj-yi)) ./ (2*pi .* ((xi-xj).^2+(yi-yj).^2));
dy = @(xi, xj, yi, yj, Sj) (Sj.*(xi-xj)) ./ (2*pi .* ((xi-xj).^2+(yi-yj).^2));

%% Initializing
xn = zeros(N,n);
yn = zeros(N,n);
        u = 0
        v = 0
xn(1,:) = xi_n;
yn(1,:) = yi_n;
t(1) = 0;
i = 1
numpt = 301
q = 2;
%% Generalized Case (Loops and loops and loops and loops and...)


x = linspace(-q*b,b*q, numpt);
y = linspace(-q*b,b*q, numpt);

[X,Y] = meshgrid(x,y);

% Time Loop enterance
for h = 1:N
    
    k1 = zeros(n, 2);
    k2 = zeros(n, 2);
    k3 = zeros(n, 2);
    k4 = zeros(n, 2);

    % k1 Values
    for s = 1:n
        for r = 1:n
            if s ~=r
                k1(s,1) = dx(xn(h,s), xn(h,r), yn(h,s),...
                    yn(h,r), S_n(r)) + k1(s,1);
                k1(s,2) = dy(xn(h,s), xn(h,r), yn(h,s),...
                    yn(h,r), S_n(r)) + k1(s,2);
            end
        end
    end    

    % k2 Values
    for s = 1:n             % looping the vortex being analyzed
        for r = 1:n         % looping the vortex inducing velocity 
            if s ~=r        % ensuring the sum doesnt consider i = j
                % x differential k
                k2(s,1) = dx(...
                    xn(h,s) + k1(s,1)*dt/2,...
                    xn(h,r) + k1(r,1)*dt/2,...
                    yn(h,s) + k1(s,2)*dt/2,... 
                    yn(h,r) + k1(r,2)*dt/2,...
                    S_n(r)) + k2(s,1);
                % y differential k
                k2(s,2) = dy(...
                    xn(h,s) + k1(s,1)*dt/2,...
                    xn(h,r) + k1(r,1)*dt/2,...
                    yn(h,s) + k1(s,2)*dt/2,...
                    yn(h,r) + k1(r,2)*dt/2,...
                    S_n(r)) + k2(s,2);
            end
        end
    end
    
    % k3 Values
    for s = 1:n
        for r = 1:n
            if s ~=r
                % x differential k
                k3(s,1) = dx(xn(h,s) + k2(s,1)*dt/2,...
                    xn(h,r) + k2(r,1)*dt/2,...
                    yn(h,s) + k2(s,2)*dt/2,...
                    yn(h,r) + k2(r,2)*dt/2,...
                    S_n(r)) + k3(s,1);
                % y differential k
                k3(s,2) = dy(xn(h,s) + k2(s,1)*dt/2,...
                    xn(h,r) + k2(r,1)*dt/2,...
                    yn(h,s) + k2(s,2)*dt/2,...
                    yn(h,r) + k2(r,2)*dt/2,...
                    S_n(r)) + k3(s,2);
            end
        end
    end
    
    % k4 Values
    for s = 1:n
        for r = 1:n
            if s ~=r
                % x differential k
                k4(s,1) = dx(xn(h,s) + k3(s,1)*dt,...
                    xn(h,r) + k3(r,1)*dt,...
                    yn(h,s) + k3(s,2)*dt,...
                    yn(h,r) + k3(r,2)*dt,...
                    S_n(r)) + k4(s,1);
                % y differential k
                k4(s,2) = dy(xn(h,s) + k3(s,1)*dt,...
                    xn(h,r) + k3(r,1)*dt,...
                    yn(h,s) + k3(s,2)*dt,...
                    yn(h,r) + k3(r,2)*dt,...
                    S_n(r)) + k4(s,2);
            end
        end
    end

    % updating x and y
    for s = 1:n
        xn(h+1,s) = xn(h,s) + dt/6*(k1(s,1)+2*k2(s,1)+2*k3(s,1)+k4(s,1));
        yn(h+1,s) = yn(h,s) + dt/6*(k1(s,2)+2*k2(s,2)+2*k3(s,2)+k4(s,2));
    end

    index = [1, 2500, 5000]
    if ismember(h,index) == 1    
        for s = 1:n
            u = dx(X,xn(h,s),Y,yn(h,s),S_n(s)) + u;
            v = dy(X,xn(h,s),Y,yn(h,s),S_n(s)) + v;
        end
        Umag = (u.^2 + v.^2).^0.5;
%% plotting Velocity
        % Velocity Contour ------------------------------------------------
        m = m + 1
        figure(1)
        subplot(3, ceil(numel(index)), m)
        contour(X,Y,(Umag),[0:.25:10], 'linewidth', 1) 
            % ln of Umag taken since Umag goes to infinity at r=0
        axis equal
        title(sprintf('t = %f seconds', h*dt))
        j = colorbar;
        set(get(j,'label'),'string','Velocity Magnitude (m/s)');
        % Quiver plot -----------------------------------------------------
        hold on
        % Reduction in quiver frequency
        qv_redux = 10;
        qv_rx = 1:qv_redux:numpt;
        qv_ry = 1:qv_redux:numpt;

        quiver(X(qv_rx, qv_ry),Y(qv_rx, qv_ry),...
            u(qv_rx, qv_ry),v(qv_rx, qv_ry), 'linewidth',1.5);
        axis equal
        axis([min(min(X)), max(max(X)),min(min(Y)), max(max(Y))])
        hold off
        xlabel('x')
        ylabel('y')
        %resetting velocity for next loop
        u = 0;
        v = 0;
        i = i + 1;
    end

end

%% plotting
% contour plot stuff\
sgtitle('Contour and Velocity Vector Field')


% Paths of Vortices
figure(3)
subplot(1, 3, M)
plot(xn(:,1:4), yn(:,1:4),'-')
hold on

plot(xn, yn,'rx');
f = figure(3);
f.Position = [100 100 1550 400];
% for s = 200:200:(N+1)
%     plot(xn(s,:), yn(s,:),'k.');
% end
hold off
title(sprintf('Vortex Path, n = %f', n))
xlabel('X')
ylabel('y')
axis([min(min(xn)), max(max(xn)),min(min(yn)), max(max(yn))]*1.1)
axis equal

end
