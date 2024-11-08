function [At, An, Cn1, Cn2, Ct1, Ct2, Cp, Ut, S] = cylPanelCoef(x, y, X, Y, ptheta, n, alpha, U)
    % Panel Geometry Coefficients
    % [X(i,j) = a(i) - b(j)] == [X = a' - b]
        dx = x'-X(1:n) 
        dy = y'-Y(1:n) 
        dptheta = ptheta' - ptheta % difference between panel angles of i and j
        pth_2diff = ptheta' - 2*ptheta

    S = sqrt((X(1)-X(2))^2+(Y(1)-Y(2))^2);

    %% Coefficients
    A = -dx.*cos(ptheta) - dy.*sin(ptheta);
    B = dx.^2 + dy.^2;
    C = sin(dptheta);
    D = cos(dptheta);
    E = dx.*sin(ptheta) - dy.*cos(ptheta);
    F = log(1 + (S^2 + 2*A*S)./B);
    G = atan2((E*S) , (B + A*S));
    P = dx.*sin(pth_2diff) + dy.*cos(pth_2diff);
    Q = dx.*cos(pth_2diff) - dy.*sin(pth_2diff);

    %% Final C coefficients
    Cn2 = D + 0.5*Q.*F/S - (A.*C + D.*E).*G/S
    Cn1 = 0.5*D.*F + C.*G - Cn2
    
    Ct2 = C + 0.5*P.*F/S + (A.*D - C.*E).*G/S
    Ct1 = 0.5*C.*F - D.*G - Ct2

    for i = 1:n
        Ct1(i,i) = pi/2;
        Ct2(i,i) = pi/2;
    end

    %% Velocity and Pressure Coef Calculations
    ang = ptheta'-alpha
    RHS = U*sin(ang);
    Ut = cos(ang);

    An = zeros(n, n+1);
    At = An;
    An(:,1) = Cn1(:,1);
    At(:,1) = Ct1(:,1);

    An(:,n+1) = Cn2(:,n);
    At(:,n+1) = Ct2(:,n);
    for i = 1:n
        for j = 2:n
            An(i,j) = Cn1(i,j) + Cn2(i,j-1);
            At(i,j) = Ct1(i,j) + Ct2(i,j-1);
        end
    end

    g = An\RHS
    
    Ut = cos(ptheta - alpha);
    
for i = 1:n
    Ut(i) = cos(ptheta(i) - alpha);
    for j = 1:n
        Ut(i) = Ut(i) + At(i,j)*g(j)
    end
    Cp(i) = 1 - Ut(i)^2
end
end


