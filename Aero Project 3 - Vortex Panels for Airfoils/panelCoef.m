function [x,y,S,th, Ut, Cp, CLD]  = panelCoef(X,Y, n, alpha)

%% Geometry Generation
for i = 1:n
    dx(i) = X(i+1) - X(i);
    dy(i) = Y(i+1) - Y(i);
    th(i) = atan2(dy(i), dx(i));
    S(i) = (dx(i)^2 + dy(i)^2)^0.5

    x(i) = mean([X(i), X(i+1)]);
    y(i) = mean([Y(i), Y(i+1)]);
    
end

RHS = sin(th'-alpha);

%% Coefficients

for i = 1:n
    for j = 1:n
        % shortcuts
        dx = x(i)-X(j);
        dy = y(i)-Y(j);

        dth = th(i) - th(j);
        dth2 = th(i) - 2*th(j);
        % Alphabet Coefficients
        A = -(dx*cos(th(j)) + dy*sin(th(j)));
        B = dx^2 + dy^2;
        C = sin(dth);
        D = cos(dth);
        E = dx*sin(th(j)) - dy*cos(th(j));
        F = log(1 + (S(j)^2 + 2*A*S(j))/B);
        G = atan2(E*S(j), B+A*S(j));

        P = dx*sin(dth2) + dy*cos(dth2);
        Q = dx*cos(dth2) - dy*sin(dth2);

        % C Coefficients
        Cn2(i,j) = D + .5*Q*F/S(j) - (A*C + D*E)*G/S(j);
        Cn1(i,j) = 0.5*D*F + C*G - Cn2(i,j);
        
        if i == j
            Ct2(i,i) = pi/2;
            Ct1(i,i) = pi/2;
        else
            Ct2(i,j) = C + 0.5*P*F/S(j) + (A*D - C*E)*G/S(j);
            Ct1(i,j) = 0.5*C*F - D*G - Ct2(i,j);
        end
    end
end

%% Influence Coefficients

An = zeros(n,n+1);
At = An;

% first Collumn
An(:,1) = Cn1(:,1);
At(:,1) = Ct1(:,1);

% Collumns 2:n
for i = 1:n
    for j = 2:n
        An(i,j) = Cn1(i,j) + Cn2(i,j-1);
        At(i,j) = Ct1(i,j) + Ct2(i,j-1);
    end
end

% Collumn n+1
An(:,n+1) = Cn2(:,n);
At(:,n+1) = Ct2(:,n);

%% Final Calculations
g = An\RHS;

Ut = cos(th - alpha);

for i = 1:n
    for j = 1:n
        Ut(i) = Ut(i) + At(i,j)*g(j);
    end
end

Cp = 1 - Ut.^2;


