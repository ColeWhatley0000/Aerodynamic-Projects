function [x,y,S,ptheta, Ut, Cp]  = panelCylCoef(X,Y, n, alpha)

%% initial calculations - control points and panel angles
for i = 1:n
    dx(i) = X(i+1) - X(i);
    dy(i) = Y(i+1) - Y(i);
    S(i) = (dx(i)^2 + dy(i)^2)^0.5

    x(i) = mean([X(i), X(i+1)]);
    y(i) = mean([Y(i), Y(i+1)]);
    
end
    

ptheta = atan2(dy, dx);
RHS = sin(ptheta'-alpha);
%% Coefficient Calculations

% initial matrices

dx = x'- X;
dy = y'- Y;

% Letter coefficients
    B = dx.^2 + dy.^2
    C = sin(ptheta' - ptheta);
    D = cos(ptheta' - ptheta);


for j = 1:n
    A(:,j) = -dx(:,j).*cos(ptheta(j)) - dy(:,j).*sin(ptheta(j));

    E(:,j) = dx(:,j).*sin(ptheta(j)) - dy(:,j).*cos(ptheta(j));
    F(:,j) = log(1 + (S(j)^2 +2*S(j)*A(:,j))./B(:,j));
    G(:,j) = atan2(E(:,j)*S(j) , B(:,j) + A(:,j)*S(j));
end

for i = 1:n
    for j = 1:n
        P(i,j) = dx(i,j)*sin(ptheta(i)-2*ptheta(j)) + dy(i,j)*cos(ptheta(i)-2*ptheta(j));
        Q(i,j) = dx(i,j)*cos(ptheta(i)-2*ptheta(j)) - dy(i,j)*sin(ptheta(i)-2*ptheta(j));
    end
end

% C coefiicients
Cn1 = zeros(n,n);
Cn2 = Cn1;
Ct1 = Cn1;
Ct2 = Cn1;

for j = 1:n
    Cn2(:,j) = D(:,j) + 0.5*Q(:,j).*F(:,j)/S(j)...
        - (A(:,j).*C(:,j) + D(:,j).*E(:,j)).*G(:,j)/S(j);

    Ct2(:,j) = C(:,j) + 0.5*P(:,j).*F(:,j)/S(j)...
        + (A(:,j).*D(:,j) - C(:,j).*E(:,j)).*G(:,j)/S(j);
end
    Cn1 = 0.5*D.*F + C.*G - Cn2;
    Ct1 = 0.5*C.*F - D.*G - Ct2;

for i = 1:n
    Ct1(i,i) = pi/2;
    Ct2(i,i) = pi/2;
end


% A coefficients

At = zeros(n,n+1);
An = At;

An(:,1) = Cn1(:,1);
An(:,n+1) = Cn2(:,n);

At(:,1) = Ct1(:,1);
At(:,n+1) = Ct2(:,n);

for j = 2:n
    An(:,j) = Cn1(:,j) + Cn2(:,j-1);
    At(:,j) = Ct1(:,j) + Ct2(:,j-1);
end

g = An\RHS;

Ut = cos(ptheta - alpha);

for i = 1:n
    Ut(i) = cos(ptheta(i) - alpha);
    for j = 1:n
        Ut(i) = Ut(i) + At(i,j)*g(j);
    end
    Cp(i) = 1 - Ut(i)^2;
end






end