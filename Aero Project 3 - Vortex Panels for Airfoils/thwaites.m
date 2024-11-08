function [theta, stagx, stagy] = thwaites (x,y,Ut,S, n)

%% finding stagnation points

% testing for stag pt at first/last index
if sign(Ut(n)) ~= sign(Ut(1))
    indH(1) = 1;
    indL(1) = n;
    j = 2;
else 
    j = 1;
end

for i = 1:n-1
    if Ut(i) == 0
        ind = i
    else
        if sign(Ut(i)) ~= sign(Ut(i+1))
            indL(j) = i;
            indH(j) = i+1
            j = j+1;
        end
    end
end

for i = 1:numel(indL)
    stagx(i) = mean([x(indL(i)), x(indH(i))])
    stagy(i) = mean([y(indL(i)), y(indH(i))])
end

%% velocity and distance vectors
Uu(1) = 0
Ul(1) = 0

su(1) = 0
su(2) = S(n)/2

sl(1) = 0
sl(2) = S(1)/2

Uu(2:n+1) = Ut(n:-1:1)
Ul(2:n+1) = Ut(1:n)

 for i = 2:n-1
     su(i+1) = su(i) + S(n+2-i)/2 + S(n+1-i)/2
     sl(i+1) = sl(i) + S(i-1)/2 + S(i)/2
 end

%% Upper separation
% upper portion {n>i>stagH}
if Uu(2) < 0
    Uu = abs(Uu)
end
uind =(n-indL(2))
Uu5 = Uu(1:numel(Uu)-1).^5
Uu6 = Uu(1:numel(Uu)-1).^6

dUu = diff(Uu(1:uind))
dsu = diff(su(1:uind))
dudxu = dUu./dsu
int_u = cumtrapz(su(1:uind-1),Uu5(1:uind-1))
K = (0.45./Uu6(1:uind-1)).*int_u.*dudxu + 0.09;

for i = 1:numel(K)-1
    if sign(K(i)) ~= sign(K(i+1))
        sepu = [i, i+1];
    end
end

sepu = n-sepu;
th = atan2(y(sepu),x(sepu));
xsep = x(sepu)
ysep = y(sepu)

sepu = n-sepu
theta(1) = -K(sepu(1))/(K(sepu(2))-K(sepu(1)))*(th(2)-th(1)) + th(1)
xs_up = -K(sepu(1))/(K(sepu(2))-K(sepu(1)))*(xsep(1)-xsep(2)) + xsep(2)
ys_up = -K(sepu(1))/(K(sepu(2))-K(sepu(1)))*(ysep(1)-ysep(2)) + ysep(2)


% lower portion {1<i<stagL}
if Ul(2) < 0
    Ul = abs(Ul)
end
uind =(indL(2))
Ul5 = Ul(1:numel(Ul)-1).^5
Ul6 = Ul(1:numel(Ul)-1).^6

dUl = diff(Ul(1:uind))
dsl = diff(sl(1:uind))
dudxl = dUl./dsl
int_l = cumtrapz(sl(1:uind-1),Ul5(1:uind-1))
K = (0.45./Ul6(1:uind-1)).*int_l.*dudxl + 0.09;

for i = 1:numel(K)-1
    if sign(K(i)) ~= sign(K(i+1))
        sepu = [i, i+1];
    end
end

th = atan2(y(sepu),x(sepu))
theta(2) = -K(sepu(1))/(K(sepu(2))-K(sepu(1)))*(th(2)-th(1)) + th(1)







end







