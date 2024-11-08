function [pu, pl] = thwaitesAirfoil (x,y,Ut,S, n)

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
            indH(j) = i+1;
            j = j+1;
        end
    end
end

for i = 1:numel(indL)
    stagx(i) = mean([x(indL(i)), x(indH(i))]);
    stagy(i) = mean([y(indL(i)), y(indH(i))]);
end

%% creating necessary vectors
% S values for integral

if indH(1) == n/2+1
    indu = ([indH(2):1:n, 1:n/2]);
    indl = ([indL(2):-1:n/2+1]);
end

if indH(1) == 1
    indu = ([1:1:n/2])
    indl = ([n:-1:n/2+1])
end

if indH(1) ~= n && indH(1) ~= n/2+1
    indu = ([indH(1):1:n/2])
    indl = ([indL(1):-1:1, n:-1:n/2+1])
end

Uu = [0, Ut(indu)];
Ul = [0, Ut(indl)];
Ut
Su = S(indu);
su = [0, Su(1)/2];
for i = 2:numel(indu)-1
    su(i+1) = su(i) + Su(i)/2 + Su(i-1)/2;
end

Sl = S(indl);
sl = [0, Sl(1)/2];

for i = 2:numel(indl)-1
    sl(i+1) = sl(i) + Sl(i)/2 + Sl(i-1)/2
end

%% upper separation
Uu = abs(Uu)
Uu5 = Uu(1:numel(Uu)-2).^5;
Uu6 = Uu(1:numel(Uu)-2).^6;
dUu = diff(Uu(1:numel(Uu)-1));
dsu = diff(su);
dudsu = dUu./dsu;
intu = cumtrapz(su(1:numel(indu)-1), Uu5);

K = (.45./Uu6) .* intu .* dudsu + 0.09;
for i = 2:numel(indu)-1;
    if sign(K(i+1)) ~= sign(K(i))
        sepu = i;
        break
    end
end

%interpolation
pu(1) = -K(sepu)/(K(sepu)-K(sepu+1)) * (x(indu(sepu+1)) - x(indu(sepu))) + x(indu(sepu));
pu(2) = -K(sepu)/(K(sepu)-K(sepu+1)) * (y(indu(sepu+1)) - y(indu(sepu))) + y(indu(sepu));

%% lower separation
Ul = abs(Ul)
Ul5 = Ul(1:numel(Ul)-2).^5;
Ul6 = Ul(1:numel(Ul)-2).^6;
dUl = diff(Ul(1:numel(Ul)-1));
dsl = diff(sl);
dudsl = dUl./dsl;
intl = cumtrapz(sl(1:numel(indl)-1), Ul5);

K = (.45./Ul6) .* intl .* dudsl + 0.09;
for i = 2:numel(indu)-1
    if sign(K(i+1)) ~= sign(K(i))
        sepu = i;
        break
    end
end

%interpolation
pl(1) = -K(sepu)/(K(sepu)-K(sepu+1)) * (x(indl(sepu+1)) - x(indl(sepu))) + x(indl(sepu));
pl(2) = -K(sepu)/(K(sepu)-K(sepu+1)) * (y(indl(sepu+1)) - y(indl(sepu))) + y(indl(sepu));


end
