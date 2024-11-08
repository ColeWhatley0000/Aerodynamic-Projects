function [x, y, X, Y, ptheta] = cylPanelGeometry(n, r)

    theta(1) = 0
    X(1) = r*cos(theta(1))
    Y(1) = r*sin(theta(1))
    
    for i = 1:n
        % edges and angles at edges of panels

        % Coordinates of the panels
        theta(i+1) = theta(i) - 2*pi/n
        X(i+1) = r*cos(theta(i+1))
        Y(i+1) = r*sin(theta(i+1))

        % panel angle
        dY(i) = Y(i+1)-Y(i);
        dX(i) = X(i+1)-X(i)
        ptheta(i) = atan2(dY(i), dX(i));

        % control points of panels
        x(i) = (X(i)+ X(i+1))/2
        y(i) = (Y(i)+ Y(i+1))/2
        % length of panels
        S(i) = (dX(i)^2 + dY(i)^2)^0.5
    end

end