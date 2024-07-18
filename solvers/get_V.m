function [V] = get_V(X,Y,amp,dia,rad,n)
    if (nargin <6)
        rad = 4;
        dia = 1;
        amp = 1;
        n   = 5;
    end
    V = zeros(size(X));
    for i=1:n
        theta = (i-1)/n*2*pi;
        x0    = rad*cos(theta);
        y0    = rad*sin(theta);
        V     = V + amp*exp(-(X-x0).^2/(2*dia^2)-(Y-y0).^2/(2*dia^2)); 
    end
end

