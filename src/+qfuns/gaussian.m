function [V] = gaussian(X, Y, amp, dia, rad, n)
    if (nargin <6)
        rad = 0;
        dia = 0.2;
        amp = 1;
        n   = 1;
    end
    V = zeros(size(X));
    for i=1:n
        theta = (i-1)/n*2*pi;
        x0    = rad*cos(theta) + 0.04;
        y0    = rad*sin(theta) + 0.03;
        V     = V + amp*exp(-(X-x0).^2/(2*dia^2)-(Y-y0).^2/(2*dia^2)); 
    end
end

