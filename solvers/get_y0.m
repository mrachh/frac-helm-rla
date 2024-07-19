function [y0] = get_y0(X,Y,V,ckb)
theta_in = pi/3;
y0 = ckb*exp(1i*X*cos(theta_in)*ckb+1i*Y*sin(theta_in)*ckb).*V;
end