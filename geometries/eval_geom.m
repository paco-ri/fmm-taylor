function [x, y, z] = eval_geom(u, v, modes)

m = modes{1};
nmodes = size(m);
n = modes{2};
D = modes{3};
rpiz = zeros(size(u));
for k = 1:nmodes
    rpiz = rpiz + D(k)*exp(-1i*(m(k)*v+n(k)*u));
end
rpiz = exp(1i*v).*rpiz;
r = real(rpiz);
z = imag(rpiz);
x = r.*cos(u);
y = r.*sin(u);
 
end