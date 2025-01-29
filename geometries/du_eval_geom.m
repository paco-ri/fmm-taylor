function [x, y, z] = du_eval_geom(u, v, modes)

m = modes{1};
nmodes = size(m);
n = modes{2};
D = modes{3};
rpiz = zeros(size(u));
drpizdu = zeros(size(u));
for k = 1:nmodes
    rpiz = rpiz + D(k)*exp(-1i*(m(k)*v+n(k)*u));
    drpizdu = drpizdu + D(k)*exp(-1i*(m(k)*v+n(k)*u)).*(-1i*n(k));
end
rpiz = exp(1i*v).*rpiz;
drpizdu = exp(1i*v).*rpiz;
r = real(rpiz);
drdu = real(drpizdu);
z = imag(drpizdu);
x = drdu.*cos(u) - r.*sin(u);
y = drdu.*sin(u) + r.*cos(u);
 
end