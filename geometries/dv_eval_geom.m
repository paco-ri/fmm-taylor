function [x, y, z] = dv_eval_geom(u, v, modes)

m = modes{1};
nmodes = size(m);
n = modes{2};
D = modes{3};
rpiz = zeros(size(u));
drpizdv = zeros(size(u));
for k = 1:nmodes
    rpiz = rpiz + D(k)*exp(-1i*(m(k)*v+n(k)*u));
    drpizdv = drpizdv + D(k)*exp(-1i*(m(k)*v+n(k)*u)).*(-1i*m(k));
end
rpiz = exp(1i*v).*rpiz;
drpizdv = exp(1i*v).*(1i*rpiz + drpizdv);
drdv = real(drpizdv);
z = imag(drpizdv);
x = drdv.*cos(u);
y = drdv.*sin(u);
 
end