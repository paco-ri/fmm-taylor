function B0 = reftaylor(ntheta,rmin,rmaj,jmag,lambda,targ)
%REFTAYLOR compute a reference Taylor state corresponding to a current ring
% docs later 

source = zeros([ntheta 3]);
current = zeros([ntheta 3]);
for i = 1:ntheta
    theta = 2*pi*(i-1)/ntheta;
    source(i,1) = rmaj + rmin*cos(theta);
    source(i,3) = rmin*sin(theta);
    current(i,1) = -jmag*rmin*sin(theta);
    current(i,3) = jmag*rmin*cos(theta);
end

helmker = zeros([ntheta 1]);
gradhelmker = zeros([ntheta 3]);
for i = 1:ntheta
    for j = 1:3
        gradhelmker(i,j) = targ(j) - source(i,j);
    end
    dist = norm(gradhelmker(i,:));
    helmker(i) = exp(1i*lambda*dist)/dist;
    gradhelmker(i,:) = gradhelmker(i,:).*exp(1i*lambda*dist)...
        .*(1i*lambda - 1/dist)./dist^2;
end

% integrand = cross(current,gradhelmker);
integrand = cross(current,gradhelmker);
B0 = (sum(helmker.*current,1) - sum(integrand,1)).*rmin./(2*ntheta);

end