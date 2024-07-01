function B0 = reftaylor(ntheta,rmin,rmaj,jmag,lambda,targ)
%REFTAYLOR compute a reference Taylor state corresponding to a current ring
% docs later 

source = zeros([3 ntheta]); % points on a ring in the xz-plane 
current = zeros([3 ntheta]); % current at these points
for i = 1:ntheta
    theta = 2*pi*(i-1)/ntheta;
    source(1,i) = rmaj + rmin*cos(theta);
    % source(2,i) = rmaj + rmin*cos(theta);
    source(3,i) = rmin*sin(theta);
    current(1,i) = -jmag*rmin*sin(theta);
    % current(2,i) = -jmag*rmin*sin(theta);   
    current(3,i) = jmag*rmin*cos(theta);
end

if lambda == 0
    lapker = zeros([1 ntheta]);
    gradlapker = zeros([3 ntheta]);
    for i = 1:ntheta
        rr = targ - source(:,i);
        lapker(i) = 1/norm(rr);
        gradlapker(:,i) = rr./norm(rr)^3;
    end

    integrand = cross(gradlapker,current);
    B0 = sum(integrand,2).*rmin./(2*ntheta);
else
    helmker = zeros([1 ntheta]);
    gradhelmker = zeros([3 theta]);
    for i = 1:ntheta
        for j = 1:3
            gradhelmker(j,i) = targ(j) - source(j,i);
        end
        dist = norm(gradhelmker(:,i));
        helmker(i) = exp(1i*lambda*dist)/dist;
        gradhelmker(:,i) = gradhelmker(:,i).*exp(1i*lambda*dist)...
            .*(1i*lambda - 1/dist)./dist^2;
    end
    
    integrand = cross(gradhelmker,current);
    B0 = (lambda*sum(helmker.*current,1) + sum(integrand,2)).*rmin./(2*ntheta);
end

end