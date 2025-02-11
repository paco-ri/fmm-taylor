function B0 = reftaylorsurffun(dom,n,nu,nv,ntheta,rmin,rmaj,jmag,lambda)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

B0 = surfacefunv(dom);
B0x = cell(nu*nv,1);
B0y = cell(nu*nv,1);
B0z = cell(nu*nv,1);
for i = 1:nu*nv
    B0x{i} = zeros(n);
    B0y{i} = zeros(n);
    B0z{i} = zeros(n);
    for j = 1:n
        for k = 1:n
            % B0eval = reftaylor(ntheta,rmin,rmaj,jmag,lambda,...
            %     [dom.x{i}(j,k); dom.y{i}(j,k); dom.z{i}(j,k)]);
            B0eval = reftaylor(ntheta,rmin,rmaj,jmag,lambda,...
                [dom.x{i}(k,j); dom.y{i}(k,j); dom.z{i}(k,j)]);
            % B0x{i}(j,k) = B0eval(1);
            % B0y{i}(j,k) = B0eval(2);
            % B0z{i}(j,k) = B0eval(3);
            B0x{i}(k,j) = B0eval(1);
            B0y{i}(k,j) = B0eval(2);
            B0z{i}(k,j) = B0eval(3);
        end
    end
end
B0.components{1} = surfacefun(B0x,dom);
B0.components{2} = surfacefun(B0y,dom);
B0.components{3} = surfacefun(B0z,dom);

end