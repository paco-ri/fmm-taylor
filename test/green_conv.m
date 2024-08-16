ns = [4 5];%12;
nus = 4:4:16;%24;
% nus = [1 2 3 4]; % for sphere

greenerrs = zeros([4 size(ns,2)*size(nus,2)]);
zk = 1.0 + 0.0i;

ii = 1;
eps = 1e-6;
for n = ns
    for nu = nus
        nv = nu*3;
        rmin = 1.0;
        rmaj = 2.0;
        
        % dom = surfacemesh.sphere(n,nu);
        dom = circulartorus(n,nu,nv,rmin,rmaj);
        % dom = surfacemesh.torus(n,nu,nv);
        
        r = @(x,y,z) sqrt(x.^2+y.^2+z.^2); 
        u = surfacefun(@(x,y,z) (exp(1i.*zk.*r(x,y,z))./4./pi)./r(x,y,z), dom);
        du = surfacefunv(@(x,y,z) (x.*exp(1i.*zk.*r(x,y,z))./4./pi).*(1i*zk./r(x,y,z).^2 - 1./r(x,y,z).^3), ...
            @(x,y,z) (y.*exp(1i.*zk.*r(x,y,z))./4./pi).*(1i.*zk./r(x,y,z).^2 - 1./r(x,y,z).^3), ...
            @(x,y,z) (z.*exp(1i.*zk.*r(x,y,z))./4./pi).*(1i.*zk./r(x,y,z).^2 - 1./r(x,y,z).^3), dom);
        % u = surfacefun(@(x,y,z) (1./4./pi)./r(x,y,z), dom);
        % du = surfacefunv(@(x,y,z) (x./4./pi).*(-1./r(x,y,z).^3), ...
        %     @(x,y,z) (y./4./pi).*(-1./r(x,y,z).^3), ...
        %     @(x,y,z) (z./4./pi).*(-1./r(x,y,z).^3), dom);
        vn = normal(dom);
        dudn = dot(vn,du);

        % sinphi = @(x,y,z) y./sqrt(x.^2 + y.^2);
        % cosphi = @(x,y,z) x./sqrt(x.^2 + y.^2);
        % sintheta = @(x,y,z) z./rmin;
        % costheta = @(x,y,z) (sqrt(x.^2 + y.^2) - rmaj)./rmin;
        % vn2 = surfacefunv(@(x,y,z) costheta(x,y,z).*cosphi(x,y,z), ...
        %                  @(x,y,z) costheta(x,y,z).*sinphi(x,y,z), ...
        %                  @(x,y,z) sintheta(x,y,z), dom);
        % dudn2 = dot(vn2,du);
        
        opts = [];
        S = surfer.surfacemesh_to_surfer(dom,opts);

        dudnvals = surfacefun_to_array(dudn,dom,S);
        uvals = surfacefun_to_array(u,dom,S);

        if zk == 0
            S0vals = lap3d.dirichlet.eval(S,dudnvals,S,eps,[1.0 0]);
            D0vals = lap3d.dirichlet.eval(S,uvals,S,eps,[0 1.0]);
            S0 = array_to_surfacefun(S0vals,dom,S);
            D0 = array_to_surfacefun(D0vals,dom,S);
        else
            Skvals = helm3d.dirichlet.eval(S,dudnvals,S,eps,zk,[1.0 0]);
            Dkvals = helm3d.dirichlet.eval(S,uvals,S,eps,zk,[0 1.0]);
            Sk = array_to_surfacefun(Skvals,dom,S);
            Dk = array_to_surfacefun(Dkvals,dom,S);
        end

        greenerrs(1,ii) = n;
        greenerrs(2,ii) = nu;
        greenerrs(3,ii) = nv;
        if zk == 0
            greenerrs(4,ii) = norm(u./2-S0+D0,inf)/norm(u,inf);
        else
            greenerrs(4,ii) = norm(u./2-Sk+Dk,inf)/norm(u,inf);
        end

        ii = ii + 1;
        fprintf('done with n=%d, nu=%d, nv=%d\n',n,nu,nv)
    end
    eps = eps*1e-4;
end

figure(1)
loglog(greenerrs(1,1:4).*greenerrs(2,1:4), greenerrs(4,1:4), 'o-')
hold on
loglog(greenerrs(1,5:8).*greenerrs(2,5:8), greenerrs(4,5:8), 'o-')
% loglog(greenerrs(1,9:12).*greenerrs(2,9:12),greenerrs(4,9:12),'o-')
hs = 20:100;
d1 = plot(hs,hs.^(-3));
d2 = plot(hs,hs.^(-4));
% d3 = plot(hs,2*20^8*hs.^(-11));
xlabel('1/h')
ylabel('relative L^\infty error')
% legend('poly. order = n = 4', 'n = 8', 'n = 12', 'O(h^{-3})', ...
    % 'O(h^{-7})', 'O(h^{-11})', 'Location', 'southwest')
legend('poly. order = p = 4', 'p = 5', 'O(h^{-3})', ...
    'O(h^{-4})', 'Location', 'southwest')



