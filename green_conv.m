ns = 4:4:8;%12;
nus = 4:4:16;%24;

greenerrs = zeros([4 size(ns,2)*size(nus,2)]);

ii = 1;
eps = 1e-6;
for n = ns
    for nu = nus
        nv = nu*3;
        rmin = 1.0;
        rmaj = 2.0;
        
        dom = circulartorus(n,nu,nv,rmin,rmaj);
        % dom = surfacemesh.torus(n,nu,nv);
        
        u = surfacefun(@(x,y,z) (1/4/pi)./sqrt(x.^2+y.^2+z.^2), dom);
        du = surfacefunv(@(x,y,z) (-x/4/pi)./sqrt(x.^2+y.^2+z.^2).^3, ...
            @(x,y,z) (-y/4/pi)./sqrt(x.^2+y.^2+z.^2).^3, ...
            @(x,y,z) (-z/4/pi)./sqrt(x.^2+y.^2+z.^2).^3, dom);
        vn = normal(dom);
        dudn = dot(vn,du);

        sinphi = @(x,y,z) y./sqrt(x.^2 + y.^2);
        cosphi = @(x,y,z) x./sqrt(x.^2 + y.^2);
        sintheta = @(x,y,z) z./rmin;
        costheta = @(x,y,z) (sqrt(x.^2 + y.^2) - rmaj)./rmin;
        vn2 = surfacefunv(@(x,y,z) costheta(x,y,z).*cosphi(x,y,z), ...
                         @(x,y,z) costheta(x,y,z).*sinphi(x,y,z), ...
                         @(x,y,z) sintheta(x,y,z), dom);
        dudn2 = dot(vn2,du);
        
        opts = [];
        % opts.iptype = 12;
        S = surfer.surfacemesh_to_surfer(dom,opts);
        [srcvals,~,~,~,~,~] = extract_arrays(S);
        vn3 = array_to_surfacefun(srcvals(10:12,:).',dom,S);

        dudnvals = surfacefun_to_array(dudn,dom,S);
        uvals = surfacefun_to_array(u,dom,S);

        S0vals = lap3d.dirichlet.eval(S,[1.0 0],dudnvals,eps);
        D0vals = lap3d.dirichlet.eval(S,[0 1.0],uvals,eps);
        S0 = array_to_surfacefun(S0vals,dom,S);
        D0 = array_to_surfacefun(D0vals,dom,S);

        greenerrs(1,ii) = n;
        greenerrs(2,ii) = nu;
        greenerrs(3,ii) = nv;
        greenerrs(4,ii) = norm(u./2-S0+D0,inf)/norm(u,inf);

        ii = ii + 1;
        fprintf('done with n=%d, nu=%d, nv=%d\n',n,nu,nv)
    end
    eps = eps*1e-3;
end

figure(1)
loglog(greenerrs(1,1:4).*greenerrs(2,1:4), greenerrs(4,1:4), 'o-')
hold on
loglog(greenerrs(1,5:8).*greenerrs(2,5:8), greenerrs(4,5:8), 'o-')
loglog(greenerrs(1,9:12).*greenerrs(2,9:12),greenerrs(4,9:12),'o-')
hs = 20:100;
d1 = plot(hs,200*hs.^(-3));
d2 = plot(hs,20^5*hs.^(-7));
d3 = plot(hs,2*20^8*hs.^(-11));
xlabel('1/h')
ylabel('relative L^\infty error')
legend('poly. order = n = 4', 'n = 8', 'n = 12', 'O(h^{-3})', ...
    'O(h^{-7})', 'O(h^{-11})', 'Location', 'southwest')



