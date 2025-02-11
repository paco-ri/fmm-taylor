ns = [9 11];
nvs = [6 8 10];

lerr = zeros(4,size(ns,2)*size(nvs,2));
lind = 1;

for n = ns
    for nv = nvs
        nu = nv*3;
        dom = prepare_torus(n,nu,nv,16,40);
        pdo = [];
        pdo.lap = 1;

        cosphi = @(x,y,z) x./sqrt(x.^2 + y.^2);
        c = surfacefun(@(x,y,z) cosphi(x,y,z)+z, dom{1});

        L = surfaceop(dom{1},pdo,c);
        L.rankdef = true;
        u = L.solve();
    end
end

plot(u)