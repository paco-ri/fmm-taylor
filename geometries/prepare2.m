function [doms, fluxnodes, fluxwts] = prepare2(ns, nus, nvs, nr, nt, np, modes)

nsurf = 2;
doms = cell(1,nsurf);
for i = 1:nsurf
    n = ns(i);
    nu = nus(i);
    nv = nvs(i);
    
    x = cell(nu*nv, 1);
    y = cell(nu*nv, 1);
    z = cell(nu*nv, 1);
    ubreaks = linspace(0, 2*pi, nu+1);
    vbreaks = linspace(0, 2*pi, nv+1);
    
    k = 1;
    for ku = 1:nu
        for kv = 1:nv
            [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
            [x{k}, y{k}, z{k}] = eval_geom(uu, vv, modes{i});
            k = k+1;
        end
    end
    
    % Try to reorder to create better merge indices:
    if ( bitand(nv, nv-1) == 0 && bitand(nu, nu-1) == 0 )
        ordering = morton(nv, nu);
        ordering = ordering(:);
        x(ordering) = x;
        y(ordering) = y;
        z(ordering) = z;
    end
    
    doms{i} = surfacemesh(x, y, z);

end

fluxnodes = cell(1,nsurf);
fluxwts = cell(1,nsurf);

[rnodes, rwts] = chebpts(nr,[0 1],1);
fluxnodes{1} = zeros(3,nt*nr);
fluxwts{1} = zeros(1,nt*nr);
for i = 1:nr
    rr = rnodes(i);
    wr = rwts(i);
    for j = 1:nt
        ij = (i-1)*nt+j;
        tt = 2*pi*(j-1)/nt;
        [go1, ~, go2] = eval_geom(0,tt,modes{1});
        [gi1, ~, gi2] = eval_geom(0,tt,modes{2});
        [dgo1, ~, dgo2] = dv_eval_geom(0,tt,modes{1});
        [dgi1, ~, dgi2] = dv_eval_geom(0,tt,modes{2});
        fluxnodes{1}(:,ij) = (1-rr)*[gi1; 0; gi2] + rr*[go1; 0; go2];
        fluxwts{1}(1,ij) = (2*pi/nt) ...
            * wr*((-gi1+go1)*((1-rr)*dgi2+rr*dgo2) ...
            - (-gi2+go2)*((1-rr)*dgi1+rr*dgo1));
    end
end

fluxnodes{2} = zeros([3 nr*np]);
fluxwts{2} = zeros([1 nr*np]);
for i = 1:nr
    rr = rnodes(i);
    wr = rwts(i);
    for j = 1:np
        ij = (i-1)*np+j;
        pp = 2*pi*(j-1)/np;
        [go1, go2] = eval_geom(pp,0,modes{1});
        [gi1, gi2] = eval_geom(pp,0,modes{2});
        [dgo1, dgo2] = du_eval_geom(pp,0,modes{1});
        [dgi1, dgi2] = du_eval_geom(pp,0,modes{2});
        fluxnodes{2}(:,ij) = (1-rr)*[gi1; gi2; 0] + rr*[go1; go2; 0];
        fluxwts{2}(1,ij) = (2*pi/np) ...
            * wr*((-gi1+go1)*((1-rr)*dgi2+rr*dgo2) ...
            - (-gi2+go2)*((1-rr)*dgi1+rr*dgo1));
    end
end

end