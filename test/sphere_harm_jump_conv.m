ns = 4:6;
nrefs = 2:4;

% spherical harmonic indices
ndeg = 4;
mdeg = 0;

lerr_sph = zeros([size(ns,1) size(nrefs,1) 8]);
nind = 1;
for n = ns
    nrefind = 1;
    for nref = nrefs
        fprintf('n = %d, nref = %d\n', n, nref)
        tic
        dom = surfacemesh.sphere(n,nref);
        S = surfer.surfacemesh_to_surfer(dom);
        [srcvals,~,~,~,~,wts] = extract_arrays(S);
        vn = normal(dom);

        f = spherefun.sphharm(ndeg,mdeg);
        sigma = surfacefun(@(x,y,z) f(x,y,z),dom);
        sigvals = surfacefun_to_array(sigma,dom,S);

        eps = 1e-7;

        targs = S.r;
        targinfoint = [];
        targinfoint.r = 0.99.*targs;
        targinfoext = [];
        targinfoext.r = 1.01.*targs;
    
        S0sigmaintvals = taylor.static.eval_gradS0(S,sigvals,eps,targinfoint);
        S0sigmaint = array_to_surfacefun(S0sigmaintvals.',dom,S);
        nS0sigmaint = dot(vn,S0sigmaint);
    
        S0sigmavals = taylor.static.eval_gradS0(S,sigvals,eps,S);
        S0sigma = array_to_surfacefun(S0sigmavals.',dom,S);
        nS0sigma = dot(vn,S0sigma);
    
        S0sigmaextvals = taylor.static.eval_gradS0(S,sigvals,eps,targinfoext);
        S0sigmaext = array_to_surfacefun(S0sigmaextvals.',dom,S);
        nS0sigmaext = dot(vn,S0sigmaext);

        res = sigma - nS0sigmaint + nS0sigmaext;
        
        lerr_sph(nind,nrefind,1) = n;
        lerr_sph(nind,nrefind,2) = nref;
        lerr_sph(nind,nrefind,3) = norm(res,1);
        lerr_sph(nind,nrefind,4) = norm(res,1)/norm(sigma,1);
        lerr_sph(nind,nrefind,5) = norm(res,2);
        lerr_sph(nind,nrefind,6) = norm(res,2)/norm(sigma,2);
        lerr_sph(nind,nrefind,7) = norm(res,inf);
        lerr_sph(nind,nrefind,8) = norm(res,inf)/norm(sigma,inf);
        toc
        nrefind = nrefind+1;
    end
    nind = nind+1;
end