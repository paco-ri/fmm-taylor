% ===== copied from lap3d matlab test ==============
nsph = 7;
nref = 3;
dom = surfacemesh.sphere(nsph,nref);
S = surfer.surfacemesh_to_surfer(dom);
[srcvals,~,~,~,~,wts] = extract_arrays(S);
vn = normal(dom);

ndeg = 4; % 1 -> 6, 2 -> 10, 3 -> 14, n -> 2(2n+1)
mdeg = 2;

f = spherefun.sphharm(ndeg,mdeg);
sigma = surfacefun(@(x,y,z) f(x,y,z),dom);
sigvals = surfacefun_to_array(sigma,dom,S);

eps = 1e-7;

% ===================================================

if true
    fprintf('nsph = %d, nref = %d\n',nsph,nref)
    % without precomputed quadrature
    Bsigma = mtxBsigma(S,dom,sigma,0,eps);
    Bsigma_ex = -(ndeg+1)/(2*ndeg+1).*sigvals + sigvals;
    Bsigma = surfacefun_to_array(Bsigma,dom,S);
    err = sqrt(sum(abs(Bsigma_ex - Bsigma).^2.*wts));
    fprintf('error w/o Q = %f\n',err)
end
   
if false
    % with precomputed quadrature
    Q = taylor.static.get_quadrature_correction(S,eps,S);
    Bsigma = mtxBsigma(S,dom,sigma,eps,S,Q);
    Bsigma = surfacefun_to_array(Bsigma,dom,S);
    err = sqrt(sum(abs(Bsigma_ex - Bsigma).^2.*wts));
    fprintf('error w/ Q = %f\n',err)
end

if false
    Q = taylor.static.get_quadrature_correction(S,eps,S);
    S0sigr = taylor.static.eval_gradS0(S,real(sigvals),eps,S,Q);
    S0sigi = taylor.static.eval_gradS0(S,imag(sigvals),eps,S,Q);
    S0sigvals = S0sigr; % + 1i.*S0sigi;
    S0sig = array_to_surfacefun(S0sigvals.',dom,S);
    nS0sig = dot(vn,S0sig);
    % figure(1)
    % plot(sigma)
    % colorbar
    % figure(2)
    % plot(nS0sig)
    % colorbar
    nS0sigvals = surfacefun_to_array(nS0sig,dom,S);
    figure(3)
    plot(real(sigvals./nS0sigvals), 'o')
    hold on
    plot(imag(sigvals./nS0sigvals), 'o')
end

if false
    % interior points
    ntarg = 16;
    targs = zeros([3 ntarg]);
    targs(1,:) = (1 - 1./2.^(1:ntarg)).*cos(pi/3);
    targs(3,:) = (1 - 1./2.^(1:ntarg)).*sin(pi/3);
    targs(:,ntarg) = [cos(pi/3); 0; sin(pi/3)];
    targinfo = [];
    targinfo.r = targs;
    S0sigma = taylor.static.eval_gradS0(S,sigvals,eps,targinfo);
    %S0sigma = lap3d.dirichlet.eval(S,[0 1.0],rhs,eps,targinfo);
    nsigma = f(cos(pi/3),0,sin(pi/3)).*[cos(pi/3); 0; sin(pi/3)];
    disp(S0sigma)
    disp(nsigma)
    %disp(S0sigma(:,end) + S0sigma(:,end-1))
    % inside - nsigma = outside
end

if false
    % smaller shell
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

    % verify jump condition
    figure(1)
    plot(sigma-nS0sigmaint+nS0sigmaext)
    colorbar

    figure(2)
    plot(nS0sigma)
    colorbar

    % (int + ext)/2 = (ndeg)/(ndeg-1) bndry?
    figure(3)
    plot((nS0sigmaint+nS0sigmaext)./2)
    colorbar

    arr1 = surfacefun_to_array(nS0sigmaint,dom,S);
    arr2 = surfacefun_to_array(nS0sigmaext,dom,S);
    arr3 = surfacefun_to_array(nS0sigma,dom,S);
    arr4 = (arr1 + arr2)./arr3;

    % conclusion: eval_gradS0 + sigma/2 is the interior limit
end
