% test \nabla_\Gamma \cdot m = i \lambda \sigma

m0 = debyem0(sigma,lambda);
m0err = div(m0) - 1i*lambda.*sigma;
norm(real(m0err))
norm(imag(m0err))
norm(m0err)

% nsph = 7;
% nref = 3;
% dom = surfacemesh.sphere(nsph,nref);
% S = surfer.surfacemesh_to_surfer(dom);
% [srcvals,~,~,~,~,wts] = extract_arrays(S);
% vn = normal(dom);
% 
% ndeg = 3; 
% mdeg = 1;
% 
% f = spherefun.sphharm(ndeg,mdeg);
% sigma = surfacefun(@(x,y,z) f(x,y,z),dom);
% sigvals = surfacefun_to_array(sigma,dom,S);
% 
% lambda = 1.0;
% m0 = debyem0(sigma,lambda);
% norm(div(m0)-1i*lambda.*sigma)
