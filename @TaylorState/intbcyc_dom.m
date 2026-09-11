function [integralb, quadpts, quadwts] = intbcyc_dom(domain, k, f)
%INTBCYC_DOM Integrate along the B-cycle of surface k of a Domain
%   Chooses between the uniform and adaptive forms of INTBCYC based
%   on whether surface k of DOMAIN carries a nontrivial quadforest
%
%   Arguments:
%     domain: Domain object
%     k: index of the surface whose B-cycle is integrated
%     f: surfacefunv living on domain.dom{k}

n  = domain.domparams(1);
nu = domain.domparams(2);
nv = domain.domparams(3);

if domain.is_adaptive(k)
    [integralb, quadpts, quadwts] = TaylorState.intbcyc(f, n, nu, nv, ...
        domain.qf{k}, domain.p2q{k});
else
    [integralb, quadpts, quadwts] = TaylorState.intbcyc(f, n, nu, nv);
end

end
