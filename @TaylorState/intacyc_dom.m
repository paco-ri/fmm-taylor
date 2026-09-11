function [integrala, quadpts, quadwts] = intacyc_dom(domain, k, f)
%INTACYC_DOM Integrate along the A-cycle of surface k of a Domain
%   Chooses between the uniform and adaptive forms of INTACYC based
%   on whether surface k of DOMAIN carries a nontrivial quadforest
%
%   Arguments:
%     domain: Domain object
%     k: index of the surface whose A-cycle is integrated
%     f: surfacefunv living on domain.dom{k}

n  = domain.domparams(1);
nu = domain.domparams(2);
nv = domain.domparams(3);

if domain.is_adaptive(k)
    [integrala, quadpts, quadwts] = TaylorState.intacyc(f, n, nu, nv, ...
        domain.qf{k}, domain.p2q{k});
else
    [integrala, quadpts, quadwts] = TaylorState.intacyc(f, n, nu, nv);
end

end
