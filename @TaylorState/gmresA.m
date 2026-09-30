function A = gmresA(s,domain,zk,epstaylor,epslh,opts,optslh)

% cell arrays: S, dom, opts, optslh
S = domain.surf;
Bsigma = TaylorState.mtxBsigma(domain,s,zk,epstaylor, ...
    epslh,S,opts,optslh);

dom = domain.dom;
if length(dom) == 1
    A = surfacefun_to_array(Bsigma{1,1},dom{1},S{1});
else
    Ao = surfacefun_to_array(Bsigma{1},dom{1},S{1});
    Ai = surfacefun_to_array(Bsigma{2},dom{2},S{2});
    A = [Ao; Ai];
end

if abs(zk) < eps
    % At zk = 0 the equilibrium density is a null vector of A11 and produces
    % no interior field; the rank-one term removes it, fixing int sigma ~ 0
    w = cellfun(@(Si) Si.wts(:), S, 'UniformOutput', false);
    w = vertcat(w{:});
    A = A + (w.'*s)/sum(w);
end

end