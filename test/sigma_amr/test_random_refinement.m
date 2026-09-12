% refine_leaves on random and clustered markings, checked by a
% Laplace-Beltrami solve on every refined mesh. Solves are carried
% out to check that quadtrees are balanced.
%
% Each trial refines every patch once, then marks 10% of the leaves below
% rmax at each further level:
%   random     leaves chosen uniformly at random, scattered over the surface
%   clustered  leaves nearest a random seed leaf, so one region is refined
%              repeatedly, as under Dorfler marking

fprintf('=== test_random_refinement start ===\n');

rmax = 4; % max quadtree depth
n = 5; nv = 3; nu = 3*nv;
d = prepare_stellarator(n,nu,nv,16,40); dom0 = d{1}; % base mesh
np0 = length(dom0.x); % number of base patches
fprintf('stellarator, %d base patches\n', np0);

nfail = 0; ntrial = 0;
% mode 1: uniformly random marks; mode 2: spatially clustered marks
for mode = 1:2
    if mode == 1, mname = 'random'; else, mname = 'clustered'; end
    for trial = 1:10
        rng(trial);
        ntrial = ntrial + 1;
        % start each trial from the uniform base mesh
        p2q = [(1:np0).', zeros(np0,2)];
        % refine one level per pass, feeding each mesh into the next pass
        for lev = 1:rmax
            npc = size(p2q,1);
            elig = find(p2q(:,2) < rmax); % leaves that can still be split
            if isempty(elig), break, end
            if lev == 1
                marked = (1:npc).';
            elseif mode == 1
                % 10% of eligible leaves at random
                k = max(1, round(0.10*numel(elig)));
                marked = elig(randperm(numel(elig), k));
            else
                % clustered: k patches nearest a random seed centroid
                k = max(1, round(0.10*numel(elig)));
                cx = cellfun(@(a) mean(a(:)), domR.x);
                cy = cellfun(@(a) mean(a(:)), domR.y);
                cz = cellfun(@(a) mean(a(:)), domR.z);
                s = elig(randi(numel(elig)));
                dd = (cx(elig)-cx(s)).^2 + (cy(elig)-cy(s)).^2 + (cz(elig)-cz(s)).^2;
                [~, ord] = sort(dd);
                marked = elig(ord(1:k));
            end
            [domR, qfR, p2qR] = surfacemesh.refine_leaves(dom0, p2q, marked, rmax);
            try
                % hodge() does two Laplace-Beltrami solves on the tangential
                % field n x phihat, where phihat is the toroidal unit vector
                vn = normal(domR);
                phihat = surfacefunv(@(x,y,z) -y./sqrt(x.^2+y.^2), ...
                                     @(x,y,z)  x./sqrt(x.^2+y.^2), ...
                                     @(x,y,z) 0.*z, domR);
                dummy = cross(vn, phihat);
                [~, ~, vH] = hodge(dummy);   %#ok<ASGLU>
            catch ME
                fprintf('FAIL mode=%s trial=%d level=%d npat=%d : %s\n', ...
                    mname, trial, lev, length(domR.x), ME.message);
                nfail = nfail + 1;
                break
            end
            p2q = p2qR;
        end
        fprintf('mode=%s trial %2d done (%d patches)\n', mname, trial, size(p2q,1));
    end
end

fprintf('\nfailed trials: %d of %d\n', nfail, ntrial);
if nfail == 0
    fprintf('VERDICT: PASS\n');
else
    fprintf('VERDICT: FAIL\n');
end
fprintf('DONE_SENTINEL\n');
