function [integrala, quadpts, quadwts] = intacyc(f, n, ~, nv, varargin)
%INTACYC integrate a function along an A-cycle
%   Detailed explanation goes here

if (nargin > 5)
    qf = varargin{1};
    p2q = varargin{2};
    adap_dom = true;
else
    adap_dom = false;
end

dom = f.domain;

% Initialize arrays for coordinates and values on the A-cycle
xv = zeros(n*nv, 3); % xv on A-cycle (derivatives or field)
x = zeros(n*nv, 3);  % x on A-cycle (coordinates)
avals = zeros(n*nv, 3); % function values on A-cycle

% Get A-cycle quadrature weights
awts = zeros(n*nv, 1);

acyc_p2q = nan(length(dom), 3);
if adap_dom
    k = 1;
    for i = 1:nv
        idx = i;

        % Get Morton codes for the current patch
        morton = qf.morton{idx};
        if ~ismember(idx, qf.tree_roots) % If the patch has no Morton codes, mark it for A-cycle integration
            acyc_p2q(k, :) = [idx, 0, 0]; % Mark patches with no Morton codes
            k = k + 1;
        end
        % Loop over Morton codes to find patches that have points on the A-cycle
        for level = 1:length(morton)
            for code = morton{level}
                [x_code, ~] = qf.deinterleave(code, level);
                if x_code == 2^(level-1)
                    acyc_p2q(k, :) = [idx, level, code];
                    k = k + 1;
                end
            end
        end
    end
    acyc_p2q = rmmissing(acyc_p2q);

    for i = 1:length(dom)
        if ismember(p2q(i, :), acyc_p2q, 'rows')
            idx = i;
            populate_patch_data(i, idx);
        end
    end
else
    for i = 1:nv
        idx = i;
        populate_patch_data(i, idx);
    end
end

axv = dot(conj(avals),xv,2);
integrala = dot(awts,axv);
quadpts = x;
quadwts = awts;

    function populate_patch_data(i, idx)
        rows = (i-1)*n+1:i*n;

        % Extract the xv (vector field) and x (coordinates) for this patch
        xv(rows, :) = [dom.xv{idx}(:,1).'; dom.yv{idx}(:,1).'; dom.zv{idx}(:,1).'].';
        x(rows, :) = [dom.x{idx}(:,1).'; dom.y{idx}(:,1).'; dom.z{idx}(:,1).'].';

        % Extract the function values on the A-cycle for this patch
        avals(rows, 1) = f.components{1}.vals{idx}(:,1).';
        avals(rows, 2) = f.components{2}.vals{idx}(:,1).';
        avals(rows, 3) = f.components{3}.vals{idx}(:,1).';

        % Compute quadrature weights using Chebyshev points
        [~, awts(rows)] = chebpts(n);
    end

end