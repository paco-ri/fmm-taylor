function [integralb, quadpts, quadwts] = intbcyc(f, n, nu, nv, varargin)
%INTACYC integrate a function along an B-cycle
%   Detailed explanation goes here

if (nargin > 4)
    qf = varargin{1};
    p2q = varargin{2};
    adap_dom = true;
else
    adap_dom = false;
end

dom = f.domain;

% Initialize arrays for coordinates and values on the B-cycle
if adap_dom
    nalloc = n*length(dom);
else
    nalloc = n*nu;
end
xu = zeros(nalloc, 3); % xu on B-cycle (derivatives or field)
x = zeros(nalloc, 3);  % x on B-cycle (coordinates)
bvals = zeros(nalloc, 3); % function values on B-cycle

% Get B-cycle quadrature weights
bwts = zeros(nalloc, 1);

if adap_dom
    bcyc_p2q = nan(length(dom), 3);
    k = 1;
    for i = 1:nu
        % Index of the first vertex in the current patch
        idx = (i-1)*nv + 1;

        % Get Morton codes for the current patch
        morton = qf.morton{idx};
        if ~ismember(idx, qf.tree_roots) % If the patch has no Morton codes, mark it for B-cycle integration
            bcyc_p2q(k, :) = [idx, 0, 0]; % Mark patches with no Morton codes
            k = k + 1;
        end
        % Loop over Morton codes to find patches that have points on the B-cycle
        for level = 1:length(morton)
            for code = morton{level}
                [~, y_code] = qf.deinterleave(code, level);
                if y_code == 0 % 2^(level-1)
                    bcyc_p2q(k, :) = [idx, level, code];
                    k = k + 1;
                end
            end
        end
    end
    bcyc_p2q = rmmissing(bcyc_p2q);

    for i = 1:length(dom)
        if ismember(p2q(i, :), bcyc_p2q, 'rows')
            idx = i;
            populate_patch_data(i, idx);
        end
    end

    nrows = length(dom)*n;
    xu = xu(1:nrows, :);
    x = x(1:nrows, :);
    bvals = bvals(1:nrows, :);
    bwts = bwts(1:nrows);
else
    for i = 1:nu
        % Index of the first vertex in the current patch
        idx = (i-1)*nv + 1;
        populate_patch_data(i, idx);
    end
end

bxu = dot(conj(bvals),xu,2);
integralb = dot(bwts,bxu);
quadpts = x;
quadwts = bwts;

    function populate_patch_data(i, idx)
        rows = (i-1)*n+1:i*n;

        % Extract the xu (vector field) and x (coordinates) for this patch
        xu(rows, :) = [dom.xu{idx}(1,:); dom.yu{idx}(1,:); dom.zu{idx}(1,:)].';
        x(rows, :) = [dom.x{idx}(1,:); dom.y{idx}(1,:); dom.z{idx}(1,:)].';

        % Extract the function values on the B-cycle for this patch
        bvals(rows, 1) = f.components{1}.vals{idx}(1,:);
        bvals(rows, 2) = f.components{2}.vals{idx}(1,:);
        bvals(rows, 3) = f.components{3}.vals{idx}(1,:);

        % Compute quadrature weights using Chebyshev points
        [~, bwts(rows)] = chebpts(n);
    end

end