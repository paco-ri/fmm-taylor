function integralb = intbcyc(f, n, nu)
%INTACYC integrate a function along an B-cycle
%   Detailed explanation goes here

dom = f.domain;
nv = length(dom.x)/nu;

% Initialize arrays for coordinates and values on the B-cycle
xu = zeros(n*nu, 3); % xu on B-cycle (derivatives or field)
x = zeros(n*nu, 3);  % x on B-cycle (coordinates)
bvals = zeros(n*nu, 3); % function values on B-cycle

% Get B-cycle quadrature weights
bwts = zeros(n*nu, 1);

for i = 1:nu
    % Index of the first vertex in the current patch
    idx = (i-1)*nv + 1;
    
    % Extract the xu (vector field) and x (coordinates) for this patch
    xu((i-1)*n+1:i*n, :) = [dom.xu{idx}(1,:); dom.yu{idx}(1,:); dom.zu{idx}(1,:)].';
    x((i-1)*n+1:i*n, :) = [dom.x{idx}(1,:); dom.y{idx}(1,:); dom.z{idx}(1,:)].';
    
    % Extract the function values on the B-cycle for this patch
    bvals((i-1)*n+1:i*n, 1) = f.components{1}.vals{idx}(1,:);
    bvals((i-1)*n+1:i*n, 2) = f.components{2}.vals{idx}(1,:);
    bvals((i-1)*n+1:i*n, 3) = f.components{3}.vals{idx}(1,:);
    
    % Compute quadrature weights using Chebyshev points
    [~, bwts((i-1)*n+1:i*n)] = chebpts(n);
end

bxu = dot(conj(bvals),xu,2);
integralb = dot(bwts,bxu);

end