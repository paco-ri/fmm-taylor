function integralb = intbcyc(f, n, nu)
%INTACYC integrate a function along an B-cycle
%   Detailed explanation goes here

dom = f.domain;
nv = length(dom.x)/nu;

% get A-cycle points
xonbpatches = cell(nu,3);
xuonbpatches = cell(nu,3);
for i = 1:nu
    % xvonapatches{i,1} = dom.xv{(i-1)*nu+1};
    idx = (i-1)*nv+1;
    xuonbpatches{i,1} = dom.xu{idx};
    xuonbpatches{i,2} = dom.yu{idx};
    xuonbpatches{i,3} = dom.zu{idx}; 
    xonbpatches{i,1} = dom.x{idx};
    xonbpatches{i,2} = dom.y{idx};
    xonbpatches{i,3} = dom.z{idx};
end

xu = zeros(n*nu,3); % xu on B-cycle
x = zeros(n*nu,3); % x on A-cycle
for i = 1:nu
    xu((i-1)*n+1:i*n,:) = [xuonbpatches{i,1}(1,:); 
        xuonbpatches{i,2}(1,:); 
        xuonbpatches{i,3}(1,:)].';
    x((i-1)*n+1:i*n,:) = [xonbpatches{i,1}(1,:);
        xonbpatches{i,2}(1,:); 
        xonbpatches{i,3}(1,:)].';
end
% 
% quiver3(x(:,1),x(:,2),x(:,3),xu(:,1),xu(:,2),xu(:,3),.6)
% hold on
% plot(dom)

% get B-cycle quad. weights
bwts = zeros(n*nu,1);
for i = 1:nu
    [~, bwts((i-1)*n+1:i*n)] = chebpts(n);
end

valsonbpatches = cell(nu,3);
for i = 1:nu
    for j = 1:3
        idx = (i-1)*nv+1;
        valsonbpatches{i,j} = f.components{j}.vals{idx};
    end
end
bvals = zeros(n*nu,3); % function values on B-cycle
for i = 1:nu
    for j = 1:3
        bvals((i-1)*n+1:i*n,j) = valsonbpatches{i,j}(1,:);
    end
end

bxu = dot(conj(bvals),xu,2);
integralb = dot(bwts,bxu);

end