function [integral] = intacyc(f, n, nu, nv)
%INTACYC integrate a function along an A-cycle
%   Detailed explanation goes here

dom = f.domain;

% get A-cycle points
xuonapatches = cell(nu,3);
for i = 1:nu
    xuonapatches{i,1} = dom.xu{(i-1)*nv+1};
    xuonapatches{i,2} = dom.yu{(i-1)*nv+1};
    xuonapatches{i,3} = dom.zu{(i-1)*nv+1}; 
end

xu = zeros(n*nu,3); % xu on A-cycle
for i = 1:nu
    xu((i-1)*n+1:i*n,:) = [xuonapatches{i,1}(1,:); 
        xuonapatches{i,2}(1,:); 
        xuonapatches{i,3}(1,:)].';
end

% get A-cycle quad. weights
awts = zeros(n*nu,1);
for i = 1:nu
    [~, awts((i-1)*n+1:i*n)] = chebpts(n);
end
%awts = awts./nu;
%awts = nu.*awts./pi;

valsonapatches = cell(nu,3);
for i = 1:nu
    for j = 1:3
        valsonapatches{i,j} = f.components{j}.vals{(i-1)*nv+1};
    end
end
avals = zeros(n*nu,3); % function values on A-cycle
for i = 1:nu
    for j = 1:3
        avals((i-1)*n+1:i*n,j) = valsonapatches{i,j}(1,:);
    end
end

axu = dot(conj(avals),xu,2);
integral = dot(awts,axu);

end