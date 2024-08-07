function integral = intacyc(f, n, nu, nv)
%INTACYC integrate a function along an A-cycle
%   Detailed explanation goes here

dom = f.domain;

% get A-cycle points
xonapatches = cell(nv,3);
xvonapatches = cell(nv,3);
for i = 1:nv
    % xvonapatches{i,1} = dom.xv{(i-1)*nu+1};
    xvonapatches{i,1} = dom.xv{i};
    xvonapatches{i,2} = dom.yv{i};
    xvonapatches{i,3} = dom.zv{i}; 
    xonapatches{i,1} = dom.x{i};
    xonapatches{i,2} = dom.y{i};
    xonapatches{i,3} = dom.z{i};
end

xv = zeros(n*nv,3); % xv on A-cycle
x = zeros(n*nv,3); % x on A-cycle
for i = 1:nv
    xv((i-1)*n+1:i*n,:) = [xvonapatches{i,1}(:,1).'; 
        xvonapatches{i,2}(:,1).'; 
        xvonapatches{i,3}(:,1).'].';
    x((i-1)*n+1:i*n,:) = [xonapatches{i,1}(:,1).';
        xonapatches{i,2}(:,1).'; 
        xonapatches{i,3}(:,1).'].';
end

% quiver3(x(:,1),x(:,2),x(:,3),xv(:,1),xv(:,2),xv(:,3))
% hold on
% plot(dom)

% normalize each row of xv
% for i = 1:nv
%     xvnorm = norm(xv(i,:));
%     xv(i,:) = xv(i,:)./xvnorm;
% end

% get A-cycle quad. weights
awts = zeros(n*nv,1);
for i = 1:nv
    [~, awts((i-1)*n+1:i*n)] = chebpts(n);
end
%awts = awts./nv;
%awts = nv.*awts./pi;

valsonapatches = cell(nv,3);
for i = 1:nv
    for j = 1:3
        valsonapatches{i,j} = f.components{j}.vals{i};
    end
end
avals = zeros(n*nv,3); % function values on A-cycle
for i = 1:nv
    for j = 1:3
        avals((i-1)*n+1:i*n,j) = valsonapatches{i,j}(:,1).';
    end
end

axv = dot(conj(avals),xv,2);
integral = dot(awts,axv);

end