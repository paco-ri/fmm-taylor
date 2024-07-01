numnus = 5;
numns = 3;
loglog(sqrt(3)*lerr(1,1:numnus).*lerr(2,1:numnus), ...
    lerr(4,1:numnus), 'o-')
hold on
for i = 1:numns-1
    cols = numnus*i+1:numnus*(i+1);
    loglog(sqrt(3)*lerr(1,cols).*lerr(2,cols), ...
        lerr(4,cols), 'o-')
end

% BIEST results, lambda = 0.5, one surface
N05 = [8.8e3 2.5e4 4.8e4 7.9e4]; % 9.8e4 2.2e5];
einf05 = [2.6e-2 2.2e-3 3.8e-4 4.9e-5]; % 2.2e-5 5.3e-7];
plot(sqrt(N05), einf05, '^:')

% lambda = 0, two surfaces
N0 = [8.8e3 2.5e4 3.5e4 6.3e4]; % 7.9e4 1.9e5];
einf0 = [2.3e-1 3.1e-2 6.6e-3 1.6e-3]; % 5.9e-4 9.4e-6];
plot(sqrt(N0), einf0, '^:')