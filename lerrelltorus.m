% load('ellipse250724.mat')

% npts = n^2*nu*nv = 3*n^2*nu^2
% h = sqrt(npts) = sqrt(3)*n*nu

numnus = 4;
numns = 3;
lerrcol = 4;
loglog(sqrt(3)*lerr(1,1:numnus).*lerr(2,1:numnus), lerr(lerrcol,1:numnus), 'o-')
hold on
xlabel('h^{-1}')
ylabel('||B - B_0||_\infty/||B_0||_\infty')
for i = 1:numns-1
    cols = numnus*i+1:numnus*(i+1);
    loglog(sqrt(3)*lerr(1,cols).*lerr(2,cols), lerr(lerrcol,cols), 'o-')
end
        
h = 250:350;
plot(h, h(1)^3*1e-3*h.^(-3), 'b--')

h = 350:500;
plot(h, h(1)^5*1e-4*h.^(-5), 'r--')

h = 450:600;
plot(h, h(1)^7*1e-5*h.^(-7), 'k--')

legend('FMM+surfacefun, p = 4', ...
    'FMM+surfacefun, p = 6', ...
    'FMM+surfacefun, p = 8', ...
    'O(h^{-3})', 'O(h^{-5})', 'O(h^{-7})', 'Location', 'southwest')

lerrcol = 10;
figure(2)
loglog(sqrt(3)*lerr(1,1:numnus).*lerr(2,1:numnus), lerr(lerrcol,1:numnus), 'o-')
xlabel('1/h')
ylabel('||B - B_0||_\infty/||B_0||_\infty on surface')
hold on
for i = 1:numns-1
    cols = numnus*i+1:numnus*(i+1);
    loglog(sqrt(3)*lerr(1,cols).*lerr(2,cols), lerr(lerrcol,cols), 'o-')
end

h = 80:120;
plot(h, h(1)^4*8e-4*h.^(-4), 'b--')

h = 120:220;
plot(h, h(1)^6*1e-5*h.^(-6), 'r--')

h = 160:260;
plot(h, h(1)^8*8e-8*h.^(-8), 'k--')

legend('FMM+surfacefun, poly. order = p = 5', ...
    'FMM+surfacefun, p = 7', ...
    'FMM+surfacefun, p = 9', ...
    'O(h^{-4})', 'O(h^{-6})', 'O(h^{-8})', 'Location', 'southwest')
