lerr = load('lerr_saved/ellipse_k0p1_140824.mat').lerr;
lerr2 = load('lerr_saved/ellipse_k0p1_140824_2.mat').lerr;
lerr3 = load('lerr_saved/ellipse_k0p1_140824_3.mat').lerr;
lerr = [lerr lerr2 lerr3];

% npts = n^2*nu*nv = 3*n^2*nu^2
% h = sqrt(npts)

numns = 3;
numnus = 6;

lerrcol = 4;
loglog(sqrt(3)*lerr(1,1:numnus).*lerr(2,1:numnus), lerr(lerrcol,1:numnus), 'o-')
hold on
xlabel('h^{-1}')
ylabel('||B - B_0||_\infty/||B_0||_\infty')
for i = 1:numns-1
    cols = numnus*i+1:numnus*(i+1);
    loglog(sqrt(3)*lerr(1,cols).*lerr(2,cols), lerr(lerrcol,cols), 'o-')
end

% % BIEST results, lambda = 0.5, one surface
% N05 = [8.8e3 2.5e4 4.8e4 7.9e4]; % 9.8e4 2.2e5];
% einf05 = [2.6e-2 2.2e-3 3.8e-4 4.9e-5]; % 2.2e-5 5.3e-7];
% plot(sqrt(N05), einf05, '^:')
% 
% % lambda = 0, two surfaces
% N0 = [8.8e3 2.5e4 3.5e4 6.3e4]; % 7.9e4 1.9e5];
% einf0 = [2.3e-1 3.1e-2 6.6e-3 1.6e-3]; % 5.9e-4 9.4e-6];
% plot(sqrt(N0), einf0, '^:')

h = 40:60;
plot(h, h(1)^2*3e-2*h.^(-2), 'b--')

h = 50:80;
plot(h, h(1)^4*5e-2*h.^(-4), 'r--')

% legend('FMM+surfacefun, poly. order = p = 5', ...
%     'FMM+surfacefun, p = 7', ...
%     'O(h^{-4})', 'O(h^{-6})', 'Location', 'southwest')

lerrcol = 6;
figure(2)
loglog(sqrt(lerr(1,1:numnus).^2.*lerr(2,1:numnus).*lerr(3,1:numnus)), ...
    lerr(lerrcol,1:numnus), 'o-')
xlabel('1/h')
ylabel('||B - B_0||_\infty/||B_0||_\infty on surface')
hold on
for i = 1:numns-1
    cols = numnus*i+1:numnus*(i+1);
    loglog(sqrt(3)*lerr(1,cols).*lerr(2,cols), lerr(lerrcol,cols), 'o-')
end

h = 40:80;
plot(h, h(1)^2*8e-2*h.^(-2), 'b--')

h = 50:100;
plot(h, h(1)^3*2e-2*h.^(-3), 'r--')

h = 60:120;
plot(h, h(1)^4*1e-2*h.^(-4), 'k--')

legend('FMM+surfacefun, poly. order = p = 4', ...
    'FMM+surfacefun, p = 5', 'FMM+surfacefun, p = 6', ...
    'O(h^{-2})', 'O(h^{-3})', 'O(h^{-4})', 'Location', 'southwest')
