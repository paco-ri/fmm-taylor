load('lerrcirctorus280624.mat')

% npts = n^2*nu*nv = 3*n^2*nu^2
% h = sqrt(npts) = sqrt(3)*n*nu

loglog(sqrt(3)*lerr(1,1:4).*lerr(2,1:4), lerr(4,1:4), 'o-')
xlabel('h^{-1}')
ylabel('||B - B_0||_\infty/||B_0||_\infty')
hold on
for i = 1:2
    loglog(sqrt(3)*lerr(1,4*i+1:4*(i+1)).*lerr(2,4*i+1:4*(i+1)), ...
        lerr(4,4*i+1:4*(i+1)), 'o-')
end

% BIEST results, lambda = 0.5, one surface
N05 = [8.8e3 2.5e4 4.8e4 7.9e4]; % 9.8e4 2.2e5];
einf05 = [2.6e-2 2.2e-3 3.8e-4 4.9e-5]; % 2.2e-5 5.3e-7];
plot(sqrt(N05), einf05, '^:')

% lambda = 0, two surfaces
N0 = [8.8e3 2.5e4 3.5e4 6.3e4]; % 7.9e4 1.9e5];
einf0 = [2.3e-1 3.1e-2 6.6e-3 1.6e-3]; % 5.9e-4 9.4e-6];
plot(sqrt(N0), einf0, '^:')

h = 80:160;
plot(h, 80^3*3e-3*h.^(-3), 'b--')

h = 120:220;
plot(h, 120^5*5e-5*h.^(-5), 'r--')

h = 150:270;
plot(h, 150^7*3e-7*h.^(-7), 'k--')

legend('FMM+surfacefun, poly. order = n = 6', ...
    'FMM+surfacefun, n = 8', ...
    'FMM+surfacefun, n = 10', ...
    'BIEST, \lambda = 0.5, one non-axisymm. surface', ...
    'BIEST, \lambda = 0, two non-axisymm. surfaces', ...
    'O(h^{-3})', 'O(h^{-5})', 'O(h^{-7})', 'Location', 'southwest')

figure(2)
loglog(sqrt(3)*lerr(1,1:3).*lerr(2,1:3), lerr(13,1:3), 'o-')
xlabel('1/h')
ylabel('||(B - B_0) \cdot n||_\infty/||B_0 \cdot n||_\infty')
hold on
for i = 1:2
    loglog(sqrt(3)*lerr(1,4*i+1:4*(i+1)).*lerr(2,4*i+1:4*(i+1)), ...
        lerr(10,4*i+1:4*(i+1)), 'o-')
end

h = 80:120;
plot(h, 80^5*8e-4*h.^(-5), 'b--')

h = 120:220;
plot(h, 120^7*1e-5*h.^(-7), 'r--')

h = 160:260;
plot(h, 160^9*8e-8*h.^(-9), 'k--')

legend('FMM+surfacefun, poly. order = n = 4', ...
    'FMM+surfacefun, n = 6', ...
    'FMM+surfacefun, n = 8', ...
    'O(h^{-3})', 'Location', 'southwest')
