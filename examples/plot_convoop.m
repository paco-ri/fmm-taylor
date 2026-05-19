% Plot convergence results from convoop.m
% This script loads the lerr array and creates convergence plots

% Load the lerr array (assumes it was saved as 'lerr' in a .mat file)
% You may need to adjust the filename or location depending on where lerr is stored
% load('lerr.mat', 'lerr', 'ns', 'nvs');
load('jupiter_data/save_convoop_lerr.mat', 'lerr')

% If ns and nvs are not available, specify them here based on your run:
ns = [5 7 9 11];
nvs = [6 8 10 12];

% Create figure
figure('Position', [100 100 900 600]);

% Plot convergence curves
loglog(sqrt(lerr(1,1:size(nvs,2)).^2.*lerr(2,1:size(nvs,2)).* ...
    lerr(3,1:size(nvs,2))), lerr(4,1:size(nvs,2)), 'o-')
hold on
for i = 1:size(ns,2)-1
    ind = i*size(nvs,2)+1:(i+1)*size(nvs,2);
    loglog(sqrt(lerr(1,ind).^2.*lerr(2,ind).* ...
        lerr(3,ind)), lerr(4,ind), 'o-')
end

% Add reference convergence lines
h1 = 50:90;
h2 = 80:120;
h3 = 100:150;
loglog(h1,1e1.*h1.^(-1),'--')
loglog(h2,1e4.*h2.^(-3),'--')
loglog(h3,1e7.*h3.^(-5),'--')

% Labels and legend
xlabel('sqrt(n^2 * nv * nu)', 'Interpreter', 'latex')
ylabel('Relative Error', 'Interpreter', 'latex')
title('Convergence Study', 'Interpreter', 'latex')
grid on
legend('p=4','p=6','p=8','O(h^{-1})','O(h^{-3})','O(h^{-5})', 'Location', 'best')
