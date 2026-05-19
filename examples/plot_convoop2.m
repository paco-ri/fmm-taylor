% Plot convergence results from convoop2.m
% This script loads the lerr array (and ns,nvs) and creates the convergence plot
% Usage: run convoop2.m to produce lerr, then save('lerr.mat','lerr','ns','nvs')
% and run this script.

% Try to load variables; prefer jupiter_data/save_convoop2_lerr.mat
% but fall back to lerr.mat if not present
if exist('lerr','var') ~= 1
    if exist('jupiter_data/save_convoop2_lerr.mat','file')
        S = load('jupiter_data/save_convoop2_lerr.mat');
    elseif exist('lerr.mat','file')
        S = load('lerr.mat');
    else
        S = struct();
    end
    if isfield(S,'lerr'), lerr = S.lerr; end
    if isfield(S,'ns'), ns = S.ns; end
    if isfield(S,'nvs'), nvs = S.nvs; end
end

if ~exist('lerr','var')
    error('lerr not found. Run convoop2.m and save lerr.mat, or put lerr in workspace.');
end

if ~exist('ns','var')
    ns = unique(lerr(1,:));
end
if ~exist('nvs','var')
    nvs = unique(lerr(2,:));
end

figure('Position',[100 100 900 600]);

% Primary plot: one curve per ns (polynomial order)
loglog(sqrt(lerr(1,1:size(nvs,2)).^2 .* lerr(2,1:size(nvs,2)) .* lerr(3,1:size(nvs,2))), ...
    lerr(4,1:size(nvs,2)), 'o-')
hold on
for i = 1:size(ns,2)-1
    ind = i*size(nvs,2)+1:(i+1)*size(nvs,2);
    loglog(sqrt(lerr(1,ind).^2 .* lerr(2,ind) .* lerr(3,ind)), lerr(4,ind), 'o-')
end

% Reference lines tuned to convoop2
h1 = 34:68;
h2 = 50:90;
h3 = 62:124;
loglog(h1,5e1 .* h1.^(-2), '--')
loglog(h2,9e4 .* h2.^(-4), '--')
loglog(h3,1e8 .* h3.^(-6), '--')

legend('p=poly. order=4','p=6','p=8','O(h^{-2})','O(h^{-4})','O(h^{-6})','location','southwest')
xlabel('h, panel size')
ylabel('||B-B_0||_\infty / ||B_0||_\infty')
title('Convergence study (convoop2)')
grid on

% Optional: save figure
% saveas(gcf,'convoop2_convergence.png')

fprintf('Plotted %d data points from lerr (size %dx%d)\n', size(lerr,2), size(lerr,1), size(lerr,2));
