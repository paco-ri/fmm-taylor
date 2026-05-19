% load B and jupiter B

load('jupiter_data/B.mat', 'B');
B_jupiter = B{1};
load('B.mat', 'B');
B = B{1};

% Create 3 subplots
figure;

%{
subplot(1, 3, 1);
plot(norm(B));
title('norm(B)');
ylabel('Norm value');

subplot(1, 3, 2);
plot(norm(B_jupiter));
title('norm(B_jupiter)');
ylabel('Norm value');

subplot(1, 3, 3);
plot(norm(B - B_jupiter));
title('norm(B - B_jupiter)');
ylabel('Norm value');
%}

plot(norm(B - B_jupiter));
colorbar;

