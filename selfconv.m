low = load('lerr_saved/circlelow260724.mat').lerr;
high = load('lerr_saved/circlehigh260724.mat').lerr;

%  7 sigma
%  8 alpha
%  9 m0
% 10 m
% 11 gradSsigma
% 12 curlSm
% 13 gradSsigmaint
% 14 curlSmint 

h = sqrt(low(1,:).^2.*low(2,:).*low(3,:));
err = low(7:14,:) - high(7:14);

ind = 8-6;

semilogy(h,abs(err(ind,:)),'o-')
hold on
semilogy(h,abs(low(6,:)),'o-')