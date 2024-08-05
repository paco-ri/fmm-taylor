CIMS = false;

if CIMS
    low = load('lerr_saved/circlelow260724.mat').lerr;
    high = load('lerr_saved/circlehigh260724.mat').lerr;
else
    low = load('circlelow310724.mat').lerr;
    high = load('circlehigh310724.mat').lerr;
    lowstat = load('circlelowstatic310724.mat').lerr;
    highstat = load('circlehighstatic310724.mat').lerr;
end

%  7 sigma
%  8 alpha
%  9 m0
% 10 m
% 11 gradSsigma
% 12 curlSm
% 13 gradSsigmaint
% 14 curlSmint 
% 15 Smint

h = sqrt(low(1,:).^2.*low(2,:).*low(3,:));
err = low(7:end,:) - high(7:end);

hstat = sqrt(lowstat(1,:).^2.*lowstat(2,:).*lowstat(3,:));
errstat = lowstat(7:end,:) - highstat(7:end);

figure(1)
loglog(h,abs(err(1,:)),'o-')
hold on
% semilogy(h,abs(low(5,:)),'o--')
for i = 2:8
    loglog(h,abs(err(i,:)),'o-')
end
hh = 40:80;
loglog(hh,1e4.*hh.^(-4))
legend('\sigma','\alpha','m_0','m','\nabla S[\sigma]', ...
    '\nabla \times S[m]','\nabla S[\sigma] int', ...
    '\nabla \times S[m] int')%,'S[m] int')

figure(2)
ind = 17-6;
loglog(hstat,abs(errstat(ind,:)),'o-')
hold on
loglog(h,abs(err(ind,:)),'o-')
% loglog(hstat,abs(lowstat(4,:)),'--')
% loglog(h,abs(low(4,:)),':')
