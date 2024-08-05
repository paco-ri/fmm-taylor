% norder ipars(1) err

% eps 1e-8 -- probably should set lower
errs = [
    6 1 0.81589e-4;
    6 2 0.20495e-5;
    6 3 0.20498e-7;
    7 1 0.25407e-4;
    7 2 0.11157e-5;
    7 3 0.12127e-8];

hs = 4.^errs(1:3,2);
loglog(hs,errs(1:3,3),'o-')
hold on
loglog(hs(1):hs(end),(hs(1):hs(end)).^(-5),'--')
hs = 4.^errs(4:6,2);
loglog(hs,errs(4:6,3),'o-')
loglog(hs(1):hs(end),(hs(1):hs(end)).^(-6),'--')