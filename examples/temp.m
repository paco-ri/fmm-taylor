% [domo, nodes, wts] = alt_stellarator(5,15,5,1);
% domi = alt_stellarator(5,15,5,.6);

n = 5;
nv = 5;
nu = 3*nv;
[doms, nodes, wts] = prepare_stellarator(n,nu,nv,n,nu,nv,1,.6,16,40,40);
domo = doms{1};
domi = doms{2};
disp(sum(wts{1}))
disp(sum(wts{2}))
plot(domo)
alpha .5
hold on
plot(domi)

% plot3(nodes{1}(1,:), nodes{1}(2,:), nodes{1}(3,:), '.')
% plot3(nodes{2}(1,:), nodes{2}(2,:), nodes{2}(3,:), '.')

d = [0.15  0.09  0.00  0.00  0.00;
     0.00  1.00  0.03 -0.01  0.00;
     0.08  4.00 -0.01 -0.02  0.00;
     0.01 -0.28 -0.28  0.03  0.02;
     0.00  0.09 -0.03  0.06  0.00;
     0.01 -0.02  0.02  0.00 -0.02];
scale = 0.6;%1.0;
d = scale.*d;
d(3,2) = 2.00;
dsum = sum(d,2);

dfun1 = @(v) 0*v;
dfun2 = @(v) 0*v;
for j=-1:4
    jj=j+2;
    dfun1 = @(v) dfun1(v) + dsum(jj)*cos(v*(1-j));
    dfun2 = @(v) dfun2(v) + dsum(jj)*cos(v*(1-j))*(1-j);
end
dfun = @(v) dfun1(v)*dfun2(v);
dchebfun = chebfun(@(v) dfun(v),[0,2*pi]);
I = sum(dchebfun);
disp(I)

dsum = sum(d);
dfun1 = @(u) 0*u;
dfun2 = @(u) 0*u;
Q = 3;
for k=-1:3
    kk=k+2;
    dfun1 = @(u) dfun1(u) + dsum(kk)*cos(u*k*Q)*cos(u);
    dfun2 = @(u) dfun2(u) + dsum(kk)*(cos(u*k*Q)*cos(u)-k*Q*sin(u*k*Q)*sin(u));
end
dfun = @(u) dfun1(u)*dfun2(u);
dchebfun = chebfun(dfun,[0,2*pi]);
I = sum(dchebfun);
disp(I)