n = 6;
nu = 3;
nv = 12;

figure(1)
[doms, qnodes, qweights] = prepare_square_torus(n, nu, nv, 1.0, 100);
plot(doms{1})
hold on
plot3(qnodes{1}(1,:), qnodes{1}(2,:), qnodes{1}(3,:), 'r.')
