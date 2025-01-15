% [domo, nodes, wts] = alt_stellarator(5,15,5,1);
% domi = alt_stellarator(5,15,5,.6);
[doms, nodes, wts] = prepare_stellarator(5,15,5,5,15,5,1,.3,4,6,6);
domo = doms{1};
domi = doms{2};
plot(domo)
alpha .5
hold on
plot(domi)
plot3(nodes{1}(1,:), nodes{1}(2,:), nodes{1}(3,:), '.')
plot3(nodes{2}(1,:), nodes{2}(2,:), nodes{2}(3,:), '.')