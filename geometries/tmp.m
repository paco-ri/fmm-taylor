fid = fopen('qas.txt');
qasmodes = textscan(fid, '%f %f %f');
fclose(fid);

fid = fopen('qascoils.txt');
coilmodes = textscan(fid, '%f %f %f');
fclose(fid);

ns = [5 5];
nus = [15 15];
nvs = [5 5];
nr = 16;
nt = 40;
np = nt;
modes = {coilmodes, qasmodes};
[doms, nodes, weights] = prepare2(ns, nus, nvs, nr, nt, np, modes);

plot(doms{1})
alpha 0.5
hold on
plot(doms{2})
plot3(nodes{1}(1,:),nodes{1}(2,:),nodes{1}(3,:),'.')
plot3(nodes{2}(1,:),nodes{2}(2,:),nodes{2}(3,:),'.')