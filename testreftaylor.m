% check that reftaylor actualy produces a Taylor state. it is. 

ntheta = 1e3;
rmin = 2.0;
rmaj = 2.0;
jmag = 1.0;
lambda = 1.0;

phi = 4*pi/7;
h = 1e-3;
center = [rmaj*cos(phi)+.2 rmaj*sin(phi)-.1 .5];
targs = [center;
    center + [h 0 0];
    center - [h 0 0];
    center + [0 h 0];
    center - [0 h 0];
    center + [0 0 h];
    center - [0 0 h]];

B0int = zeros(size(targs));
for i = 1:7
    B0int(i,:) = reftaylor(ntheta,rmin,rmaj,jmag,lambda,targs(i,:));
end

curlB0 = [B0int(6,2)-B0int(7,2)-B0int(4,3)+B0int(5,3); % DyBz-DzBy
    B0int(2,3)-B0int(3,3)-B0int(6,1)+B0int(7,1); % DzBx-DxBz
    B0int(4,1)-B0int(5,1)-B0int(2,2)+B0int(3,2)]; % DxBy-DyBx
curlB0 = curlB0./(2*h);

lambdaB0 = lambda.*B0int(1,:).';
disp(norm(curlB0-lambdaB0))
