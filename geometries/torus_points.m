function X = torus_points(u, v, s)
%TORUS_POINTS Points on the scale-s member of the prepare_torus family.
%   X = TORUS_POINTS(U, V, S) returns a 3 x numel(U) array of the points with
%   toroidal angle U, poloidal angle V on the torus prepare_torus builds at
%   scale S. U, V, S are arrays of the same size, or S is a scalar. Every
%   term except the fixed major radius scales with S, so S < 1 lies inside
%   the scale-1 torus and ai < S < ao lies in the shell between two of them.
%
%   The coefficients must match prepare_torus.m.

d = [0.17 0.11  0    0;
     0    1     0.01 0;
     0    2     0    0;
     0   -0.25  0    0];

u = u(:).';
v = v(:).';
s = s(:).' + 0*u;
R = 2 + 0*u;
z = 0*u;
for i = -1:2
    for j = -1:2
        if (i == 1 && j == 0) || d(i+2,j+2) == 0, continue, end
        a = (1-i)*v + j*u;
        R = R + s*d(i+2,j+2).*cos(a);
        z = z + s*d(i+2,j+2).*sin(a);
    end
end
X = [R.*cos(u); R.*sin(u); z];
end
