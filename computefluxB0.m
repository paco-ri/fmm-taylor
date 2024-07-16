function flux = computefluxB0(domrmaj,domrmin,nr,nt,ntheta,rmin,rmaj,jmag,lambda)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

flux = 0;
shift = [domrmaj; 0; 0];
[rr,wr] = chebpts(nr,[0 domrmin]);
tt = 2*pi.*(0:nt-1)./nt;
wt = 2*pi/nt;
for i = 1:nr
    for j = 1:nt
        targpt = [rr(i)*cos(tt(j)); 0; rr(i)*sin(tt(j))];
        targpt = targpt + shift;
        if true
            plot3(targpt(1),targpt(2),targpt(3),'ob')
        end
        B0eval = reftaylor(ntheta,rmin,rmaj,jmag,lambda,targpt);
        flux = flux - B0eval(2)*rr(i)*wt*wr(i);
    end
end

end