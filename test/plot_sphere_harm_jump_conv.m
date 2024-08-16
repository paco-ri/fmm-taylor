% load('lerr_sph.mat')

for nind = 1:size(lerr_sph,1)
    for nrefind = 1:size(lerr_sph,2)
        loglog(lerr_sph(nind,nrefind,1)^2*4^lerr_sph(nind,nrefind,2), lerr_sph(nind,nrefind,8), 'ok')
        hold on
    end
end
