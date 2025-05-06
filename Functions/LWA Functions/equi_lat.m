
function phi_e = equi_lat(gph,lat,lon,z_k)
    t = gph(1:end-1,1:end-1,:) <= z_k;
    A = abs(lon(1:end-1) - lon(2:end));
    B = abs(lat(1:end-1) - lat(2:end)).*(cos(lat(1:end-1))+cos(lat(2:end)));
    P = (A * B' / 2) .* t;
    phi_i = squeeze(sum(P,[1 2]));
    phi_e = asin(1-phi_i/(2*pi));
end

