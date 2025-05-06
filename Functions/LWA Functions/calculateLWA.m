
function LWA = calculateLWA(gph,lat_deg,lon_deg)
    % gph can be 2D with dims lat lon, or 3D with dims lat lon time
    % lat_deg/lon_deg are in degrees

    LWA = zeros(size(gph,1),size(gph,2),size(gph,3),3);
    phi_i = pi/2;

    lat = radians(lat_deg); 
    lon = radians(lon_deg);
    
    for z = 4800:10:6000
        [phi_e,LWA_tot, LWA_a, LWA_c] = LWA_z(gph,lat,lon,z);
        diff = permute(repmat(phi_i-phi_e,[1 size(LWA,1) size(LWA,2)]),[2 3 1]);
        LWA(:,:,:,1) = LWA(:,:,:,1) + LWA_tot.*diff; 
        LWA(:,:,:,2) = LWA(:,:,:,2) + LWA_a.*diff;
        LWA(:,:,:,3) = LWA(:,:,:,3) + LWA_c.*diff ;
        phi_i = phi_e;
    end
end
