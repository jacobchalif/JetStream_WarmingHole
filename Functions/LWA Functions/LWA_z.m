function [phi_e,LWA_tot, LWA_a, LWA_c] = LWA_z(gph,lat,lon,z_k)
    % lat_deg/lon_deg are in radians
    % gph can be 2D with dims lat lon, or 3D with dims lat lon time

    phi_e = equi_lat(gph,lat,lon,z_k);
    r = 6.375e6; % radius of earth  (m)

    z_hat = gph - z_k;

    
    lt = permute(repmat(abs(lat) <= phi_e',[1 1 size(z_hat,1)]),[3 1 2]);
    cyc = z_hat <= 0 & lt;

    lt = permute(repmat(abs(lat) >= phi_e',[1 1 size(z_hat,1)]),[3 1 2]);
    anticyc = z_hat >= 0 & lt;
        
    latT = abs(lat(1:end-1) - lat(2:end)) .* (cos(lat(1:end-1))+cos(lat(2:end)))/2;
    LWA_ci = sum(cyc(:,1:end-1,:) .* latT',2);
    LWA_ai = sum(anticyc(:,1:end-1,:) .* latT',2);
    
    LWA_ai = repmat(LWA_ai,[1, size(gph,2), 1]);
    LWA_ci = repmat(LWA_ci,[1, size(gph,2), 1]);

    LWA_a = zeros(size(gph)); % anticyclonic wave activity
    LWA_c = zeros(size(gph)); % cyclonic wave activity
    
    C = permute(repmat(cos(phi_e),[1 size(gph,1) size(gph,2)]),[2 3 1]);
    LWA_a(anticyc) = r * z_hat(anticyc) .*  LWA_ai(anticyc);
    LWA_a = LWA_a ./ C;
    LWA_c(cyc) = r * -z_hat(cyc) .*  LWA_ci(cyc);
    LWA_c = LWA_c ./ C;

    LWA_tot = LWA_a + LWA_c;
end