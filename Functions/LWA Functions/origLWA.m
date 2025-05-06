function LWA = origLWA(gph,lat,lon,dayi)
    LWA = zeros(size(gph,1),size(gph,2),3);
    phi_i = pi/2;

    latr = radians(lat); 
    lonr = radians(lon);

    for z = 4800:10:6000
        [phi_e,LWA_t,LWA_a,LWA_k] = calcLWA(gph(:,:,dayi),latr,lonr,z);
        diff = (phi_i-phi_e);
        LWA(:,:,1) = LWA(:,:,1) + LWA_t*diff ; 
        LWA(:,:,2) = LWA(:,:,2) + LWA_a*diff ; 
        LWA(:,:,3) = LWA(:,:,3) + LWA_k*diff ;
        phi_i = phi_e;
    end
end

function [phi_e, LWA_tot, LWA_a, LWA_k] = calcLWA(gph,lat,lon,z_k)
    phi_e = equi_lat(gph,lat,lon,z_k);
    r = 6.375e6; % radius of earth  (m)
    
    LWA_ai = zeros(size(gph,1),1); % anticyclonic wave activity
    LWA_ki = zeros(size(gph,1),1); % cyclonic wave activity
    z_hat = zeros(size(gph,1),size(gph,2));

    for loni = 1:size(gph,1)
        for lati = 1:size(gph,2)-1
            z_hat(loni,lati) = gph(loni,lati) - z_k;
            if (z_hat(loni,lati) >= 0) && (lat(lati) >= phi_e)
                LWA_ai(loni) = LWA_ai(loni) + abs(lat(lati) - lat(lati + 1)) * (cos(lat(lati))+cos(lat(lati+1)))/2;
            elseif (z_hat(loni,lati) <= 0) && (lat(lati) <= phi_e)
                LWA_ki(loni) = LWA_ki(loni) + abs(lat(lati) - lat(lati + 1)) * (cos(lat(lati))+cos(lat(lati+1)))/2;
            end
        end
    end

    LWA_a = zeros(size(gph,1),size(gph,2)); % anticyclonic wave activity
    LWA_k = zeros(size(gph,1),size(gph,2)); %cyclonic wave activity
    for loni = 1:size(gph,1)
        for lati = 1:size(gph,2)-1
            if (z_hat(loni,lati) >= 0) && (lat(lati) >= phi_e)
                LWA_a(loni,lati) = z_hat(loni,lati) * r * LWA_ai(loni) / cos(phi_e);
            elseif (z_hat(loni,lati) <= 0) && (lat(lati) <= phi_e)
                LWA_k(loni,lati) = -z_hat(loni,lati) * r * LWA_ki(loni) / cos(phi_e);
            end
        end
    end
    LWA_tot = LWA_a + LWA_k;
end


function phi_e = equi_lat(gph,lat,lon,z_k)
    phi_i = 0;

    for loni = 1:size(gph,1)-1
        for lati = 1:size(gph,2)-1
            if gph(loni,lati) <= z_k
                phi_i = phi_i + abs(lat(lati) - lat(lati + 1))*abs(lon(loni) - lon(loni + 1))*(cos(lat(lati))+cos(lat(lati+1)))/2;
            end
        end
    end
    phi_e = asin(1-phi_i/(2*pi));
end

function r = radians(degrees)
    r= degrees*pi/180;
end