function [GDOP,HDOP,PDOP,TDOP] = genmatrices(num_sat, sat_coords, receiver, UERE)

    all_coords = [];
    
    for i = 1:num_sat
        all_coords = [all_coords; sat_coords(i,1:3)];
    end
    
    Ri = all_coords - [0 0 0];
    
    S_R = all_coords - receiver;
    mag_ei = sqrt(sum(Ri.*Ri));
    ei = S_R./mag_ei;
    
    Gu = ones(num_sat,3);
    Gu = [(Gu + ei) ones(num_sat,1)];
    
    cov = inv(transpose(Gu)*Gu);
    
    GDOP = sqrt(cov(1,1) + cov(2,2) + cov(3,3) + cov(4,4))/UERE;
    PDOP = sqrt(cov(1,1) + cov(2,2) + cov(3,3))/UERE;
    HDOP = sqrt(cov(1,1) + cov(2,2))/UERE;
    TDOP = sqrt(cov(4,4))/UERE;
end