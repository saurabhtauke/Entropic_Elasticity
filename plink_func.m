function [ Plink ] = plink_func( link_dist)

T = 1;
k =  1000*T;
Plink = 0.5*k*(((link_dist)-1)^2);
end

