function [ Flink ] = flink_func( link_dist )

T = 1;
k = 100*T;
Flink = (-1)*k*(abs(link_dist)-1);

end

