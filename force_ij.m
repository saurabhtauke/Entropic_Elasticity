function [ Fij ] = force_ij( distance )

Fij = 4*(12*(power(0.015,13))/(power(distance,13)));

end

