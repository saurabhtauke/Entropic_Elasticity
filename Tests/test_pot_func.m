function [ spring_force ] = test_pot_func( dist )

    eq_length = 2;
    k = 1;
    spring_force = -1*k*(dist - eq_length);
end

