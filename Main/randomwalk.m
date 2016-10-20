%% 3D random walk 

pos = zeros(10000,3);
r = zeros(9999,1);

for i = 2:10000;
    
    u1 = (rand()*2) -1;
    u2 = (rand()*2) -1;
    
    while (u1^2 + u2^2) >= 1;
        u1 = (rand()*2) -1;
        u2 = (rand()*2) -1;
    end
    
    pos(i,1) = pos(i-1,1) + (2*u1 * sqrt(1 - (u1^2 + u2^2)));
    pos(i,2) = pos(i-1,2) + (2*u2 * sqrt(1 - (u1^2 + u2^2)));
    pos(i,3) = pos(i-1,3) + (1 - 2*(u1^2 + u2^2));
    
    p = pos(i,1)^2 + pos(i,2)^2 + pos(i,3)^2;
    q = pos(i-1,1)^2 + pos(i-1,2)^2 + pos(i-1,3)^2;
    
    r(i-1) = p-q;
end

plot3(pos(:,1),pos(:,2),pos(:,3),'.')