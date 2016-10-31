clc

% SaurabhTauke

%% initialization Values

N_monomer = 10000;               % Number of monomer units

pos = zeros(N_monomer,3);        %array of position coordinates stored as X,Y,Z
velocity = zeros(N_monomer,3);           %array of velocities stored as VX,VY,VZ
force = zeros(N_monomer-1,3);              %force between any two non consecutive particles


f_link = zeros(N_monomer-1,3);             %force between two consecutive particles

dt = 0.0001;                                         % time step
iter = 10000;                                        % number of iterations in simulation

total_energy = zeros(iter,1);    % pre-allocated ... can be allocated in the force iteration loop tho.

%% 3D random walk 

r = zeros(N_monomer-1,1);

for i = 2 : N_monomer;
    
    u1 = (rand()*2) -1;
    u2 = (rand()*2) -1;
    
    while (u1^2 + u2^2) >= 1;
        u1 = (rand()*2) -1;
        u2 = (rand()*2) -1;
    end
    
    pos(i,1) = pos(i-1,1) + (2*u1 * sqrt(1 - (u1^2 + u2^2)));
    pos(i,2) = pos(i-1,2) + (2*u2 * sqrt(1 - (u1^2 + u2^2)));
    pos(i,3) = pos(i-1,3) + (1 - 2*(u1^2 + u2^2));
    
    r(i-1) = (pos(i,1)-pos(i-1,1))^2 + (pos(i,2)-pos(i-1,2))^2 + (pos(i,3)-pos(i-1,3))^2 ;
    
    
end

plot3(pos(:,1),pos(:,2),pos(:,3),'.')
