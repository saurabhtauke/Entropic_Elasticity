clc

% testing the potential and force between two objects

%% Initialize values

% x = 0;           %array of position coordinates stored as X,Y,Z
% v = 0;           %array of velocities stored as VX,VY,VZ
force = 0;              %force between any two non consecutive particles

%radius = 0.2;
%sigma = 1;


dt = 0.001;
iter = 10000;

position = zeros((dt*iter),2);
velocity = zeros((dt*iter),2);
kinetic = zeros((dt*iter),2);
potential = zeros((dt*iter),2);
tot_E = zeros((dt*iter),2);
%% initialize position

 x = 0;
 %plot3(position(:,1), position(:,2), position(:,3),'o-b')

%% provide a velocity 

vel_arb = 5;
v = vel_arb;

%% calculate the force
%velocity verlet

t=0;
f0 = test_pot_func(x);                  

x_nxt = x + (v*dt) + (0.5*f0*dt*dt);    % STEP 1

f_nxt = test_pot_func(x_nxt);           % STEP 2

v_nxt = v + (0.5*(f0 + f_nxt) *dt);     % STEP 3

x = x_nxt;
v = v_nxt;


%% Run the iteration

count = 0;

for i = 1:iter;
    
    count = count + 1;
    
    % verlet begins
    f0 = test_pot_func(x);                  

    x_nxt = x + (v*dt) + (0.5*f0*dt*dt);    % STEP 1

    f_nxt = test_pot_func(x_nxt);           % STEP 2

    v_nxt = v + (0.5*(f0 + f_nxt) *dt);     % STEP 3
    
    
    position(count,1) = abs(x-2);
    velocity(count,1) = v;
    
    kinetic(count,1) = 0.5* (velocity(count,1)^2);
    potential(count,1) = 0.5 *(position(count,1)^2);
    tot_E(count,1) = kinetic(count,1) + potential(count,1);
    
    % update position and velocity
    x = x_nxt;
    v = v_nxt;
    
    position(count,2) = t;
    velocity(count,2) = t;
    kinetic(count,2) = t;
    potential(count,2) = t;
    tot_E(count,2) = t;
    
    
    t = t + dt;
end


hold on
%plot (position(:,2),position(:,1),'b+')
%plot (velocity(:,2),velocity(:,1),'rO')
plot (kinetic(:,2),kinetic(:,1),'g*')
plot (potential(:,2),potential(:,1),'r-')
plot (tot_E(:,2),tot_E(:,1),'b-')
hold off
