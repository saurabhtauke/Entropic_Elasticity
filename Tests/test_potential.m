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

position = zeros((dt*iter),1);
velocity = zeros((dt*iter),1);
kinetic = zeros((dt*iter),1);
potential = zeros((dt*iter),1);
tot_E = zeros((dt*iter),1);
time = zeros(dt*iter,1);

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
    
    
    position(count) = abs(x-2);
    velocity(count) = v;
    
    kinetic(count) = 0.5* (velocity(count)^2);
    potential(count) = 0.5 *(position(count)^2);
    tot_E(count) = kinetic(count) + potential(count);
    
    % update position and velocity
    x = x_nxt;
    v = v_nxt;
   
    time(count) = t;
    t = t + dt;
end


hold on
%plot (position, time,'b+')
%plot (velocity,time,'rO')
plot (kinetic,time,'g*')
plot (potential,time,'r-')
plot (tot_E,time,'b-')
hold off
