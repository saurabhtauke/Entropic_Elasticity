clc


%% MD parameters

N_1D = 5;                                          % number of particles in 1-D
N = power(N_1D,3);                                  % number of particles
Density = 0.5;                                      % density of system
L = power(N/Density,1/3);                           % edge L of cube
T = 2;                                              % temperature of the system
dt = 0.001;                                         % time step
iter = 1000;                                        % number of simulations


position = zeros(N,3);
force = zeros(N,3);
velocity = zeros(N,3);

%potential = 0;
kinetic = 0;
%% initialisation of particle positions

h = L/(N_1D + 2); % spacing between atoms

count = 1;

for i = 1 : N_1D
    for j = 1 : N_1D
        for k = 1 : N_1D
            if (count <= N)
                position(count,1)= i*h;
                position(count,2)= j*h;
                position(count,3)= k*h;
            end
            count = count + 1;
        end
    end
end

 %plot3(position(:,1), position(:,2), position(:,3),'b*')
 %return
%% initialisation of particle velocities


for i = 1 : 3*N
    velocity(i) = (12^0.5)*(rand-0.5);
end

vx = sum(velocity(:,1));
vy = sum(velocity(:,2));
vz = sum(velocity(:,3));

for i = 1 : N
    velocity(i,1) =  velocity(i,1) - (vx/N);
    velocity(i,2) =  velocity(i,2) - (vy/N);
    velocity(i,3) =  velocity(i,3) - (vz/N);
end

%% evaluation of forces

dr = zeros(3,1);
potential = 0;
F0 = force_ij(2.5);

for i = 1 : N
    for j = i+1 : N
        
        dist = 0.0;
        
        for k = 1 : 3
            
            dr(k) = position(i,k) - position(j,k);
            
            if(dr(k) > L/2)
                dr(k) = dr(k) - L;
            end
            
            if(dr(k) < -L/2)
                dr(k) = dr(k) + L;
            end
            
            dist = dist + dr(k)*dr(k);
        end
        
        dist = power(dist,0.5);
        
        if (dist <= 2.5)
            
            potential = potential + potential_ij(dist) + F0*dist;
            F = force_ij(dist) - F0;
            
            for k = 1 : 3
                force(i,k) = force(i,k) + F*dr(k)/dist;
                force(j,k) = force(j,k) - F*dr(k)/dist;
            end
        end
        
    end
end

%% iterations

total_E = zeros(iter,1);

for z = 1 : iter

    z;
old_force = force;

for i = 1 : N
    
    for j = 1 : 3
        
        position(i,j) = position(i,j) + velocity(i,j)*dt + 0.5*force(i,j)*dt*dt;     
        
  
        if (position(i,j) > L)
            position(i,j) = position(i,j) - L;
        end
        
        if (position(i,j) < 0)
            position(i,j) = position(i,j) + L;
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr = zeros(3,1);
potential = 0;
force = zeros(N,3);

for i = 1 : N
    for j = i+1 : N
        
        dist = 0.0;
        
        for k = 1 : 3
            
            dr(k) = position(i,k) - position(j,k);
            
            if(dr(k) > L/2)
                dr(k) = dr(k) - L;
            end
            
            if(dr(k) < -L/2)
                dr(k) = dr(k) + L;
            end
            
            dist = dist + dr(k)*dr(k);
        end
        
        dist = power(dist,0.5);
        
        if (dist <= 2.5)
            
            potential = potential + potential_ij(dist) + F0*dist;
            F = force_ij(dist) - F0;
            
            for k = 1 : 3
                force(i,k) = force(i,k) + F*dr(k)/dist;
                force(j,k) = force(j,k) - F*dr(k)/dist;
            end
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kinetic = 0;

for i = 1 : N
        for j = 1 : 3
            
            velocity(i,j) = velocity(i,j) + 0.5*(force(i,j) + old_force(i,j))*dt;
            kinetic = kinetic + 0.5*power(velocity(i,j),2);
        end     
end

total_E(z) = potential + kinetic;

  if(rem(i,5) == 0)
      
      plot3(position(:,1), position(:,2), position(:,3),'b*')
      pause(0.00005)
      
  end

end
figure
plot(total_E/N)

