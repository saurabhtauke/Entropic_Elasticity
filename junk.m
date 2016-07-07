clc


%% initialization Values

N1d = 10;                        %number of particles in one side of the square. This is initializing 100 particles in 2D square of side 10.

N = N1d*N1d;                     %number of particles
Density = 1;                     % density of system
L = power(N/Density,1/2);        % side L of square

position = zeros(N,3);           %array of position coordinates stored as X,Y,Z
velocity = zeros(N,3);           %array of velocities stored as VX,VY,VZ
force = zeros(N,3);              %force between any two non consecutive particles

pos2 = zeros(N,3);               %array of position coordinates stored as X,Y,

f_link = zeros(N,3);             %force between two consecutive particles

dt = 0.001;                                         % time step
iter = 1000;                                        % number of simulations


%% initialized positions

h = L/(N1d); % spacing between atoms

count = 1;

for i = 1 : N1d
    for j = 1 : N1d
            if (count <= N)
                position(count,1)= i*h;
                position(count,2)= j*h;
            end
             count = count + 1;
    end
end

 %plot(position(:,1), position(:,2),'*b') 
 %return
 
 count = 1;
 
 for r = 1:N1d          %to make the polymer snakelike
     
     if (mod(r,2)==1)
         
         for i = 1:N1d
     pos2((r*N1d)-(N1d-i),1) = position((r*N1d)-(N1d-i),1);
     pos2((r*N1d)-(N1d-i),2) = position((r*N1d)-(N1d-i),2);
        
         end
     end
     
     if (mod(r,2)==0)
         
        for i = 1:N1d 
     pos2((r*N1d)-(N1d-i),1) = position((r*N1d)-(i-1),1);
     pos2((r*N1d)-(N1d-i),2) = position((r*N1d)-(i-1),2);
     
        end
     end
 end
 
 %position
 %pos2
 %plot(pos2(:,1), pos2(:,2)) 
 %return
 
 %% Initialize Velocity

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

 % for force between two  non-linked molecules (i and i+2) 
dr = zeros(3,1);
potential = 0;

F0 = force_ij(1.22); % !!!!!!!!!!!!!!!check this
%.... have %commented out
%the F0 terms in the following part they represent the cutoff. since we
%have a E-12 decaying potential, th ecutoff is really not neccecary.

for i = 1 : N
    for j = i+2 : N
        
        dist = 0.0;
        
        for k = 1 : 3
            
            dr(k) = pos2(i,k) - pos2(j,k);
            
           if(dr(k) > L/2)               %PBC
               dr(k) = dr(k) - L;
           end
           
           if(dr(k) < -L/2)
               dr(k) = dr(k) + L;
           end
            
            dist = dist + dr(k)*dr(k);
        end
        
        dist = power(dist,0.5);
        
        if (dist <= 101)
            
            potential = potential + potential_ij(dist); % + F0*dist;
            F = force_ij(dist); % - F0;
           
            for k = 1 : 3
                force(i,k) = force(i,k) + F*dr(k)/dist;
                force(j,k) = force(j,k) - F*dr(k)/dist;   % because Fji = -Fij
            end
        end
        
    end
end

% force between two immediate linked molecules
dr = zeros(3,1);
potential = 0;
F2 = flink_func(1.22);

for i = 1 : N-1
    nxt = i+1;
    dist = 0.0;
        
        for k = 1 : 3
            
            dr(k) = pos2(i,k) - pos2(nxt,k);
            
           if(dr(k) > L/2)                   %PBC
               dr(k) = dr(k) - L;
           end
           
           if(dr(k) < -L/2)
               dr(k) = dr(k) + L;
           end
            
            dist = dist + dr(k)*dr(k);
        end
        
        dist = power(dist,0.5);
        
        if (dist <= 101)
            
            potential = potential + plink_func(dist); %+ F2*dist;
            F = flink_func(dist) ;% - F2;
           
            for k = 1 : 3
                f_link(i,k) = f_link(i,k) + F*dr(k)/dist;
                f_link(nxt,k) = f_link(nxt,k) - F*dr(k)/dist;   % because Fji = -Fij
            end
        end
        
end      

%force
%f_link

%% Force Iteration

total_E = zeros(iter,1);

for z = 1 : iter

    z;
old_force = force;

for i = 1 : N
    
    for j = 1 : 3
        
        pos2(i,j) = pos2(i,j) + velocity(i,j)*dt + 0.5*force(i,j)*dt*dt;       
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%force calculation for non-linked units
dr = zeros(3,1);
potential = 0;
force = zeros(N,3);

%.... have %commented out
%the F0 terms in the following part they represent the cutoff. since we
%have a E-12 decaying potential, th ecutoff is really not neccecary.

for i = 1 : N                          % CAREFUL ! we have not yet merged the two force arrays !!!!!
    for j = i+2 : N                  % have to be i+2 because we have i+1 linked with different restoring force
                                      % we also are not calculating the
                                      % force between two successive units
        dist = 0.0;
        
        for k = 1 : 3
            
            dr(k) = pos2(i,k) - pos2(j,k);
            
            if(dr(k) > L/2)                  %PBC
               dr(k) = dr(k) - L;
           end
           
           if(dr(k) < -L/2)
               dr(k) = dr(k) + L;
           end
           
            dist = dist + dr(k)*dr(k);
        end
        
        dist = power(dist,0.5);
        
        if (dist <= 101)     %101 was just a cutoff for earlier calculation. it can be set to any value now depending on your cutofff interaction distance
            
            potential = potential + potential_ij(dist); %+ F0*dist;
            F = force_ij(dist); %- F0;
            
            for k = 1 : 3
                force(i,k) = force(i,k) + F*dr(k)/dist;
                force(j,k) = force(j,k) - F*dr(k)/dist;
            end
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%force calculation for linked units
dr = zeros(3,1);
potential = 0;
F2 = flink_func(1.22);

for i = 1 : N-1
    nxt = i+1;
    dist = 0.0;
        
        for k = 1 : 3
            
            dr(k) = pos2(i,k) - pos2(nxt,k);
            
            if(dr(k) > L/2)                  %PBC
               dr(k) = dr(k) - L;
           end
           
           if(dr(k) < -L/2)
               dr(k) = dr(k) + L;
           end
            
            dist = dist + dr(k)*dr(k);
        end
        
        dist = power(dist,0.5);
        
        if (dist <= 101)
            
            potential = potential + plink_func(dist); %+ F2*dist;
            F = flink_func(dist) ;% - F2;
           
            for k = 1 : 3
                f_link(i,k) = f_link(i,k) + F*dr(k)/dist;
                f_link(nxt,k) = f_link(nxt,k) - F*dr(k)/dist;   % because Fji = -Fij
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
      
      plot3(pos2(:,1), pos2(:,2), pos2(:,3))
      pause(0.00005)
      
  end

end
figure
plot(total_E/N)
