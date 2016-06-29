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
 
 for r = 1:N1d
     
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
 plot(pos2(:,1), pos2(:,2)) 
 %return
 
 %% Initialize Velocity
 %% evaluation of forces

 % for force between two  non-linked molecules (i and i+2) 
dr = zeros(3,1);
potential = 0;

F0 = force_ij(1.22); % !!!!!!!!!!!!!!!check this

%.... have %commented out
%the F0 terms in the following part they represent the cutoff. since we
%have a E-12 decaying potential, the cutoff is really not neccecary.

for i = 2 : N-1       %this is becasue 
    for j = i+2 : N
        
        dist = 0.0;
        
        for k = 1 : 3
            
            dr(k) = pos2(i,k) - pos2(j,k);
            
           if(dr(k) > L/2)                 %PBC
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
            
           if(dr(k) > L/2)         %PBC
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
%% Plot
