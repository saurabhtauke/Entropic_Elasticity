clc

% SaurabhTauke

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

dt = 0.0001;                                         % time step
iter = 100000;                                        % number of simulations

total_energy = zeros(iter,1);    % had to  pre-allocate ... can be allocated in the force iteration loop tho.
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
 %plot(pos2(:,1), pos2(:,2)) 
 %return
 
 %% Initialize Velocity
 
for i = 1 : 3*N
    velocity(i) = (5)*(rand-0.5);
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

F0 = force_ij(1.33); % !!!!!!!!!!!!!!!check this
%F2 = flink_func(1.22);

%.... have %commented out
%the F0 terms in the following part they represent the cutoff. since we
%have a E-12 decaying potential, th ecutoff is really not neccecary.




for i = 1 : N
    for j = 1 : N
        
        dist = 0.0;  %  calculating the distance
        
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
        
        dindex= abs(i-j);
        
        if (dindex==0)
            force(i,k) =  0;
            force(j,k) =  0;
                
        elseif (dindex==1)            % force between two immediate linked molecules
            
            potential = potential + plink_func(dist); 
            F = flink_func(dist) ;
           
            for k = 1 : 3
                f_link(i,k) = f_link(i,k) + F*dr(k)/dist;
                f_link(j,k) = f_link(j,k) - F*dr(k)/dist;   % because Fji = -Fij
            end
        elseif (dindex > 1 && dist <= 100)       % for force between two  non-linked molecules (i and i+2)  and   %The cutoff distance  
            
            potential = potential + potential_ij(dist) + F0*dist;
            F = force_ij(dist) - F0;
           
            for k = 1 : 3
                force(i,k) = force(i,k) + F*dr(k)/dist;
                force(j,k) = force(j,k) - F*dr(k)/dist;   % because Fji = -Fij
            end
        end
        
    end
end

%force
%f_link

%% Force Iteration

for z = 1 : iter
    
        z;
    old_force = force;

    for i = 1 : N

        for j = 1 : 3

            pos2(i,j) = pos2(i,j) + velocity(i,j)*dt + 0.5*force(i,j)*dt*dt;       

        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    dr = zeros(3,1);
    potential = 0;

    force = zeros(N,3);

    for i = 1 : N
        for j = 1 : N


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

            dist = power(dist,0.5);        %this is the distance between ith and jth molecules

            dindex= abs(i-j);

            % same force calculation code as above but THE FORCE ARRAY AND
            % FLINK ARRAY HAVE BEEN MERGED INTO FORCE ARRAY

            if (dindex==0)
                force(i,k) =  0;
                force(j,k) =  0;

            elseif (dindex==1)            % force between two immediate linked molecules
                
                F2 = flink_func(1.22);
                potential = potential + plink_func(dist); 
                F = flink_func(dist) ;% - F2;

                for k = 1 : 3
                    force(i,k) = force(i,k) + F*dr(k)/dist;
                    force(j,k) = force(j,k) - F*dr(k)/dist;   % because Fji = -Fij
                end
                

%                 for k = 1 : 3
%                     force(i,k) = 0;
%                     force(j,k) = 0;   % because Fji = -Fij
%                 end

            elseif (dindex > 1 && dist <= 1.33)       % for force between two  non-linked molecules (i and i+2)  and   %The cutoff distance  
                
                F0 = force_ij(1.33);
                potential = potential + potential_ij(dist) + F0*dist;
                F = force_ij(dist) - F0;

                for k = 1 : 3
                    force(i,k) = force(i,k) + F*dr(k)/dist;
                    force(j,k) = force(j,k) - F*dr(k)/dist;   % because Fji = -Fij
                end

%                   for k = 1 : 3
%                     force(i,k) = 0;
%                     force(j,k) = 0;
%                   end

            end
        end
    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %  kinetic = 0;
    %
     
    
v_sqr = 0.0;

for i = 1 : N
    vtemp = 0.0;
            for j = 1 : 3

                velocity(i,j) = velocity(i,j) + 0.5*(force(i,j) + old_force(i,j))*dt;
	    	    vtemp = vtemp + (velocity(i,j)*velocity(i,j));
            end  
    	v_sqr = v_sqr + vtemp;   
end

temperature = v_sqr/(N*3.0);
	

if (mod(z,500) == 0)
		
                velocity = velocity .* (sqrt(1/temperature));

                v_sqr = 0.0;
            for i = 1 : N
				vtemp = 0.0;

            		for j = 1 : 3
						vtemp = vtemp + (velocity(i,j)*velocity(i,j));
            		end 
 
				v_sqr = v_sqr + vtemp;   
            end
                z
                temperature = v_sqr/(N*3.0)
end


total_energy(z) = (potential + (0.5*v_sqr))/(N);

%      % for making the movie
%      if(rem(i,5) == 0)
% 
%           plot3(pos2(:,1), pos2(:,2), pos2(:,3))
%           pause(0.000005)
% 
%      end

end

figure
plot(total_energy/N)
