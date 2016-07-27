program md
	Implicit none
	Integer:: n_part, i,j, iter, k, counter, jj
	real*8:: t ,dt, A,avg_v, lx,ly,lz, dx, dt_sqby2, rc, dist_x, dist_y, dist_z, f0, v0, poten, vprime0, total_energy
	real*8:: potenprime, poten_dp, sumvx, sumvy, sumvz, avgv2, temp, dist, sig, v_sqr, total_poten, dr, r
	real*8, dimension(:), allocatable:: x,y,z, vx,vy,vz, fx, fy, fz, poten1, x2,y2,z2, vx2,vy2,vz2, fx2, fy2, fz2, avg_v_max_boltz
	real*8, dimension(:,:), allocatable:: distance, force1
	real*8, dimension(:,:), allocatable:: max_boltz, avg_n_r
	!real*8, dimension(:,:,:), allocatable :: n_r
	Character(len=15) :: name1, name2
	!*******************Initialise********************

	print*, 'Enter no of particles < 1000'
	read*, n_part
	Print*, 'enter dt'
	read*, dt
	Print*, 'Enter no. of iterations > 15000'
	read*, iter
	dt_sqby2 = (dt**2)/2
	sig = 1.0d0
	lx = 15.0d0
	ly = 15.0d0
	lz = 15.0d0
	rc = 2.5d0
	dr = 0.2
	
	allocate(x(n_part),y(n_part),z(n_part), vx(n_part),vy(n_part),vz(n_part), x2(n_part), vx2(n_part))
	allocate(distance(n_part,n_part), force1(n_part, n_part), y2(n_part), z2(n_part),vy2(n_part), vz2(n_part))
	allocate(fx(n_part), fy(n_part), fz(n_part), fx2(n_part), fy2(n_part), fz2(n_part), poten1(n_part))
	allocate(max_boltz(iter,n_part), avg_v_max_boltz(n_part))
	
	open(90, file = 'position.dat')
	open(91, file = 'temp.dat')
	open(92, file = 'force.dat')
	open(93, file = 'energy.dat')
	open(94, file = 'momentum.dat')
	open(95,file = 'mbstat.dat')
	!name1 = 'struc.dat'
	!call addnumtostring(name1,n_part)
	!open(96, file = name1)
	
	max_boltz(:,:) = 0.0d0
	!n_r(:,:,:) = 0
	
	f0 = 4*(12*(sig**12)*1.0d0/(rc**13) - 6*(sig**6)*1.0d0/(rc**7))
	v0 = 4*((sig**12)*1.0d0/(rc**12) - (sig**6)*1.0d0/(rc**6))
	vprime0 = v0 + f0*rc

	do i= 1,n_part
		x(i) = mod(int((i)),int(lx))
		y(i) = mod(int((i)/lx),int(ly))
		z(i) = mod(int((i)/(ly*lx)),int(lz))
		write(90,*) x(i), y(i), z(i)
	End do


	A = dsqrt(12.0d0)
	
	sumvx = 0.0d0; sumvy = 0.0d0; sumvz = 0.0d0

	Do i = 1, n_part
		vx(i) = A*(rand() - 0.5d0) ; sumvx = sumvx + vx(i)
		vy(i)= A*(rand() - 0.5d0) ; sumvy = sumvy + vy(i)
		vz(i)= A*(rand() - 0.5d0) ; sumvz = sumvz + vz(i)
	end do

	if (sumvx /= 0.0d0) then
		do i = 1, n_part
			vx(i) = vx(i) - sumvx/real(n_part)
		end do
	end if

	if (sumvy /= 0.0d0) then
		do i = 1, n_part
			vy(i) = vy(i)  - sumvy/real(n_part)
		end do
	end if

	if (sumvz /= 0.0d0) then
		do i = 1, n_part
			vz(i) = vz(i)  - sumvz/real(n_part)
		end do
	end if

	sumvx = 0; sumvy = 0; sumvz = 0

	Do i = 1, n_part
		sumvx = sumvx + vx(i)
		sumvy = sumvy + vy(i)
		sumvz = sumvz + vz(i)
	End do

	t = 0.0d0

	Print*, sumvx, sumvy, sumvz
	write(94,*), t,sumvx, sumvy, sumvz

	poten1 = 0.0d0
	total_poten = 0.0d0
	fx = 0.0d0
	fy = 0.0d0
	fz = 0.0d0
    force1 = 0.0d0
    distance = 0.0d0
	v_sqr = 0.0d0
	
	Do i= 1,n_part
		v_sqr = v_sqr + vx(i)**2 + vy(i)**2 + vz(i)**2
	End do

	temp = v_sqr/(3*n_part)
	write(91,*) t, temp

	total_poten = 0.0d0

	do i=1,n_part-1
		do j=i+1, n_part
			dist_x=x(j)-x(i);dist_y=y(j)-y(i); dist_z=z(j)-z(i)
			if(abs(dist_x)>lx*0.5d0) then
				if(dist_x>0) then
					dist_x = abs(dist_x) -lx 
				else
					dist_x = lx - abs(dist_x)
				endif
			endif
			if(abs(dist_y)>ly*0.5d0) then
				if(dist_y>0) then
					dist_y = abs(dist_y) -ly 
				else
					dist_y = ly - abs(dist_y)
				endif
			endif
			if(abs(dist_z)>lz*0.5d0) then
				if(dist_z>0) then
					dist_z = abs(dist_z) -lz 
				else
					dist_z = lz - abs(dist_z)
				endif
			endif
			dist= dsqrt((dist_x)**2 + (dist_y)**2+ (dist_z)**2) !r square
			distance(i,j) = dist
			distance(j,i) = dist
			if(dist <= rc) then
				
				force1(i,j) = 4*(12*(sig**12)*1.0d0/(dist**13) - 6*(sig**6)*1.0d0/(dist**7)) - f0
				potenprime= 4*((sig**12)*1.0d0/(dist**12) - (sig**6)*1.0d0/(dist**6)) +f0*dist -vprime0
				total_poten = total_poten + potenprime
				fx(i) = fx(i) - force1(i,j)*(dist_x)/dist ; fx(j) = fx(j) + force1(i,j)*(dist_x)/dist
				fy(i) = fy(i) - force1(i,j)*(dist_y)/dist ; fy(j) = fy(j) + force1(i,j)*(dist_y)/dist
				fz(i) = fz(i) - force1(i,j)*(dist_z)/dist ; fz(j) = fz(j) + force1(i,j)*(dist_z)/dist
			endif		
		enddo
	enddo
	
	
	total_energy = total_poten + (0.5d0*v_sqr)

	write(93,*) t,(total_energy)/dfloat(n_part), (0.5d0*v_sqr)/dfloat(n_part), total_poten/dfloat(n_part)

	do k = 1,n_part
		Write(92,*) fx(k), fy(k), fz(k)
	end do

	Do k = 1,iter

		poten1 = 0.0d0
		force1 = 0.0d0

		t = t+dt
		! Update position
		Do i = 1, n_part

           x(i) =x(i) + vx(i)*dt + 0.5d0*fx(i)*(dt**2)
			y(i) =y(i) + vy(i)*dt + 0.5d0*fy(i)*(dt**2)
			z(i) =z(i) + vz(i)*dt + 0.5d0*fz(i)*(dt**2)
			if(x(i)>=lx) x(i) = x(i)-lx ; if(x(i)<0) x(i) = x(i)+lx; !periodic boundary conditions
			if(y(i)>=ly) y(i) = y(i)-ly ; if(y(i)<0) y(i) = y(i)+ly;
			if(z(i)>=lz) z(i) = z(i)-lz ; if(z(i)<0) z(i) = z(i)+lz;


        End do


        fx2 = 0.0d0
        fy2 = 0.0d0
        fz2 = 0.0d0


		!Force recalculation
		total_poten = 0.0d0
		
       do i=1,n_part-1
		do j=i+1, n_part
			dist_x=x(j)-x(i);dist_y=y(j)-y(i); dist_z=z(j)-z(i)
			if(abs(dist_x)>lx*0.5d0) then
				if(dist_x>0) then
					dist_x = abs(dist_x) -lx 
				else
					dist_x = lx - abs(dist_x)
				endif
			endif
			if(abs(dist_y)>ly*0.5d0) then
				if(dist_y>0) then
					dist_y = abs(dist_y) -ly 
				else
					dist_y = ly - abs(dist_y)
				endif
			endif
			if(abs(dist_z)>lz*0.5d0) then
				if(dist_z>0) then
					dist_z = abs(dist_z) -lz 
				else
					dist_z = lz - abs(dist_z)
				endif
			endif
			dist= dsqrt((dist_x)**2 + (dist_y)**2+ (dist_z)**2) !r square
			if(dist <= rc) then
				
				force1(i,j) = 4*(12*(sig**12)*1.0d0/(dist**13) - 6*(sig**6)*1.0d0/(dist**7)) - f0
				potenprime= 4*((sig**12)*1.0d0/(dist**12) - (sig**6)*1.0d0/(dist**6)) +f0*dist -vprime0
				total_poten = total_poten + potenprime
				fx2(i) = fx2(i) - force1(i,j)*(dist_x)/dist ; fx2(j) = fx2(j) + force1(i,j)*(dist_x)/dist
				fy2(i) = fy2(i) - force1(i,j)*(dist_y)/dist ; fy2(j) = fy2(j) + force1(i,j)*(dist_y)/dist
				fz2(i) = fz2(i) - force1(i,j)*(dist_z)/dist ; fz2(j) = fz2(j) + force1(i,j)*(dist_z)/dist
			endif		
		enddo
	enddo

		sumvx = 0; sumvy = 0; sumvz = 0
		v_sqr = 0
        Do i = 1, n_part !velocity update
            vx2(i) = vx(i) + dt*0.5d0*(fx2(i) + fx(i)) ; sumvx = sumvx + vx2(i)
			vy2(i) = vy(i) + dt*0.5d0*(fy2(i) + fy(i)) ; sumvy = sumvy + vy2(i)
			vz2(i) = vz(i) + dt*0.5d0*(fz2(i) + fz(i)) ; sumvz = sumvz + vz2(i)
			v_sqr = v_sqr + (vx2(i))**2 + (vy2(i))**2 + (vz2(i))**2
        End do
		
		avgv2 = v_sqr/dfloat(n_part)
      	temp = avgv2/3.0d0
		
		!Thermostat

        If (mod(k,500) == 0) then
                vx2 = vx2*(dsqrt(3.0d0/avgv2))
                vy2 = vy2*(dsqrt(3.0d0/avgv2))
                vz2 = vz2*(dsqrt(3.0d0/avgv2))
                v_sqr = 0.0d0
                Do i = 1, n_part
                    v_sqr = v_sqr + vx2(i)**2 + vy2(i)**2 + vz2(i)**2
                End do
                avgv2 = v_sqr/dfloat(n_part)
                temp = avgv2/3.0d0
        End if
        
        fx = fx2
        fy = fy2
        fz = fz2
		vx = vx2
		vy = vy2
		vz = vz2

      
      	write(91,*), t, temp

      	total_energy = (total_poten + (0.5d0*v_sqr))/dfloat(n_part)

		write(93,*) t,total_energy,((0.5d0*v_sqr))/dfloat(n_part), total_poten/dfloat(n_part)  !, (sum(poten1))/dfloat(n_part)  !
      	write(94,*), t,sumvx, sumvy, sumvz
      	
      	! Maxwell Boltzmann statistics
        
		
		
        If (k >= 8000 .and. mod(k,100) == 0) then
            do i = 1,n_part
            	max_boltz(k,i) = dsqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
            end do
            r = 0.0d0
            !do i = 1,n_part
            !	do jj = 1,75
            !	   do j = 1,n_part-1	
           ! 			If (distance(i,j)>r .and. distance(i,j)<r+dr) then
           ! 				n_r(k,jj,i) = n_r(k,jj,i) + 1
            !			End if
            !		End do
            !		r = r+dr
            !	End do
            !End do
        End if
		
					
End do
		do i = 1,n_part
			counter = 0
			do k = 1,iter
			avg_v_max_boltz(i) = avg_v_max_boltz(i) +  (max_boltz(k,i))
			if (abs(max_boltz(k,i)) > 0.0d0) counter = counter + 1
			end do
			avg_v_max_boltz(i) = avg_v_max_boltz(i)/counter
			write(95,*), avg_v_max_boltz(i)
		End do
		
		!!avg_n_r = 0
		!do i = 1,n_part
		!	counter = 0
		!	do jj = 1,75
		!		do k = 1,iter
		!			avg_n_r(i,jj) = avg_n_r(i,jj) + n_r(k,jj,i)
		!			
		!		end do
		!	End do
		!End do
		
	deallocate(x,y,z, vx,vy,vz, x2, vx2)
	deallocate(distance, force1, y2, z2,vy2, vz2)
	deallocate(fx, fy, fz, fx2, fy2, fz2, poten1)
	deallocate(max_boltz, avg_v_max_boltz)
		
		
End Program





!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine addnumtostring(string,number)
implicit none
integer :: i,strlen,number,nodig,num,snum
character*(*) ::  string

      snum=number
      do i=len(string),1,-1
       if (string(i:i) .ne. ' ') goto 10
      enddo
   10 strlen=i


        nodig=int(log10(1.0*snum+0.1))+1
        do i=nodig,1,-1
       num=snum/10**(i-1)
       string(strlen+1:strlen+1)=char(48+num)
       strlen=strlen+1
       snum=snum-num*10**(i-1)
        enddo

        return
end subroutine addnumtostring


















