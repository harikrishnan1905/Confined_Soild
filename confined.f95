use omp_lib
implicit none
integer :: m,t,tmax,i,j,iseed,n,tequl,mtot,nwall,m_start,m_end
real (kind = 8) :: gasdev,ran2,vp,dm,a
real (kind = 8), dimension (:) , allocatable :: xx,xy,ax,ay,the
real (kind =8) :: L,dt,temp,rho,d,pe
external :: gasdev,ran2
open (400,file='test1.txt') 
!open (100,file='density150_second.txt') 

!####################################################################################################
!************************************** PARAMETERS **************************************************
!####################################################################################################

dt = 0.00001D0  
n = 64				
m = n*n 						! total number of particles
dm = 1.12246205						! diameter of the particle
Pe = 150.00d0
iseed = -887  
temp = 1.0d0 
vp = Pe  						! self propulsion velocity
d = 1.0!temp/gam	      				! translational diffusion coefficient
tmax = 5
tequl = 1000000
rho = 0.80d0						! density
l = sqrt(m/rho) 				! simulation box is square with side l
a = l/n							! a specifies the base of the triangle	
nwall = int(l/dm) + 1
m_start = nwall + 1
m_end = nwall + m
mtot = m + 2*nwall					! total number of particles including the barrier ones
print*, l,a,nwall,mtot,m_start,m_end
allocate(xx(mtot),xy(mtot),the(mtot))				! 'the(i)' specifies the direction of the i-th particle	
allocate(ax(mtot),ay(mtot))
the = 0.D0 ; ax = 0.D0 ; ay = 0.D0 
xx = 0.D0 ; xy = 0.D0 

!####################################################################################################
!************************************** MAIN ********************************************************
!####################################################################################################

call initial(xx,xy,the,m,mtot,n,l,dt,iseed,dm)	
do t = 1,tmax
	call force(ax,ay,xx,xy,mtot,L,m_start,m_end)	
	call integrate(ax,ay,xx,xy,the,m,mtot,dt,L,iseed,d,vp,m_start,m_end)	
	if(t .eq. tmax) then	
		do j = 1,mtot
			write(400,*) xx(j),xy(j)
		end do
	end if	
end do	
end

!####################################################################################################
!**************************************** INITIALIZATION ********************************************
!####################################################################################################

subroutine initial(xx,xy,the,m,mtot,n,l,dt,iseed,d)	
implicit none
integer ,intent(in) :: iseed,m,n,mtot
integer :: i,j,k,p,g,mw,m_start,m_end,nwall
real (kind = 8),intent(inout) :: xx(mtot),xy(mtot),the(mtot)
real (kind = 8)::gasdev,a,ran2
real (kind = 8), intent (in) :: dt,l,d
REAL (kind = 8), PARAMETER :: Pi = 3.14159265
external :: gasdev,ran2					
nwall = (mtot-m)/2
a = l/n
! lower barrier
g = 1
do i = 0,nwall-1
	xx(g) = i*d
	xy(g) = 0.0d0
	the(g) = 0.0d0	
	g = g + 1
end do
! m particle in triangular lattice in lxl box 
g = nwall+1	
mw = m + nwall		
do i =0,m-1						
	if (g.le.mw) then
		xx(g) = a*i				
		xy(g) = 2.0d0
		p = int((i)/n)
		k = mod(p,2)
		if (i.gt.(n-1)) then			
			xx(g) = xx(g-n) + a*0.5D0	
			if (k.eq. 0) then		
				xx(g) = xx(g) - a*1.0d0
			end if
			xy(g) = xy(g-n) + a		
		end if
	end if	
	g = g + 1
end do
g = nwall + 1
do i=1,m
	the(g) = 2.0d0*Pi*ran2(iseed)	
	g = g + 1
end do
! upper barrier
g = m + nwall + 1
xx(m + nwall + 1) = 0.0d0
xy(m + nwall + 1) = xy(m+nwall+1-n) + 2.0d0
do i = 0, nwall-1
	xx(g) = xx(m + nwall + 1) + i*d
	xy(g) = xy(m+nwall+1) 
	the(g) = 0.0d0	
	g = g + 1 
end do
do i=1,mtot
	write(50,*) xx(i),xy(i)
	write(60,*) the(i)
end do
end subroutine initial

!####################################################################################################
!************************************** FORCE CALCULATION ******************************************* 
!####################################################################################################

subroutine force(ax,ay,xx,xy,mtot,L,m_start,m_end)
implicit none
integer , intent (in) :: mtot,m_start,m_end
integer :: i,j
real (kind = 8), intent (inout) :: ax(mtot),ay(mtot)
real (kind = 8), intent (in) :: xx(mtot),xy(mtot),L
real (kind = 8) :: rsqd,ff,r2i,r6i,rcut,rcut2,virij,rx,ry
rx = 0.D0 ; ry = 0.D0 
ax = 0.D0 ; ay = 0.D0 
rcut = 1.12246205 !2**(1.0d0/6.0d0)			! cut off distance for the WCA potential is 2^(1/6)
rcut2 = rcut*rcut !1.259921054
call omp_set_num_threads(10)
!$omp parallel do private(rx,ry,rsqd,r2i,r6i,virij,ff,j) schedule(dynamic) 
do i=1,mtot-1
	do j=i+1,mtot
		rsqd = 0.D0
		rx = xx(i) - xx(j)
		rx = rx - l*anint(rx/l)
		ry = xy(i) - xy(j)	
		!if (((i.ge.m_end).and.(i.le.mtot)).or.((j.ge.m_end).and.(j.le.mtot)).or.((i.le.(m_start-1))).or.((j.le.(m_start-1)))) then
		!	ry = ry 
		!else 
		!	ry = ry - l*anint(ry/l)					
		!end if
		rsqd = rsqd + rx*rx +ry*ry 
	if (rsqd.le.rcut2) then;
		r2i = 1/rsqd
		r6i = r2i**3
		virij = 48*(r6i*r6i-0.5D0*r6i) 		
	        ff = virij*r2i	
	        if (ff .gt. 10000) then	
			print*, 'I BETRAYED YOU !!!!!!!! ', i,j
			EXIT
		end if
		!$omp atomic
		ax(i) = ax(i) + rx*ff		
		!$omp atomic				
		ay(i) = ay(i) + ry*ff	
		!$omp atomic									
		ax(j) = ax(j) - rx*ff
		!$omp atomic
		ay(j) = ay(j) - ry*ff					
	end if			
	end do	
end do
!$omp end parallel do
do i = 1,m_start-1
	ax(i) = 0.0d0
	ay(i) = 0.0d0
end do
do i = m_end + 1 , mtot
	ax(i) = 0.0d0
	ay(i) = 0.0d0
end do
end subroutine force

!####################################################################################################
!************************************** INTEGRATION *************************************************
!####################################################################################################

subroutine integrate(ax,ay,xx,xy,the,m,mtot,dt,L,iseed,d,vp,m_start,m_end)
implicit none
integer , intent (in) :: m,iseed,mtot,m_start,m_end
integer :: i,j,ini,fin
real (kind = 8), intent (in) :: ax(mtot),ay(mtot),dt,L,d,vp  ! d and dr translational and rotational diffusion coefficients
real (kind = 8), intent (inout) :: xx(mtot),xy(mtot),the(mtot)
real (kind = 8) :: g1,g2,g3,gasdev,dr
external :: gasdev
dr = 3.0d0*d						! in the low Reynolds number regime
ini = m_start 
fin = m_end
do i=1,mtot
	if((i .ge. ini).and.(i .le. fin)) then	
		g1 = gasdev(iseed)				
		xx(i) = xx(i) + dt*(vp*cos(the(i))+ ax(i))*d + sqrt(2*d*dt)*g1	
		g2 = gasdev(iseed)				
		xy(i) = xy(i) + dt*(vp*sin(the(i))+ ay(i))*d + sqrt(2*d*dt)*g2        
		g3 = gasdev(iseed)
		the(i) = the(i) + sqrt(2*dr*dt)*g3			
		xx(i) = xx(i) - L*anint((xx(i)/L)-(0.5D0))				
		!xy(i) = xy(i) - L*anint((xy(i)/L)-(0.5D0))
	else
		xx(i) = xx(i)		
	end if	
end do
end subroutine integrate
