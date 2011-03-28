! Concentration fields

module fields
use global
implicit none
save

real, allocatable :: S1P_conc(:,:,:), S1P_grad(:,:,:,:)
real, allocatable :: RANKL_conc(:,:,:), RANKL_grad(:,:,:,:), RANKLinflux(:,:,:)

contains

!----------------------------------------------------------------------------------------
! Check to see if (x,y,z) is outside the grid
!----------------------------------------------------------------------------------------
logical function outside(x,y,z)
integer :: x, y, z
outside = .true.
if (x < 1 .or. x > NX) return
if (y < 1 .or. y > NY) return
if (z < 1 .or. z > NZ) return
outside = .false.
end function

!----------------------------------------------------------------------------------------
! Since it is assumed that the cells in the marrow have an insignificant effect on the 
! cytokine concentration, and also that the concentration in the blood is relatively
! constant, it is valid to consider the cytokine concentration field in the marrow
! to be unchanging.  Therefore we can solve for the steady-state concentration field
! once, compute the gradient field from it, and use this for all chemotaxis calculations.
!
! This all works, but there is now serious doubt as to whether there is such a thing as
! S1P chemotaxis.  According to Irina Grigorova, S1P1 on a T cell enables it to cross
! the endothelial boundary into a sinus (high-S1P), but the cell must reach the sinus
! by random motion.
!----------------------------------------------------------------------------------------
subroutine init_fields

if (use_RANK) then
	call init_RANKL
endif
if (S1P_chemotaxis) then
	call init_S1P
endif
end subroutine

!----------------------------------------------------------------------------------------
! The flux of S1P into the gridcell (x,y,z) is proportional to influx(x,y,z)
! The first task is to determine the equilibrium S1P concentration field.
! One way is to solve the diffusion equation over a time long enough to reach steady-state.
! S1P both diffuses and decays.
! Treat the grid spacing as the unit of distance, therefore gridcell volume = 1.
! The flux at the endothelial boundary is a proportionality constant (permeability) Kperm
! times the area (which is influx(x,y,z)) times the concentration drop across the
! endothelium.  Normalize blood concentration to 1.  Decay rate is Kdecay.
! Note that with gridcell volume = 1 the mass of S1P = concentration.
! Note that influx(x,y,z) also flags the region of diffusion of S1P, since 
! influx(x,y,z) >= 0 for gridcells within the region.
!----------------------------------------------------------------------------------------
subroutine init_S1P
integer :: x, y, z, xx, yy, zz, k, nb
real, allocatable :: influx(:,:,:)
real, allocatable :: dCdt(:,:,:)
real, parameter :: dt = 1.0
real, parameter :: Kperm = 0.05		! determines rate of influx of S1P
real :: g(3), gamp, gmax
logical :: steady = .true.

write(logmsg,*) 'Initializing S1P'
call logger(logmsg)
allocate(S1P_conc(NX,NY,NZ))
allocate(S1P_grad(3,NX,NY,NZ))
allocate(influx(NX,NY,NZ))

influx = -1

do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			nb = 0
			if (occupancy(x,y,z)%region /= MARROW) cycle
			do k = 1,6
				xx = x + neumann(1,k)
				yy = y + neumann(2,k)
				zz = z + neumann(3,k)
				if (outside(xx,yy,zz)) cycle
				if (occupancy(xx,yy,zz)%region == BLOOD) then
					nb = nb + 1
				endif
			enddo
			influx(x,y,z) = Kperm*nb
		enddo
	enddo
enddo

if (steady) then
	call steadystate(influx,S1P_KDIFFUSION,S1P_KDECAY,S1P_conc)
else
!	allocate(dCdt(NX,NY,NZ))
!	dCdt = 0
!	do it = 1,nt
!		if (mod(it,100) == 0) write(*,*) 'it: ',it,S1P_conc(25,25,25)
!		do z = 1,NZ
!			do y = 1,NY
!				do x = 1,NX
!					dC = 0
!					if (influx(x,y,z) < 0) cycle
!					C = S1P_conc(x,y,z)
!					dC = -Kdecay*C
!					if (influx(x,y,z) > 0) then
!						dC = dC + influx(x,y,z)*(1 - C)
!					endif
!					sum = 0
!					nb = 0
!					do k = 1,6
!						xx = x + neumann(1,k)
!						yy = y + neumann(2,k)
!						zz = z + neumann(3,k)
!						if (outside(xx,yy,zz)) cycle
!						if (influx(xx,yy,zz) < 0) cycle
!						nb = nb + 1
!						sum = sum + S1P_conc(xx,yy,zz)
!					enddo
!					sum = sum - nb*C
!					dC = dC + Kdiffusion*sum
!					dCdt(x,y,z) = dC
!				enddo
!			enddo
!		enddo
!		S1P_conc = S1P_conc + dt*dCdt
!	enddo
!	deallocate(dCdt)
endif
! Now compute the gradient field.
call gradient(influx,S1P_conc,S1P_grad)
!write(*,*) 'S1P gradient: ',S1P_grad(:,25,25,25)
deallocate(influx)
gmax = 0
do x = 1,NX
	do y = NBY+1,NY
		do z = 1,NZ
			g = S1P_grad(:,x,y,z)
			gamp = sqrt(dot_product(g,g))
			gmax = max(gamp,gmax)
		enddo
	enddo
enddo
write(logmsg,*) 'Max S1P gradient: ',gmax
call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine init_RANKL
integer :: isignal, site(3), x, y, z, it, nt, isig
real :: dt, sum
real :: g(3), gamp, gmax
!real, parameter :: Kdecay = 0.00001, Kdiffusion = 0.001

write(logmsg,*) 'Initializing RANKL: KDECAY: ', RANKL_KDECAY
call logger(logmsg)
allocate(RANKL_conc(NX,NY,NZ))
allocate(RANKL_grad(3,NX,NY,NZ))
allocate(RANKLinflux(NX,NY,NZ))

RANKLinflux = -1
do x = 1,NX
	do y = NBY+1,NY
		do z = 1,NZ
			if (occupancy(x,y,z)%region == MARROW .or. occupancy(x,y,z)%region == LAYER) RANKLinflux(x,y,z) = 0
		enddo
	enddo
enddo
!do isignal = 1,nsignal
!	site = signal(isignal)%site
!	influx(site(1),site(2)+1,site(3)) = RANK_BONE_RATIO*signal(isignal)%intensity
!enddo
y = NBY+1
do x = 1,NX
	do z = 1,NZ
		RANKLinflux(x,y,z) = RANK_BONE_RATIO*surface(x,z)%signal
	enddo
enddo
call steadystate(RANKLinflux,RANKL_KDIFFUSION,RANKL_KDECAY,RANKL_conc)
call gradient(RANKLinflux,RANKL_conc,RANKL_grad)
!write(*,*) 'RANKL gradient: ',RANKL_grad(:,25,25,25)
write(logmsg,*) 0,RANKL_conc(NX/2,NBY+2,NZ/2)
call logger(logmsg)
! Testing convergence
nt = 10
dt = 100*DELTA_T
do it = 1,nt
	call evolve(RANKLinflux,RANKL_KDIFFUSION,RANKL_KDECAY,RANKL_conc,dt)
	write(logmsg,*) it,RANKL_conc(NX/2,NBY+2,NZ/2)
	call logger(logmsg)
	sum = 0
	do isig = 1,nsignal
		site = signal(isig)%site
		sum = sum + RANKL_conc(site(1),NBY+1,site(3))
	enddo
	write(logmsg,*) 'Mean patch RANKL: ',sum/nsignal
	call logger(logmsg)
enddo
gmax = 0
do x = 1,NX
	do y = NBY+1,NY
		do z = 1,NZ
			g = RANKL_grad(:,x,y,z)
			gamp = sqrt(dot_product(g,g))
			gmax = max(gamp,gmax)
		enddo
	enddo
enddo
write(logmsg,*) 'Max RANKL gradient: ',gmax
call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine EvolveRANKL(nt,totsig)
integer :: nt
real :: totsig
!real, allocatable :: influx(:,:,:)
integer :: x, y, z, site(3), isig
real :: dt, sum

dt = nt*DELTA_T
totsig = 0
y = NBY+1
do x = 1,NX
	do z = 1,NZ
		RANKLinflux(x,y,z) = RANK_BONE_RATIO*surface(x,z)%signal
		totsig = totsig + RANKLinflux(x,y,z)
	enddo
enddo
call evolve(RANKLinflux,RANKL_KDIFFUSION,RANKL_KDECAY,RANKL_conc,dt)
call gradient(RANKLinflux,RANKL_conc,RANKL_grad)
sum = 0
do isig = 1,nsignal
	site = signal(isig)%site
	sum = sum + RANKL_conc(site(1),NBY+1,site(3))
enddo
write(logmsg,*) 'Mean patch RANKL: ',sum/nsignal, totsig
call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------------
! Solve for steady-state concentration field for a constituent that is defined by:
!	influx(:,:,:)	if > 0, the rate of influx into grid cell (x,y,z)
!					if < 0, (x,y,z) is not MARROW
!	Kdiffusion		diffusion coefficient (um^2.min^-1)
!	Kdecay			decay coefficient (min^-1)
!----------------------------------------------------------------------------------------
subroutine steadystate(influx,Kdiffusion,Kdecay,C)
real :: influx(:,:,:), C(:,:,:)
real :: Kdiffusion, Kdecay
real :: dx2diff, total, maxchange, dC, sum, dV
real, parameter :: alpha = 0.5
real, parameter :: tol = 1.0e-6		! max change in C at any site as fraction of average C
integer :: x, y, z, xx, yy, zz, nb, nc, k, it
integer :: nt = 10000
real, allocatable :: Ctemp(:,:,:)

dx2diff = DELTA_X**2/Kdiffusion
dV = DELTA_X**3
allocate(Ctemp(NX,NY,NZ))
maxchange = 1.0e10
total = 0
nc = 1
C = 0
do it = 1,nt
	if (mod(it,100) == 0) then
		x = NX/2
		z = NZ/2
!		write(*,'(i6,20e10.3)') it,(C(x,y,z),y=NBY+1,NY/2-NBY,4),maxchange,total/nc
	endif
	if (maxchange < tol*total/nc) then
		write(logmsg,*) 'Convergence reached: it: ',it
		call logger(logmsg)
		exit
	endif
	maxchange = 0
	total = 0
	nc = 0
	do z = 1,NZ
		do y = NBY+1,NY
			do x = 1,NX
				if (influx(x,y,z) < 0) cycle
				if (influx(x,y,z) > 0) then
					dC = influx(x,y,z)*dx2diff
				else
					dC = 0
				endif
				sum = 0
				nb = 0
				do k = 1,6
					xx = x + neumann(1,k)
					yy = y + neumann(2,k)
					zz = z + neumann(3,k)
					if (outside(xx,yy,zz)) cycle
					if (influx(xx,yy,zz) < 0) cycle
					nb = nb + 1
					sum = sum + C(xx,yy,zz)
				enddo
!				Ctemp(x,y,z) = alpha*(sum + dC)/(Kdecay*dx2diff + dC + nb) + (1-alpha)*C(x,y,z)
				Ctemp(x,y,z) = alpha*(DELTA_X*Kdiffusion*sum + influx(x,y,z))/(Kdecay*dV + nb*DELTA_X*Kdiffusion) + (1-alpha)*C(x,y,z)
				dC = abs(Ctemp(x,y,z) - C(x,y,z))
				maxchange = max(dC,maxchange)
				nc = nc + 1
				total = total + Ctemp(x,y,z)
			enddo
		enddo
	enddo
	C = Ctemp
enddo
deallocate(Ctemp)

end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine gradient(influx,C,grad)
real :: influx(:,:,:), C(:,:,:), grad(:,:,:,:)
integer :: x, y, z, xx, yy, zz, x1, x2, y1, y2, z1, z2, i, k
real :: g(3)
logical :: missed
real, parameter :: MISSING_VAL = 1.0e10

grad = 0
do z = 1,NZ
	do y = NBY+1,NY
		do x = 1,NX
			if (influx(x,y,z) < 0) cycle
			x1 = x - 1
			x2 = x + 1
			if (x1 < 1 .or. x2 > NX) then
				g(1) = 0
			elseif (influx(x1,y,z) >= 0 .and. influx(x2,y,z) >= 0) then
				g(1) = (C(x2,y,z) - C(x1,y,z))/(2*DELTA_X)
			else
				g(1) = MISSING_VAL
			endif
			y1 = y - 1
			y2 = y + 1
			if (y1 < NBY .or. y2 > NY) then
				g(2) = 0
			elseif (influx(x,y1,z) >= 0 .and. influx(x,y2,z) >= 0) then
				g(2) = (C(x,y2,z) - C(x,y1,z))/(2*DELTA_X)
			else
				g(2) = MISSING_VAL
			endif
			z1 = z - 1
			z2 = z + 1
			if (z1 < 1 .or. z2 > NZ) then
				g(3) = 0
			elseif (influx(x,y,z1) >= 0 .and. influx(x,y,z2) >= 0) then
				g(3) = (C(x,y,z2) - C(x,y,z1))/(2*DELTA_X)
			else
				g(3) = MISSING_VAL
			endif
			grad(:,x,y,z) = g
		enddo
	enddo
enddo
do z = 1,NZ
	do y = NBY+1,NY
		do x = 1,NX
			if (influx(x,y,z) < 0) cycle
			do i = 1,3
				if (grad(i,x,y,z) == MISSING_VAL) then
					missed = .true.
					grad(i,x,y,z) = 0
					do k = 1,6
						xx = x + neumann(1,k)
						yy = y + neumann(2,k)
						zz = z + neumann(3,k)
						if (outside(xx,yy,zz)) cycle
						if (influx(xx,yy,zz) < 0) cycle
						if (grad(i,xx,yy,zz) /= MISSING_VAL) then
							grad(i,x,y,z) = grad(i,xx,yy,zz)
							missed = .false.
							exit
						endif
					enddo
					if (missed) then
!						write(*,*) 'Missing gradient at: ',x,y,z,grad(:,x,y,z)
					endif
				endif
			enddo
		enddo
	enddo
enddo
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine evolve(influx,Kdiffusion,Kdecay,C,dt)
real :: influx(:,:,:), C(:,:,:)
real :: Kdiffusion, Kdecay, dt
real :: dx2diff, total, maxchange, C0, dC, sum, dV, dMdt
real, parameter :: alpha = 0.99
integer :: x, y, z, xx, yy, zz, nb, nc, k, it
real, allocatable :: Ctemp(:,:,:)

dV = DELTA_X**3
allocate(Ctemp(NX,NY,NZ))
do z = 1,NZ
	do y = NBY+1,NY
		do x = 1,NX
			if (influx(x,y,z) < 0) cycle
			C0 = C(x,y,z)
			sum = 0
			nb = 0
			do k = 1,6
				xx = x + neumann(1,k)
				yy = y + neumann(2,k)
				zz = z + neumann(3,k)
				if (outside(xx,yy,zz)) cycle
				if (influx(xx,yy,zz) < 0) cycle
				nb = nb + 1
				sum = sum + C(xx,yy,zz)
			enddo
			dMdt = Kdiffusion*DELTA_X*(sum - nb*C0) - Kdecay*C0*dV + influx(x,y,z)
			Ctemp(x,y,z) = (C(x,y,z)*dV + dMdt*dt)/dV
!			if (x == NX/2 .and. y == NBY+2 .and. z == NZ/2) write(*,*) 'dMdt: ',dMdt
!			Ctemp(x,y,z) = dt*(alpha*(sum + dC)/(Kdecay*dx2diff + dC + nb) +
		enddo
	enddo
enddo
C = Ctemp
deallocate(Ctemp)

end subroutine
			
end module