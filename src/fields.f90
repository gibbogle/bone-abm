! Concentration fields

module fields
use global
use chemokine

implicit none
save

real, allocatable :: S1P_conc(:,:,:), S1P_grad(:,:,:,:)
real, allocatable :: CXCL12_conc(:,:,:), CXCL12_grad(:,:,:,:), CXCL12_influx(:,:,:)

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
! Convert halflife in hours to a decay rate /min
!----------------------------------------------------------------------------------------
real function DecayRate(halflife)
real :: halflife

DecayRate = log(2.0)/(halflife*60)
end function

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine setupChemo

chemo(CXCL12)%name = "CXCL12"
chemo(CXCL12)%diff_coef = CXCL12_KDIFFUSION
chemo(CXCL12)%decay_rate = DecayRate(CXCL12_HALFLIFE)
chemo(S1P)%name = "S1P"
chemo(S1P)%diff_coef = S1P_KDIFFUSION
chemo(S1P)%decay_rate = DecayRate(S1P_HALFLIFE)
end subroutine

!----------------------------------------------------------------------------------------
! Since it is assumed that the cells in the marrow have an insignificant effect on the 
! cytokine concentration, and also that the concentration in the blood is relatively
! constant, it is valid to consider the cytokine concentration field in the marrow
! to be unchanging.  Therefore we can solve for the steady-state concentration field
! once, compute the gradient field from it, and use this for all chemotaxis calculations.
!
! This all works, but there is now serious doubt as to whether there is such a thing as
! S1P chemotaxis.  According to Irina Grigorova, S1PR1 on a T cell enables it to cross
! the endothelial boundary into a sinus (high-S1P), but the cell must reach the sinus
! by random motion.
!----------------------------------------------------------------------------------------
subroutine init_fields

allocate(CXCL12_conc(NX,NY,NZ))
allocate(CXCL12_grad(3,NX,NY,NZ))
allocate(CXCL12_influx(NX,NY,NZ))
if (use_CXCL12) then
	call init_CXCL12
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
	call steadystate(influx,chemo(S1P)%diff_coef,chemo(S1P)%decay_rate,S1P_conc)
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
subroutine init_CXCL12
integer :: isignal, site(3), x, y, z, it, nt, isig, iblast
real :: dt, sum
real :: g(3), gamp, gmax
type(osteoblast_type), pointer :: pblast
!real, parameter :: Kdecay = 0.00001, Kdiffusion = 0.001

write(logmsg,*) 'Initializing CXCL12: KDECAY: ', chemo(CXCL12)%decay_rate
!call formulate(CXCL12_KDIFFUSION,CXCL12_KDECAY)

call logger(logmsg)
!allocate(CXCL12_conc(NX,NY,NZ))
!allocate(CXCL12_grad(3,NX,NY,NZ))
!allocate(CXCL12_influx(NX,NY,NZ))

CXCL12_grad = 0
CXCL12_influx = -1
do x = 1,NX
	do y = NBY+1,NY
		do z = 1,NZ
			if (occupancy(x,y,z)%region == MARROW .or. occupancy(x,y,z)%region == LAYER) CXCL12_influx(x,y,z) = 0
		enddo
	enddo
enddo
y = NBY+1
!do x = 1,NX
!	do z = 1,NZ
!		CXCL12_influx(x,y,z) = RANK_BONE_RATIO*surface(x,z)%signal
!	enddo
!enddo
sum = 0
do iblast = 1,nblast
	pblast => blast(iblast)
	x = pblast%site(1)
	z = pblast%site(3)
	CXCL12_influx(x,y,z) = BlastSignal(pblast)
	sum = sum + CXCL12_influx(x,y,z)
	write(logmsg,*) iblast,x,z,CXCL12_influx(x,y,z)
	call logger(logmsg)
enddo
if (sum == 0) return

call steadystate(CXCL12_influx,chemo(CXCL12)%diff_coef,chemo(CXCL12)%decay_rate,CXCL12_conc)
call gradient(CXCL12_influx,CXCL12_conc,CXCL12_grad)
!write(*,*) 'CXCL12 gradient: ',CXCL12_grad(:,25,25,25)
write(logmsg,*) 0,CXCL12_conc(NX/2,NBY+2,NZ/2)
call logger(logmsg)
! Testing convergence
nt = 0
dt = 100*DELTA_T
do it = 1,nt
	call evolve(CXCL12_influx,chemo(CXCL12)%diff_coef,chemo(CXCL12)%decay_rate,CXCL12_conc,dt)
	write(logmsg,*) it,CXCL12_conc(NX/2,NBY+2,NZ/2)
	call logger(logmsg)
	sum = 0
	do isig = 1,nsignal
		site = signal(isig)%site
		sum = sum + CXCL12_conc(site(1),NBY+1,site(3))
	enddo
	write(logmsg,*) 'Mean patch CXCL12: ',sum/nsignal
	call logger(logmsg)
enddo
y = NBY + 3
gmax = 0
do x = 1,NX
!	do y = NBY+1,NY
		do z = 1,NZ
			g = CXCL12_grad(:,x,y,z)
			gamp = sqrt(dot_product(g,g))
			gmax = max(gamp,gmax)
		enddo
!	enddo
enddo
write(logmsg,*) 'Max CXCL12 gradient at y=NBY+3: ',gmax
call logger(logmsg)
CXCL12_GRADLIM = gmax
CXCL12_initialized = .true.
end subroutine

!----------------------------------------------------------------------------------------
! Need to suppress osteocyte signal from surface sites covered by OCs 
!----------------------------------------------------------------------------------------
subroutine EvolveCXCL12(nt,totsig)
integer :: nt
real :: totsig
!real, allocatable :: influx(:,:,:)
!real, allocatable :: factor(:,:)
integer :: x, y, z, site(3), isig, i, iblast
real :: dt, sum, sig
!type(osteoclast_type), pointer :: pclast
type(osteoblast_type), pointer :: pblast

dt = nt*DELTA_T
!allocate(factor(NX,NZ))
!factor = 1
!do iclast = 1,nclast
!	pclast => clast(iclast)
!	if (pclast%status == DEAD) cycle
!	do i = 1,pclast%npit
!		site = pclast%site + pclast%pit(i)%delta
!		factor(site(1),site(3)) = 0
!	enddo
!enddo
totsig = 0
y = NBY+1
!do x = 1,NX
!	do z = 1,NZ
!		CXCL12_influx(x,y,z) = RANK_BONE_RATIO*factor(x,z)*surface(x,z)%signal
!		totsig = totsig + CXCL12_influx(x,y,z)
!	enddo
!enddo
CXCL12_influx(:,y,:) = 0
do iblast = 1,nblast
	pblast => blast(iblast)
	if (pblast%status /= ALIVE) cycle
	x = pblast%site(1)
	z = pblast%site(3)
	sig = BlastSignal(pblast)
	write(*,'(a,4i6,f8.3)') 'sig: ',iblast,nblast,x,z,sig
	if (sig > OB_SIGNAL_THRESHOLD) then
		CXCL12_influx(x,y,z) = sig
	endif
	totsig = totsig + CXCL12_influx(x,y,z)
enddo
write(logmsg,*) 'evolveCXCL12: totsig: ',istep,totsig
call logger(logmsg)
!if (totsig > 0) then
	call evolve(CXCL12_influx,chemo(CXCL12)%diff_coef,chemo(CXCL12)%decay_rate,CXCL12_conc,dt)
	call gradient(CXCL12_influx,CXCL12_conc,CXCL12_grad)
!endif
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine steadystate(influx,Kdiffusion,Kdecay,C)
real :: influx(:,:,:), C(:,:,:)
real :: Kdiffusion, Kdecay
real :: dx2diff, total, maxchange, maxchange_par, total_par
real, parameter :: alpha = 0.5
real, parameter :: tol = 1.0e-5		! max change in C at any site as fraction of average C
integer :: nc, nc_par, k, it, kpar, dz
integer :: zlim(2,16), z1, z2, zfr, zto, n	! max 16 threads
real, allocatable :: C_par(:,:,:)
integer :: nt = 10000
real, allocatable :: Ctemp(:,:,:)

dz = NZ/Mnodes + 0.5
do k = 1,Mnodes
	zlim(1,k) = (k-1)*dz + 1
	zlim(2,k) = k*dz
enddo
zlim(2,Mnodes) = NZ
C = 0
do it = 1,nt
	maxchange = 0
	total = 0
	nc = 0
	!$omp parallel do private(z1,z2,n,C_par,maxchange_par,total_par,nc_par)
	do kpar = 0,Mnodes-1
		z1 = zlim(1,kpar+1)
		z2 = zlim(2,kpar+1)
		n = z2 - z1 + 1
		allocate(C_par(NX,NY,n))
		call par_steadystate(C,Kdiffusion,Kdecay,C_par,influx,z1,z2,maxchange_par,total_par,nc_par,kpar)
		C(:,:,z1:z2) = C_par(:,:,1:n)
		deallocate(C_par)
		nc = nc + nc_par
		total = total + total_par
		maxchange = max(maxchange,maxchange_par)
	enddo
	if (maxchange < tol*total/nc) then
		write(logmsg,*) 'Convergence reached: it: ',it
		call logger(logmsg)
		exit
	endif
enddo
end subroutine	

!----------------------------------------------------------------------------------------
! Note that the array Ctemp() is defined with z:z1-z2, i.e. with
! a range of n = z2 - z1 + 1 values.  The read accesses of C() and influx() are shared 
! in this version.  Does this have a big penalty?
!----------------------------------------------------------------------------------------
subroutine par_steadystate(C,Kdiffusion,Kdecay,Ctemp,influx,z1,z2,maxchange,total,nc,kpar)
real :: C(:,:,:), Ctemp(:,:,:), influx(:,:,:)
integer :: z1, z2, kpar
real :: Kdiffusion, Kdecay
real :: dx2diff, total, maxchange, dC, sum, dV
real, parameter :: alpha = 0.7
integer :: x, y, z, xx, yy, zz, nb, nc, k, zpar

dx2diff = DELTA_X**2/Kdiffusion
dV = DELTA_X**3

maxchange = 0
total = 0
nc = 0
do zpar = 1,z2-z1+1
	z = zpar + z1-1
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
!			Ctemp(x,y,z) = alpha*(sum + dC)/(Kdecay*dx2diff + dC + nb) + (1-alpha)*C(x,y,z)
			Ctemp(x,y,zpar) = alpha*(DELTA_X*Kdiffusion*sum + influx(x,y,z))/(Kdecay*dV + nb*DELTA_X*Kdiffusion) + (1-alpha)*C(x,y,z)
			dC = abs(Ctemp(x,y,zpar) - C(x,y,z))
			maxchange = max(dC,maxchange)
			nc = nc + 1
			total = total + Ctemp(x,y,zpar)
		enddo
	enddo
enddo

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
			if (isnan(g(1)) .or. isnan(g(2)) .or. isnan(g(3))) then
				write(logmsg,'(a,3i6,3f8.4)') 'gradient: ',x,y,z,g
				call logger(logmsg)
				write(logmsg,'(a,5i6,2f8.4)') 'x: ',x,y,z,x1,x2,C(x1,y,z),C(x2,y,z)
				call logger(logmsg)
				write(logmsg,'(a,5i6,2f8.4)') 'y: ',x,y,z,y1,y2,C(x,y1,z),C(x,y2,z)
				call logger(logmsg)
				write(logmsg,'(a,5i6,2f8.4)') 'z: ',x,y,z,z1,z2,C(x,y,z1),C(x,y,z2)
				call logger(logmsg)
				stop
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
! This is a crude explicit solver.  If dt is too large it goes unstable.
! Should be replaced by a higher order solver, e.g. Runge-Kutta-Chebychev
!----------------------------------------------------------------------------------------
subroutine evolve(influx,Kdiffusion,Kdecay,C,dt)
real :: influx(:,:,:), C(:,:,:)
real :: Kdiffusion, Kdecay, dt
real :: dx2diff, total, maxchange, C0, dC, sum, dV, dMdt, Cmin, Cmax
real, parameter :: alpha = 0.99
integer :: x, y, z, xx, yy, zz, nb, nc, k, it
real, allocatable :: Ctemp(:,:,:)

dV = DELTA_X**3
allocate(Ctemp(NX,NY,NZ))
Cmin = 1.0e20
Cmax = -Cmin
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
			Cmin = min(Cmin,Ctemp(x,y,z))
			Cmax = max(Cmax,Ctemp(x,y,z))
!			if (x == NX/2 .and. y == NBY+2 .and. z == NZ/2) write(*,*) 'dMdt: ',dMdt
!			Ctemp(x,y,z) = dt*(alpha*(sum + dC)/(Kdecay*dx2diff + dC + nb) +
		enddo
	enddo
enddo
C = Ctemp
deallocate(Ctemp)
!write(logmsg,'(a,2e12.3)') 'Cmin, Cmax: ',Cmin,Cmax
!call logger(logmsg)
end subroutine


end module