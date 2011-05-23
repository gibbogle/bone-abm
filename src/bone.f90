! There is an important design issue in the treatment of signals.
! Currently no distinction is made between signal from the bone surface (osteocytes)
! and signal from osteoblasts (RANKL).
! One very crude approach:
! (a) There is an osteocyte-generated signal from the bone that varies with the depth of bone to be excavated,
!     i.e. initially the signal is strongest at the centre of the remodelling patch, reducing to zero at the edge.
!     As the OCs remove bone, the signal is reduced, until finally it is zero everywhere.
! (b) OB-generated RANKL signal at a grid site is taken to be proportional to the bone signal.  Initially the bone 
!     signal is assumed to have been sustained for some time, and steady-state RANKL field is created.  As bone
!     is removed, the bone signal reduces, and correspondingly RANKL secretion decreases, finally to zero.
! (c) The bone signal exerts a chemotactic influence on the OCs, ensuring that bone removal is mainly restricted to
!     the remodelling patch.

module bone_mod
!!DEC$ ATTRIBUTES DLLEXPORT :: BONE_MOD
use, intrinsic :: ISO_C_binding
use global
use fields
use motion
use omp_lib

implicit none 
save

integer(c_int),BIND(C) :: success = 12345
!DEC$ ATTRIBUTES DLLEXPORT :: success

contains


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine rng_initialisation()
integer, allocatable :: zig_seed(:)
integer :: i
integer :: npar, grainsize = 32

npar = Mnodes
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
deallocate(zig_seed)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine omp_initialisation(ok)
logical :: ok
integer :: npr, nth

ok = .true.
if (Mnodes == 1) return
!!DEC$ IF ( DEFINED (_OPENMP) .OR. DEFINED (IBM))
#if defined(OPENMP) || defined(_OPENMP)
write(logmsg,'(a,i2)') 'Requested Mnodes: ',Mnodes
call logger(logmsg)
!close(nflog)
!ok = .false.
!return
npr = omp_get_num_procs()
write(logmsg,'(a,i2)') 'Machine processors: ',npr
call logger(logmsg)

nth = omp_get_max_threads()
write(logmsg,'(a,i2)') 'Max threads available: ',nth
call logger(logmsg)
if (nth < Mnodes) then
    Mnodes = nth
    write(logmsg,'(a,i2)') 'Setting Mnodes = max thread count: ',nth
	call logger(logmsg)
endif

call omp_set_num_threads(Mnodes)
!$omp parallel
nth = omp_get_num_threads()
write(logmsg,*) 'Threads, max: ',nth,omp_get_max_threads()
call logger(logmsg)
!$omp end parallel
#endif
!!DEC$ END IF
!write(*,*) 'Threads: ',nth
!if (nth > npr) then
!    write(*,*) 'Number of threads exceeds CPUs'
!    stop
!endif

!stop

!call set_affinity
!write(*,*) 'set affinity mappings'

!if (track_DCvisits) then
!    write(*,*) 'To track DC visits %DClist(:) must be allocated for cells'
!    stop
!endif
call logger('did omp_initialisation')
end subroutine


!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine setup(infile,ok)
logical :: ok
character*(*) :: infile
integer :: x, y, z, del, dx, dy, dz, ix, i, site(3), ndiv, icap, k, na, nbr(3)
integer, allocatable :: templist(:,:)
real :: xc, yc, zc, fac, x2, y2, z2, v(3)
integer :: kpar = 0
real(8) :: R
real :: alfa

ok = .true.
inputfile = infile
call read_inputfile(ok)
if (.not.ok) return

Mnodes = ncpu
#if defined(OPENMP) || defined(_OPENMP)
    call logger("OPENMP defined")
    call omp_initialisation(ok)
    if (.not.ok) return
#else
    call logger("OPENMP NOT defined")
    if (Mnodes > 1) then
        write(logmsg,'(a)') 'No OpenMP, using one thread only'
        call logger(logmsg)
        Mnodes = 1
    endif
#endif

call rng_initialisation
PI = 4*atan(1.0)
call make_reldir
call MakeClumpOrder

write(logmsg,*) 'NX,NY,NZ: ',NX,NY,NZ
call logger(logmsg)
allocate(occupancy(NX,NY,NZ))
occupancy%region = MARROW
occupancy%indx = 0
!occupancy%signal = 0
!occupancy%intensity = 0
occupancy%bone_fraction = 0
do y = 1,NBY
	occupancy(:,y,:)%region = BONE
	occupancy(:,y,:)%bone_fraction = 1.0
enddo
occupancy(:,NBY+1,:)%region = LAYER
allocate(surface(NX,NZ))
surface%signal = 0
!MAX_PIT_DEPTH = MAX_PIT_DEPTH/DELTA_X

! Create a test capillary
ncap = 1
allocate(capillary(ncap))
capillary(1)%radius = CAPILLARY_DIAMETER/2
capillary(1)%pos1 = (/ 0.5, NY/2., NZ/2. /)
capillary(1)%pos2 = (/ NX+0.5, NY/2., NZ/2. /)
v = capillary(1)%pos2 - capillary(1)%pos1
capillary(1)%length = rnorm(v)
capillary(1)%surface_area = PI*capillary(1)%radius**2*capillary(1)%length

! Set up capillary sites - this is currently VERY CRUDE
! If the centre of a site (cube) falls inside the capillary tube,
! the site is tagged BLOOD.
! Note that (xc,yc,zc), which lies on the capillary centreline,
! is in (x,y,z) coord values, offset from site index values by 0.5
! E.g. if site(:) = (/2, 3, 4/), corresponding (x,y,z) = (/1.5, 2.5, 3.5/)
! This has the origin of the axis system at the lower left rear corner
! of the region.  (Note that in OpenGL/VTK Z axis is out of the screen.)

ndiv = 100
do icap = 1,ncap
	del = capillary(icap)%radius + 2
	do k = 1,ndiv
		alfa = (k-1.0)/(ndiv-1.0)
		xc = (1-alfa)*capillary(icap)%pos1(1) + alfa*capillary(icap)%pos2(1)
		yc = (1-alfa)*capillary(icap)%pos1(2) + alfa*capillary(icap)%pos2(2)
		zc = (1-alfa)*capillary(icap)%pos1(3) + alfa*capillary(icap)%pos2(3)
		do dx = -del,del
			x = xc + dx
			if (x < 1 .or. x > NX) cycle
			do dy = -del,del
				y = yc + dy
				if (y <= NBY .or. y > NY) cycle
				do dz = -del,del
					z = zc + dz
					if (z < 1 .or. z > NZ) cycle
!					x2 = (x - 0.5 - xc)*(x - 0.5 - xc)
!					y2 = (y - 0.5 - yc)*(y - 0.5 - yc)
!					z2 = (z - 0.5 - zc)*(z - 0.5 - zc)
					x2 = (x - xc)*(x - xc)
					y2 = (y - yc)*(y - yc)
					z2 = (z - zc)*(z - zc)
					! These are the squared distances, in three axis directions, from the
					! site (cube) midpoint to the capillary centreline location (xc,yc,zc)
					if (x2+y2+z2 <= (capillary(icap)%radius)**2) then
						occupancy(x,y,z)%region = BLOOD
!						if (k <= 10) then
!							write(nflog,'(4f6.1,3i4)') xc,yc,zc,sqrt(x2+y2+z2),x,y,z
!						endif
					endif
				enddo
			enddo
		enddo
	enddo
enddo
! Now set up a list of all marrow sites that are adjacent to a blood site
allocate(templist(3,NX*NZ*10))
na = 0
do x = 2,NX-1
	do y = NBY+1,NY-1
		do z = 2,NZ-1
			if (occupancy(x,y,z)%region /= MARROW) cycle
			site = (/x,y,z/)
			do k = 1,6
				nbr = site + neumann(:,k)
				if (occupancy(nbr(1),nbr(2),nbr(3))%region == BLOOD) then
					na = na + 1
					templist(:,na) = site
					exit
				endif
			enddo
		enddo
	enddo
enddo
allocate(entrysite(3,na))
entrysite(:,1:na) = templist(:,1:na)
deallocate(templist)
nentrysites = na

call CreatePatch

NMONO_INITIAL = (NX*(NY-NBY)*NZ*DELTA_X**3/1.0e9)*MONO_PER_MM3	! domain as fraction of 1 mm3 x rate of monocytes
!NSTEM = (PI*NX*CAPILLARY_DIAMETER*DELTA_X**2/1.0e6)*STEM_PER_MM2	! capillary surface area as fraction of 1 mm2 x rate of stem cells
NSTEM = (NX*(NY-NBY)*NZ*DELTA_X**3/1.0e9)*STEM_PER_MM3	! domain as fraction of 1 mm3 x rate of monocytes
NBLAST_INITIAL = (patch%volume*DELTA_X**3)*OB_PER_UM3
!write(*,*) 'Lacuna volume: ',patch%volume,NBLAST_INITIAL
!write(logmsg,*) 'NSTEM, NMONO_INITIAL: ',NSTEM,NMONO_INITIAL
call logger(logmsg)

allocate(mono(MAX_MONO))
allocate(clast(MAX_CLAST))
allocate(blast(MAX_BLAST))
call InitialBlastPlacement
call logger('did InitialBlastPlacement')
!call Initiation

RANKSIGNAL_decayrate = log(2.0)/(RANKSIGNAL_HALFLIFE)    ! rate/min
CXCL12_KDECAY = log(2.0)/(CXCL12_HALFLIFE)    ! rate/min
call init_fields
call logger('Did init_fields')

nclast = 0
nliveclast = 0
nmono = 0
nborn = 0
nleft = 0
mono_cnt = 0
do while (nmono < NMONO_INITIAL)
	x = random_int(1,NX,kpar)
	y = random_int(NBY+1,NY,kpar)
	z = random_int(1,NZ,kpar)
	if (occupancy(x,y,z)%region /= MARROW) cycle
	site = (/x,y,z/)
	call addMono('initial',site)
enddo

allocate(stem(NSTEM))
i = 0
do while (i < NSTEM)
	x = random_int(1,NX,kpar)
	y = random_int(NBY+1,NY,kpar)
	z = random_int(1,NZ,kpar)
	if (occupancy(x,y,z)%region /= MARROW) cycle
	i = i + 1
	stem(i)%ID = i
	stem(i)%site = (/x,y,z/)
!	call random_number(R)
	R = par_uni(kpar)
	stem(i)%dividetime = R*STEM_CYCLETIME	! stem cells are due to divide at random times
	occupancy(x,y,z)%species = STEMCELL
enddo

nclump = 0
stuck = .false.
initiated = .false.
!call setSignal(1,ON,ok)
!if (.not.ok) return

end subroutine

!------------------------------------------------------------------------------------------------
! An elliptical patch of bone surface is identified by osteocyte-generated bone signal.
! Returns the total volume (in grid cells) of the expected lacuna. 
!------------------------------------------------------------------------------------------------
subroutine CreatePatch
real :: vol
integer :: x, z
real :: c

call logger('CreatePatch:')
patch%a = LACUNA_A
patch%b = LACUNA_B
if (HALF_ELLIPSE) then
	patch%x0 = NX/10
else
	patch%x0 = NX/2
endif
patch%z0 = NZ/2
nsignal = 0
vol = 0
do x = 1,NX
	do z = 1,NZ
		surface(x,z)%signal = 0
		surface(x,z)%depth = 0
		surface(x,z)%target_depth = 0
		surface(x,z)%seal = 1
		if (HALF_ELLIPSE .and. x < patch%x0) cycle
		c = ((x-patch%x0)/patch%a)**2 + ((z-patch%z0)/patch%b)**2
		if (c < 1) then
			surface(x,z)%target_depth = (1 - c/2)*MAX_PIT_DEPTH
			surface(x,z)%signal = GetSignal(x,z)
			nsignal = nsignal + 1
			signal(nsignal)%site = (/ x,NBY,z /)
!			write(logmsg,'(2i5,f8.4)') x,z,surface(x,z)%target_depth
!			call logger(logmsg)
			vol = vol + surface(x,z)%target_depth
		endif
	enddo
enddo
vol = vol*1000
patch%volume = vol
write(logmsg,*) 'Patch volume: ',vol
call logger(logmsg)
call logger('did CreatePatch:')
end subroutine

!------------------------------------------------------------------------------------------------
! The initial placement of OBs is determined by: 
!  attraction to signal, 
!  mutual repulsion,
!  randomness.
!------------------------------------------------------------------------------------------------
subroutine InitialBlastPlacement
integer :: xmin, xmax, zmin, zmax, x, z, ib, site(3), dx
real(8) :: R
real :: sig, d2, dfactor, prob, prob0 = 0.5
integer :: kpar=0
real :: d2min = 10	!0.4*OB_SIGNAL_RADIUS**2
real, parameter :: MINIMUM_SIGNAL = 0.7
logical, parameter :: REGULAR = .false.

!NBLAST_INITIAL = 1
write(logmsg,*) 'InitialBlastPlacement: ',NBLAST_INITIAL
call logger(logmsg)
if (HALF_ELLIPSE) then
	xmin = patch%x0
else
	xmin = patch%x0 - patch%a
endif
xmax = patch%x0 + patch%a
zmin = patch%z0 - patch%b
zmax = patch%z0 + patch%b
nblast = 0
if (NBLAST_INITIAL == 1) then
	site = (/ patch%x0, NBY+1, patch%z0 /)
	call CreateOsteoblast(site)
	return
endif
if (REGULAR) then
	dx = 5
	do ib = 1,NBLAST_INITIAL
		x = xmin + 1 + (ib-1)*dx
		z = patch%z0
		sig = surface(x,z)%signal
		if (sig < MINIMUM_SIGNAL) exit
		site = (/ x, NBY+1, z /)
		call CreateOsteoblast(site)
		write(logmsg,'(a,4i4)') 'Placed OB: ',nblast,site
		call logger(logmsg)
	enddo
else
	do
		R = par_uni(kpar)
		x = xmin + R*(xmax-xmin) + 0.5
		R = par_uni(kpar)
		z = zmin + R*(zmax-zmin) + 0.5
		sig = surface(x,z)%signal
		if (sig < MINIMUM_SIGNAL) cycle
		dfactor = 0
		do ib = 1,nblast
			d2 = (x-blast(ib)%site(1))**2 + (z-blast(ib)%site(3))**2
			if (d2 == 0) then
				dfactor = 1.0
				exit
			else
				dfactor = dfactor + d2min/d2
			endif
		enddo
		if (dfactor >= 1) cycle
		prob = prob0*(sig**2)*(1-dfactor)
		R = par_uni(kpar)
		if (R < prob) then
			site = (/ x, NBY+1, z /)
			call CreateOsteoblast(site)
			write(logmsg,'(a,4i4)') 'Placed OB: ',nblast,site
			call logger(logmsg)
			if (nblast == NBLAST_INITIAL) exit
		endif
	enddo
endif
write(logmsg,*) 'did InitialBlastPlacement'
call logger(logmsg)
end subroutine

!------------------------------------------------------------------------------------------------
! To initiate reformation, the surface seal near one OB is removed.
! This could be:
! the highest signal OB
! the leftmost OB
!------------------------------------------------------------------------------------------------
subroutine Initiation
integer :: x, z, xmin, iblast, ibinit, bsite(3), dx, dz
real :: d2, r2, sig, sigmax
logical, parameter :: HIGHEST = .false.

call logger('Initiation')
r2 = OB_SIGNAL_RADIUS**2
!xmin = 1.0e10
!do iblast = 1,nblast
!	bsite = blast(iblast)%site
!	if (bsite(1) < xmin) then
!		xmin = bsite(1)
!		ibinit = iblast
!	endif
!enddo
if (HIGHEST) then
	sigmax = 0
	do iblast = 1,nblast
		bsite = blast(iblast)%site
		sig = surface(bsite(1),bsite(3))%signal
		if (sig > sigmax) then
			sigmax = sig
			ibinit = iblast
		endif
	enddo
else
	xmin = 1000
	do iblast = 1,nblast
		bsite = blast(iblast)%site
		if (bsite(1) < xmin) then
			xmin = bsite(1)
			ibinit = iblast
		endif
	enddo
endif
bsite = blast(ibinit)%site
do dx = -OB_SIGNAL_RADIUS,OB_SIGNAL_RADIUS
	x = bsite(1) + dx
	if (x < 1) cycle
	do dz = -OB_SIGNAL_RADIUS,OB_SIGNAL_RADIUS
		z = bsite(3) + dz
		d2 = dx**2 + dz**2
		if (d2 > r2) cycle
		if (surface(x,z)%target_depth == 0) cycle
		surface(x,z)%seal = 0
	enddo
enddo
!blast%status = DORMANT
blast(ibinit)%status = ALIVE
CXCL12_chemotaxis = .true.
initiated = .true.
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine CreateOsteoblast(site)
integer :: site(3)

nblast = nblast + 1
if (nblast > MAX_BLAST) then
	call logger('CreateOsteoblast: array dimension exceeded')
	stop
endif
blast(nblast)%ID = nblast
blast(nblast)%status = DORMANT
blast(nblast)%site = site
blast(nblast)%iclump = 0
!write(*,*) 'blast: ',nblast,blast(nblast)%site
end subroutine

!------------------------------------------------------------------------------------------------
! Determine rates of resorption at the various sites (pits) under an osteoclast.
! Modified to use only lateral (x,y) distance.  Pit site now an offset from cm.
! Note: cheating by setting rate = 0 when signal = 0
!------------------------------------------------------------------------------------------------
subroutine pitrates(pclast)
type(osteoclast_type), pointer :: pclast
integer :: ipit, site(3), x, z
real :: d, v(3)

do ipit = 1,pclast%npit
	v = pclast%pit(ipit)%delta
	x = pclast%site(1) + v(1)
	z = pclast%site(3) + v(3)
	if (surface(x,z)%signal > 0) then
!		v(2) = 0
!		d = sqrt(dot_product(v,v))
		pclast%pit(ipit)%rate = pclast%pit(ipit)%fraction*resorptionRate(pclast%count,pclast%npit)
	else
		pclast%pit(ipit)%rate = 0
	endif
enddo
end subroutine

!------------------------------------------------------------------------------------------------
! For now the rate of entry of osteoclast-precursor monocytes is an input variable, i.e. constant.
! In future it should be dependent on the level of osteoblast signalling, i.e. on work required.
! The inflow rate is specified as cells/hour, in_per_hour
!------------------------------------------------------------------------------------------------
subroutine influx
real :: rate, dr
real(8) :: R
integer :: i, n, kpar = 0

rate = in_per_hour*DELTA_T/60 
n = rate
dr = rate - n
do i = 1,n
	call mono_entry
enddo
if (dr > 0) then
!	call random_number(R)
	R = par_uni(kpar)
	if (R < dr) then
		call mono_entry
	endif
endif
end subroutine

!------------------------------------------------------------------------------------------------
! Choose a random location on the side of a capillary for monocyte initial location
!------------------------------------------------------------------------------------------------
subroutine old_mono_entry
integer :: icap, kpar=0
real(8) :: p(MAX_CAP), R
real :: d, v(3), alfa

do icap = 1,ncap
	p(icap) = capillary(icap)%surface_area
enddo
p(1:ncap) = p(1:ncap)/sum(p(1:ncap))

icap = random_selection(p,ncap)
R = par_uni(kpar)
d = 4 + R*(capillary(icap)%length - 8)		! leave out the ends
alfa = d/capillary(icap)%length
v = (1-alfa)*capillary(icap)%pos1 + alfa*capillary(icap)%pos2
end subroutine

!------------------------------------------------------------------------------------------------
! Choose a random location on the side of a capillary for monocyte initial location
!------------------------------------------------------------------------------------------------
subroutine mono_entry
integer :: i, site(3), kpar=0

do
	i = random_int(1,nentrysites,kpar)
	site = entrysite(:,i)
	if (occupancy(site(1),site(2),site(3))%indx == 0) exit
enddo
call addMono('blood',site)
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine update_S1P1(pmono)
type(monocyte_type), pointer :: pmono
real :: S

S = pmono%S1P1
S = S + rate_S1P1(S)*DELTA_T
pmono%S1P1 = min(S,1.0)
end subroutine

!------------------------------------------------------------------------------------------------
! This has been changed to allow RANK signalling only when the monocyte is in contact with an OB.
!------------------------------------------------------------------------------------------------
subroutine update_RANK(pmono)
type(monocyte_type), pointer :: pmono
real :: S, signal
integer :: site(3), status, iblast

S = pmono%RANKSIGNAL
site = pmono%site
status = pmono%status
!C = CXCL12_conc(site(1),site(2),site(3))
!S = (1-RANKSIGNAL_decayrate)*S + rate_RANKSIGNAL(S,C)*DELTA_T
iblast = NearOB(pmono)
if ((iblast > 0 .and. status <= STICKY) .or. status == CLUMPED) then
	signal = rate_RANKSIGNAL(S)
else
	signal = 0
endif
S = (1-RANKSIGNAL_decayrate)*S + signal*DELTA_T
pmono%RANKSIGNAL = min(S,1.0)
if (status == MOTILE .and. pmono%RANKSIGNAL >= ST1) then
	pmono%status = CHEMOTACTIC
!	call logger('CHEMOTACTIC')
endif
if (status == CHEMOTACTIC .and. pmono%RANKSIGNAL < ST1/2) then
	pmono%status = MOTILE		! DEAD
endif
if (status == STICKY .and. pmono%RANKSIGNAL < ST2/2) then
	pmono%status = CHEMOTACTIC	! DEAD
endif
if (status == CHEMOTACTIC .and. pmono%RANKSIGNAL > ST2 .and. NearOB(pmono) > 0) then
	pmono%status = STICKY
endif
if (pmono%status == STICKY) then
	pmono%stickiness = pmono%RANKSIGNAL
endif
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine updater
real :: S
integer :: i, j, k, iclump, iclast, iblast, irel, dir, region, kcell, site(3), res, kpar=0
integer :: ir, dx, dz, x, z
real :: tnow, stickysum, rmargin
real(8) :: R
type(monocyte_type), pointer :: pmono
type(osteoclast_type), pointer :: pclast
type(osteoblast_type), pointer :: pblast
type(clump_type), pointer :: pclump
logical :: on_surface

tnow = istep*DELTA_T
! Monocytes enter from the blood
call influx
! Monocyte S1P1 level grows, RANKL signal accumulates
do i = 1,nmono
	pmono => mono(i)
	if (pmono%status == LEFT) cycle
	if (pmono%region /= MARROW) cycle
	if (pmono%status < 1 .or. pmono%status > FUSING) cycle
	call update_S1P1(pmono)
	call update_RANK(pmono)
enddo
do iclump = 1,nclump
	pclump => clump(iclump)
	if (pclump%status < 0) cycle
!	if (pclump%status < FUSED) then
!		call consolidate_clump(pclump)
!	endif
!	if (nclump > 1) then
!		call separate_clump(pclump)
!	endif
	if (pclump%status < FUSING) then
		stickysum = 0
		do i = 1,pclump%ncells
			stickysum = stickysum + mono(pclump%list(i))%stickiness
		enddo
		if (stickysum/pclump%ncells < ST2) then	! clump breaks up
			write(logmsg,*) 'Clump breaks up: ',iclump
			call logger(logmsg)
			do i = 1,pclump%ncells
				mono(pclump%list(i))%status = CHEMOTACTIC
			enddo
			call RemoveClump(pclump)
			cycle
		endif
		call detacher(pclump)
		if (pclump%ncells > CLUMP_THRESHOLD) then
			pclump%status = FUSING
			pclump%starttime = tnow
			pclump%fusetime = tnow + FUSING_TIME
			do i = 1,pclump%ncells
				mono(pclump%list(i))%status = FUSING
			enddo
			write(logmsg,*) 'Started FUSING: ',iclump
			call logger(logmsg)
		endif
	elseif (pclump%status == FUSING) then
		if (tnow >= pclump%fusetime) then
!			write(logmsg,*) '***fuse clump: ',iclump,pclump%fusetime,tnow,pclump%cm 
!			call logger(logmsg)
			call fuse_clump(pclump)
			write(logmsg,*) 'FUSED: ',iclump
			call logger(logmsg)
		endif
	elseif (pclump%status == FUSED) then
		call lower_clump(pclump,on_surface)
		if (on_surface) then
			write(logmsg,*) 'On surface: ',iclump
			call logger(logmsg)
			call createOsteoclast(pclump)
		endif
	endif
enddo

! Stem cells divide
do i = 1,NSTEM
	if (tnow > stem(i)%dividetime) then
!		call random_number(R)
		R = par_uni(kpar)
		irel = nreldir*R
		do k = 1,nreldir
			irel = irel + 1
			if (irel > nreldir) irel = 1
			dir = reldir(1,irel)
			site = stem(i)%site + jumpvec(:,dir)
			if (free_site(site,region,kcell)) then
				call addMono('stem',site)
				stem(i)%dividetime = tnow + STEM_CYCLETIME
				nborn = nborn + 1
				exit
			endif
		enddo
	endif
enddo

! Osteoclasts complete fusing, dissolve bone, move, or die.
do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status == DEAD) cycle
	if (tnow > pclast%dietime) then
		call ClastDeath(iclast)
		cycle
	endif
!	if (pclast%status == FUSING) then
!		if (tnow >= pclast%entrytime) then
!			call completeFusing(iclast)
!			cycle
!		endif
!	endif

	! The seal on uneroded bone under and near the OC is gradually removed.
	rmargin = pclast%radius + OC_MARGIN
	ir = rmargin + 0.5
!	write(*,*) 'eraode seal: ',pclast%radius,OC_MARGIN,rmargin,ir
	do dx = -ir, ir
		do dz = -ir,ir
			if (dx**2 + dz**2 > rmargin**2) cycle
			x = pclast%site(1) + dx
			if (x < 1 .or. x > NX) cycle
			z = pclast%site(3) + dz
			if (z < 1 .or. z > NZ) cycle
!			if (surface(x,z)%target_depth == 0) cycle
			if (surface(x,z)%depth > 0) then
				surface(x,z)%seal = 0
			else
!				write(*,*) dx,dz,x,z
				surface(x,z)%seal = max(0.0,surface(x,z)%seal - SEAL_REMOVAL_RATE*DELTA_T)
			endif
		enddo
	enddo
	
	if (tnow > pclast%movetime) then
!		write(*,*) 'Move osteoclast: ',tnow,iclast
		call MoveClast(pclast,res)
		if (res == 0) then			! no move
			pclast%movetime = tnow + CLAST_DWELL_TIME
			pclast%blocktime = 0
			pclast%status = RESORBING
		elseif (res == 1) then		! moved
			pclast%movetime = tnow + CLAST_DWELL_TIME
			pclast%blocktime = 0
			pclast%status = RESORBING
		elseif (res == 2) then		! osteoclast blocked	
			call unblocker
			pclast%movetime = tnow + CLAST_DWELL_TIME
			pclast%blocktime = 0
			pclast%status = RESORBING

!			if (pclast%blocktime == 0) then
!				write(logmsg,*) 'OC blocked: ',iclast,tnow
!				call logger(logmsg)
!				pclast%blocktime = tnow
!			elseif (tnow - pclast%blocktime > CLAST_STOP_TIME) then
!				write(logmsg,*) 'OC blocked too long: ',iclast,tnow
!				call logger(logmsg)
!				call clastDeath(iclast)
!				cycle
!			endif
!			! osteoclast needs to be checked more often for possible move
!			pclast%movetime = tnow + CLAST_DWELL_TIME/10
!			pclast%status = ALIVE
			
		elseif (res == 3) then		
			call logger('No signal left, OC dies')
			call ClastDeath(iclast)
			cycle
		endif
	endif
	call resorber(pclast)
enddo
! Osteoblast state changes and motion
do iblast = 1,nblast
	pblast => blast(iblast)
	site = pblast%site
	if (pblast%status == DORMANT .and. surface(site(1),site(3))%seal < 1) then
		pblast%status = ALIVE
	endif
	if (pblast%status == DEAD) cycle
!	if (tnow > pblast%movetime) then
!		call MoveBlast(pblast,res)
!		pblast%movetime = tnow + BLAST_DWELL_TIME
!	endif
!	if (mod(istep,60*12) == 0) then
!		write(logmsg,'(a,i4,f8.4)') 'blastsignal: ',iblast,BlastSignal(pblast)
!		call logger(logmsg)
!	endif
enddo
end subroutine

!------------------------------------------------------------------------------------------------
! Now use surface(:,:)%depth and %signal instead of pclast%pit and bone_fraction.
! If the OC has never received any bone signal, no resorption.
!------------------------------------------------------------------------------------------------
subroutine resorber(pclast)
type(osteoclast_type), pointer :: pclast
!type(occupancy_type), pointer :: pbone
integer :: k, x, z, site(3), kk
real :: totsig

if (pclast%status == MOVING) then
	return
elseif (pclast%status /= RESORBING) then
	site = pclast%site
	totsig = TotalSignal(pclast,site)
	if (totsig == 0) then
		return
	endif
	pclast%status = RESORBING
endif
call pitrates(pclast)
do k = 1,pclast%npit
!	site = pclast%pit(k)%site
!	pbone => occupancy(site(1),site(2),site(3))
!	if (pbone%bone_fraction > 0) then
!		pbone%bone_fraction = pbone%bone_fraction - pclast%pit(k)%rate*DELTA_T
!		if (pbone%bone_fraction <= 0) then
!			pbone%bone_fraction = 0
!			pbone%region = PIT
!			pbone%indx = 0
!			if (site(2) > 1) then
!				pclast%pit(k)%site(2) = site(2) - 1
!			endif
!		endif
!	endif
	x = pclast%site(1) + pclast%pit(k)%delta(1)
	z = pclast%site(3) + pclast%pit(k)%delta(3)
	surface(x,z)%depth =  surface(x,z)%depth + pclast%pit(k)%rate*DELTA_T
	surface(x,z)%signal = max(0.0,GetSignal(x,z))
	if (surface(x,z)%depth > 1.1*surface(x,z)%target_depth) then
		write(logmsg,*) 'target_depth exceeded: ',pclast%ID,k,x,z,surface(x,z)%depth,surface(x,z)%target_depth,surface(x,z)%signal
		call logger(logmsg)
!		do kk = 1,pclast%npit
!			x = pclast%cm(1) + pclast%pit(kk)%delta(1) + 0.5
!			z = pclast%cm(3) + pclast%pit(kk)%delta(3) + 0.5
!			surface(x,z)%depth =  surface(x,z)%depth + pclast%pit(kk)%rate*DELTA_T
!			surface(x,z)%signal = max(0.0,(surface(x,z)%target_depth - surface(x,z)%depth)/MAX_PIT_DEPTH)
!			write(logmsg,'(3i4,3f8.3,f8.4)') kk,x,z,surface(x,z)%depth,surface(x,z)%target_depth,surface(x,z)%signal,pclast%pit(kk)%rate
!			call logger(logmsg)
!		enddo
!		stop
	endif
	if (surface(x,z)%depth > 1.1*MAX_PIT_DEPTH) then
		write(logmsg,*) 'MAX_PIT_DEPTH exceeded: ',x,z,surface(x,z)%depth,surface(x,z)%target_depth,surface(x,z)%signal
		call logger(logmsg)
		stop
	endif
enddo
end subroutine

!------------------------------------------------------------------------------------------------
! After fusing is complete the clump settles down onto the bone surface and becomes an osteoclast.
! Currently no more monocytes can join a clump after it has fused.
!------------------------------------------------------------------------------------------------
subroutine fuse_clump(pclump)
type(clump_type) :: pclump
integer :: i, iblast

!call logger('fusing')
pclump%status = FUSED
do i = 1,pclump%ncells
	mono(pclump%list(i))%status = FUSED
enddo
!icm = pclump%cm + 0.5
!pclump%cm = icm
!iblast = pclump%iblast
!blast(iblast)%status = DORMANT
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine addMono(source,site)
character*(*) :: source
integer :: site(3)
integer :: kpar = 0

if (nmono == MAX_MONO) then
	call logger('Error: addMono: Exceeded monocyte array, increase MAX_MONO')
	stop
endif
nmono = nmono + 1
mono_cnt = mono_cnt + 1
mono(nmono)%ID = nmono
mono(nmono)%site = site
mono(nmono)%region = MARROW
mono(nmono)%status = MOTILE
mono(nmono)%lastdir = random_int(1,6,kpar)
mono(nmono)%iclump = 0
mono(nmono)%S1P1 = 0
mono(nmono)%RANKSIGNAL = 0
mono(nmono)%stickiness = 0
mono(nmono)%lastmovestep = istep
!nullify(mono(nmono)%clump)
occupancy(site(1),site(2),site(3))%species = MONOCYTE
!if (.not.use_TCP) then
!	write(*,*) 'added monocyte: ',nmono,site,'  ',source
!else
!	call logger('added monocyte  '//source)
!endif
end subroutine

!------------------------------------------------------------------------------------------------
! S1P1 (S1P receptor) expression is currently assumed to grow at a steady rate.
!------------------------------------------------------------------------------------------------
real function rate_S1P1(S)
real :: S
rate_S1P1 = S1P1_BASERATE
end function

!------------------------------------------------------------------------------------------------
! The rate of RANK signalling depends on the RANKL concentration (C) and the current integrated
! signal level (S).
!------------------------------------------------------------------------------------------------
real function rate_RANKSIGNAL1(S,C)
real :: S, C
rate_RANKSIGNAL1 = RANKSIGNAL_rateconstant*(1-S)*C
end function

!------------------------------------------------------------------------------------------------
! The rate of RANK signalling depends on the current integrated signal level (S),
! and possibly on the OB signal strength.  For now just on S - saturating when S = 1.
!------------------------------------------------------------------------------------------------
real function rate_RANKSIGNAL(S)
real :: S
rate_RANKSIGNAL = RANKSIGNAL_rateconstant*(1-S)
end function

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine completeFusing(iclast)
integer :: iclast
integer :: i
real :: tnow

tnow = istep*DELTA_T
write(logmsg,*) 'completed fusing: ',iclast,tnow
call logger(logmsg)
clast(iclast)%status = ALIVE
clast(iclast)%movetime = tnow + CLAST_DWELL_TIME 
clast(iclast)%dietime = tnow + clastLifetime()
!do i = 1,clast(iclast)%count
!	mono(clast(iclast)%mono(i))%status = FUSED
!enddo
end subroutine

!------------------------------------------------------------------------------------------------
! Change pit sites to relative displacements from osteoclast site
! Need to redo the pits.  It makes more sense to:
! (a) Make OCs circular
! (b) Precompute the dissolving rate factor as a function of distance from the OC centre,
!     and accounting for partial coverage of a grid site.  A site is included if the 
!     centre is within the radius.
! (c) Remove seal from under OCs initial location
!------------------------------------------------------------------------------------------------
subroutine createOsteoclast(pclump)
type(clump_type), pointer :: pclump
integer :: bonesite(3,100), i, j, n, npit, dx, dz, x, z
type(osteoclast_type), pointer :: pclast
real :: tnow, r, fract(100), fsum
integer :: kpar = 0

tnow = istep*DELTA_T
nclast = nclast + 1
nliveclast = nliveclast + 1
write(logmsg,*) 'createOsteoclast: ',nclast,tnow
call logger(logmsg)
pclast => clast(nclast)
pclast%ID = nclast
!pclast%site = pclump%cm + 0.5
pclast%site = pclump%site
pclast%site(2) = NBY + 1
pclast%lastdir = random_int(1,8,kpar)
pclast%entrytime = tnow
pclast%status = ALIVE
pclast%movetime = tnow 
pclast%dietime = tnow + clastLifetime()
pclast%count = pclump%ncells
do i = 1,pclump%ncells
	j = pclump%list(i)
	mono(j)%iclast = nclast
	mono(j)%status = OSTEO
enddo
! Now need to create the pit list
pclast%radius = clast_size(pclast%count)
n = pclast%radius + 1
npit = 0
fsum = 0
do dx = -n,n
	do dz = -n,n
		r = sqrt(real(dx*dx + dz*dz))
		if (r < pclast%radius) then
			npit = npit + 1
			bonesite(:,npit) = pclast%site + (/dx,0,dz/)
			fract(npit) = 1 - (1-0.5)*r/pclast%radius
			if (r > pclast%radius - 0.5) then
				fract(npit) = fract(npit)*(pclast%radius + 0.5 - r)
			endif
			fsum = fsum + fract(npit)
		endif
	enddo
enddo
pclast%npit = npit
allocate(pclast%pit(npit))
! Make pit location an offset from the centre, set dissolving rate fraction
do i = 1,npit
	pclast%pit(i)%delta = bonesite(:,i) - pclast%site
	pclast%pit(i)%fraction = fract(i)*(npit/fsum)		! normalize %fraction to sum to npit
	x = bonesite(1,i)
	z = bonesite(3,i)
	surface(x,z)%seal = 0
enddo
call RemoveClump(pclump)
call pitrates(pclast)
write(logmsg,'(a,7i4)') 'clast site, count, npit: ',pclast%site,pclast%count,pclast%npit
call logger(logmsg)
call UpdateSurface
end subroutine


!------------------------------------------------------------------------------------------------
! The bone resorption rate at a given pit site (x,z) depends on:
!	Nm = the number of monocytes that fused to make the osteoclast
!   Np = number of pit sites that the osteoclast covers
!	(d = the distance of the target bone site from the osteoclast centre
!   The depth factor df decreases linearly to zero as d goes from 0 to MAX_RESORPTION_D)
! The volume rate of resorption per pit site is:
!   (MAX_RESORPTION_RATE/MAX_RESORPTION)*(Nm/Np) um^3/day
! which is converted to grids/min, /(24*60*DELTA_X^3)
! Note: currently not using d dependence, using %fraction
!------------------------------------------------------------------------------------------------
real function resorptionRate(Nm,Np)
integer :: Nm, Np
!real :: d
!real :: df

!if (d >= MAX_RESORPTION_D) then
!	df = 0
!else
!	df = 1 - d/MAX_RESORPTION_D
!endif
resorptionRate = (MAX_RESORPTION_RATE*Nm)/(MAX_RESORPTION_N*Np*24*60*DELTA_X**3)
!write(*,*) 'resorptionRate: ',Nm,Np,resorptionRate
!write(*,*) MAX_RESORPTION_RATE,Nm,MAX_RESORPTION_N,Np,24*60,DELTA_X 
end function

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine ClastDeath(i)
integer :: i
integer :: k, site(3), kcell
integer :: x,y,z

!do k = 1,clast(i)%npit
!	if (clast(i)%pit(k)%fraction > 0 .and. clast(i)%pit(k)%fraction < 0.5) then
!		site = clast(i)%pit(k)%site
!		occupancy(site(1),site(2),site(3))%region = PIT
!		occupancy(site(1),site(2),site(3))%indx = 0
!	endif
!enddo
write(logmsg,*) 'clast death: ',i
call logger(logmsg)
clast(i)%status = DEAD
deallocate(clast(i)%pit)
!do k = 1,clast(i)%count
!	kcell = clast(i)%mono(k)
!	site = mono(kcell)%site
!	mono(kcell)%status = DEAD
!	occupancy(site(1),site(2),site(3))%indx = 0
!enddo
nliveclast = nliveclast - 1
call UpdateSurface
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
real function clastLifetime()

clastLifetime = CLAST_LIFETIME
end function

!------------------------------------------------------------------------------------------------
! Determine approximate normal to the bone surface at site().
! Method: inspect all neighbours of site, average the (x,y,z) of locations in marrow.
! The normal vector is from site to this average point.
!------------------------------------------------------------------------------------------------
subroutine boneNormal(site,v)
integer :: site(3)
real :: v(3)
integer :: k, irel, cnt, site1(3)

v = 0
cnt = 0
do irel = 1,nreldir
	k = reldir(1,irel)
	site1 = site + jumpvec(:,k)
	if (occupancy(site1(1),site1(2),site1(3))%region == MARROW) then
		cnt = cnt + 1
		v = v + site1
	endif
enddo
v = v/cnt - site
v = v/rnorm(v)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(NX_dim,NY_dim,NZ_dim,NBY_dim) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: NX_dim,NY_dim,NZ_dim,NBY_dim

NX_dim = NX
NY_dim = NY
NZ_dim = NZ
NBY_dim = NBY
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_summary(summaryData) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*)

summaryData(1:4) = (/istep,mono_cnt,nborn,nleft/)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine CheckBlast
integer :: iblast, site(3)
type(osteoblast_type), pointer :: pblast

do iblast = 4,7,3
	pblast => blast(iblast)
	site = pblast%site
	write(logmsg,'(a,3i4,2f8.4)') 'OB: ',iblast,pblast%status,surface(site(1),site(3))%iclast,BlastSignal(pblast),surface(site(1),site(3))%signal
	call logger(logmsg)
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine CheckMono
integer :: imono, stat, cnt(9)
integer :: clumpcnt(100), iclump, icmax

cnt = 0
clumpcnt = 0
icmax = 0
do imono = 1,nmono
	stat = mono(imono)%status
	if (stat < 1) cycle
	cnt(stat) = cnt(stat) + 1
	if (stat == CLUMPED) then
		iclump = mono(imono)%iclump
		clumpcnt(iclump) = clumpcnt(iclump) + 1
		icmax = max(icmax,iclump)
	endif
	if (stat == FUSING) then
		iclump = mono(imono)%iclump
		clumpcnt(iclump) = clumpcnt(iclump) + 1
		icmax = max(icmax,iclump)
	endif
	if (stat == FUSED) then
		iclump = mono(imono)%iclump
		clumpcnt(iclump) = clumpcnt(iclump) + 1
		icmax = max(icmax,iclump)
	endif
enddo
write(logmsg,'(a,i8,6i6,2f8.4,i6)') 'mono: ',istep,icmax,cnt(3:7)
call logger(logmsg)
do iclump = 1,icmax
	if (clumpcnt(iclump) /= clump(iclump)%ncells) then
		write(logmsg,*) 'ERROR: CheckMono: bad clumpcnt: ',iclump,clumpcnt(iclump),clump(iclump)%ncells
		call logger(logmsg)
		stop
	endif
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_scene(ncap_list,cap_list,nmono_list,mono_list,npit_list,pit_list, &
					nclast_list,clast_list,nblast_list,blast_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene
use, intrinsic :: iso_c_binding
integer(c_int) :: nmono_list, mono_list(*), ncap_list, npit_list, nclast_list, nblast_list, blast_list(*)
real(c_float) :: cap_list(*), pit_list(*), clast_list(*)
integer :: k, j, kcell, iclump, site(3), nfused, x, y, z, status, mono_state, icap, iclast, iblast
real :: tnow, t1, t2, fraction, ypit, size, lastjump(3)
type(osteoclast_type), pointer :: pclast
type(osteoblast_type), pointer :: pblast
real :: clast_diam = 0.9
real :: mono_diam = 0.5

ncap_list = 0
nmono_list = 0
npit_list = 0
nclast_list = 0 
nblast_list = 0

tnow = istep*DELTA_T
! Monocyte section
if (nmono > 0) then
	nfused = 0
	k = 0
    do kcell = 1,nmono
		status = mono(kcell)%status
		if (status == DEAD .or. status == LEFT) cycle
		if (mono(kcell)%region /= MARROW) cycle
		if (status == CROSSING) then
			mono_state = 1
		elseif (status == CLUMPED) then
			mono_state = 2
		elseif (status == FUSING) then
			iclump = mono(kcell)%iclump
			t2 = clump(iclump)%fusetime
			t1 = clump(iclump)%starttime
			fraction = min(1.0,(tnow-t1)/(t2-t1))
			mono_state = 2 + fraction*98
!			if (fraction < 0) then
!				write(logmsg,*) 'fraction = 0: ',istep,kcell
!				call logger(logmsg)
!			endif
		elseif (status == FUSED) then
			mono_state = 100
			nfused = nfused + 1
		elseif (status == OSTEO) then
			cycle
		elseif (status == STICKY) then	
			mono_state = 2
		else
            mono_state = 0
		endif
        site = mono(kcell)%site
        if (FAST_DISPLAY .and. mono_state < 1) cycle
		k = k+1
		j = 5*(k-1)
		mono_list(j+1) = kcell-1
		mono_list(j+2:j+4) = site
		mono_list(j+5) = mono_state
!		write(logmsg,*) k,mono_state
!		call logger(logmsg)
    enddo
endif
nmono_list = k

! Capillary section
!do icap = 1,ncap
!	write(nfpos,'(a,6f6.1,f5.2)') 'C ',capillary(icap)%pos1,capillary(icap)%pos2,capillary(icap)%radius + 0.25
!enddo
do icap = 1,ncap
	j = (icap-1)*7
	cap_list(j+1:j+3) = capillary(icap)%pos1
	cap_list(j+4:j+6) = capillary(icap)%pos2
	cap_list(j+7) = capillary(icap)%radius+0.25
enddo
ncap_list = ncap

! Pit section
y = NBY
k = 0
do x = 1,NX
	do z = 1,NZ
!		do y = 1,NBY
!			if (occupancy(x,y,z)%region == PIT .or. &
!			   (occupancy(x,y,z)%region == BONE .and. occupancy(x,y,z)%bone_fraction < 1.0)) then
!			   ypit = y - 0.5 + occupancy(x,y,z)%bone_fraction
!			   k = k+1
!			   j = 4*(k-1)
!			   pit_list(j+1:j+3) = (/x,y,z/)
!			   pit_list(j+4) = ypit
!				write(nfpos,'(a,3i4,f7.3)') 'P ',x,y,z,ypit
!				exit
!			endif
!		enddo

!		if (surface(x,z)%depth > 0) then
!			ypit = y + 0.5 - surface(x,z)%depth
		if (surface(x,z)%target_depth > 0) then
!			if (surface(x,z)%signal == 0) then
!				ypit = 99
!			else
!				ypit = y + 0.5 - surface(x,z)%target_depth
!			endif
			ypit = surface(x,z)%signal
			k = k+1
			j = 4*(k-1)
			pit_list(j+1:j+3) = (/x,y,z/)
			pit_list(j+4) = ypit
!			write(nflog,'(a,3i4,f7.3)') 'P ',x,y,z,ypit
		endif
	enddo
enddo
npit_list = k

! Osteoclast section
k = 0
do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status == DEAD) cycle
	k = k+1
	size = pclast%radius
	j = 7*(k-1)
	clast_list(j+1:j+3) = pclast%site
	clast_list(j+2) = clast_list(j+2) - 0.5
!	clast_list(j+4:j+6) = (/1,0,0/)		! direction unit vector
	lastjump = dir2D(:,pclast%lastdir)
	clast_list(j+4:j+6) = lastjump/sqrt(dot_product(lastjump,lastjump))
	clast_list(j+7) = size
enddo
nclast_list = k

! Osteoblast section
k = 0
do iblast = 1,nblast
	pblast => blast(iblast)
	if (pblast%status == DEAD) cycle
	k = k+1
	j = 5*(k-1)
	blast_list(j+1) = pblast%ID-1
	blast_list(j+2:j+4) = pblast%site
	blast_list(j+5) = pblast%status
!	write(nflog,*) 'OB: ',blast_list(j+1:j+5)
enddo
nblast_list = k
!write(logmsg,'(a,5i6)') '# of mono, cap, pit, OC, OB: ',nmono_list, ncap_list, npit_list, nclast_list, nblast_list
!call logger(logmsg)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_pause
logical :: paused

inquire(file=pausefile,exist=paused)
if (paused) then
	call logger('Pause order received')
	do while (paused)
        call sleeper(1)   ! Too coarse!
		inquire(file=pausefile,exist=paused)
	enddo
	call logger('Resuming ...')
endif
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine read_inputfile(ok)
logical :: ok
real :: TC_AVIDITY_MEAN,TC_AVIDITY_SHAPE,TC_STIM_RATE_CONSTANT,TC_STIM_HALFLIFE,divide_mean1,divide_shape1
real :: DC_LIFETIME_MEAN,DC_LIFETIME_SHAPE
real :: IL2_THRESHOLD,ACTIVATION_THRESHOLD,FIRST_DIVISION_THRESHOLD,DIVISION_THRESHOLD,EXIT_THRESHOLD,STIMULATION_LIMIT
real :: chemo_radius,chemo_K_exit,chemo_K_DC
real :: S1P1_RISETIME

ok = .false.
open(nfinp,file=inputfile,status='old')

!read(nfinp,*) TC_AVIDITY_MEAN				! mean of avidity distribution (only if fix_avidity = false)
!read(nfinp,*) TC_AVIDITY_SHAPE			    ! shape -> 1 gives normal dist with small variance
!read(nfinp,*) TC_STIM_RATE_CONSTANT			! rate const for TCR stimulation (-> molecules/min)
!read(nfinp,*) TC_STIM_HALFLIFE				! halflife of T cell stimulation (hours)
!read(nfinp,*) divide_mean1
!read(nfinp,*) divide_shape1
read(nfinp,*) MONOCYTE_DIAMETER				! monocyte diameter (um) = 10	! um
read(nfinp,*) BETA							! speed: 0 < beta < 1		(0.65)
read(nfinp,*) RHO							! persistence: 0 < rho < 1	(0.95)
read(nfinp,*) S1P_CHEMOLEVEL
read(nfinp,*) S1P_KDIFFUSION
read(nfinp,*) S1P_KDECAY
read(nfinp,*) S1P_GRADLIM
read(nfinp,*) S1P1_THRESHOLD
read(nfinp,*) S1P1_RISETIME

!read(nfinp,*) DC_LIFETIME_MEAN				! days
!read(nfinp,*) DC_LIFETIME_SHAPE 			! days

!read(nfinp,*) IL2_THRESHOLD					! stimulation needed to initiate IL-2/CD25 production
!read(nfinp,*) ACTIVATION_THRESHOLD			! stimulation needed for activation
!read(nfinp,*) FIRST_DIVISION_THRESHOLD		! activation level needed for first division
!read(nfinp,*) DIVISION_THRESHOLD			! activation level needed for subsequent division
!read(nfinp,*) EXIT_THRESHOLD				! activation level below which exit is permitted
!read(nfinp,*) STIMULATION_LIMIT				! maximum activation level

read(nfinp,*) X_SIZE						! size of bone region (square)
read(nfinp,*) Y_SIZE						! thickness of slice

read(nfinp,*) CAPILLARY_DIAMETER			! capillary diameter (um) (= 3)
read(nfinp,*) MONO_PER_MM3					! initial (equil) number of monocytes/mm3 (= 2000)
read(nfinp,*) IN_PER_HOUR					! rate of influx of monocytes from the blood
read(nfinp,*) STEM_PER_MM3					! number of stem cells/mm3
read(nfinp,*) STEM_CYCLETIME				! stem cell division cycle time (hours) (= 6*60	! 6 hours)
read(nfinp,*) CROSSING_TIME					! time taken for a monocyte to cross into the capillary (mins)

read(nfinp,*) FUSING_TIME					! time taken by monocytes fusing into an osteoclast	(120) (mins)
read(nfinp,*) CLAST_LIFETIME				! lifetime of an osteoclast (4) days -> mins
read(nfinp,*) CLAST_DWELL_TIME				! time an osteoclast spends in one spot (180) (mins)
read(nfinp,*) MAX_RESORPTION_RATE			! maximum bone removal rate (/grid cell) (0.02) (um/min)
read(nfinp,*) MAX_RESORPTION_D				! maximum pit depth (for scaling resorption rate) (50) (um)
read(nfinp,*) MAX_RESORPTION_N				! number of monos in osteoclast corresponding to MAX_RESORPTION_RATE (30)

!read(nfinp,*) SIGNAL_RADIUS					! radius of influence of bone signal (um -> grids) (10)
!read(nfinp,*) SIGNAL_THRESHOLD				! defines the high-signal region, near the source (0.14)
!read(nfinp,*) SIGNAL_AFACTOR				! field amplification factor (0.4)
!read(nfinp,*) MTHRESHOLD					! number of monocytes in the high-signal region that triggers fusing (25)

read(nfinp,*) cross_prob					! probability (/timestep) of monocyte egress to capillary
!read(nfinp,*) chemo_radius					! radius of chemotactic influence (sites)
!read(nfinp,*) chemo_K_exit					! level of chemotactic influence towards exits
!read(nfinp,*) chemo_K_DC					! level of chemotactic influence towards DCs

read(nfinp,*) days							! number of days to simulate
read(nfinp,*) seed(1)						! seed vector(1) for the RNGs
read(nfinp,*) seed(2)						! seed vector(2) for the RNGs
read(nfinp,*) ncpu							! number of threads, not used currently
read(nfinp,*) NT_GUI_OUT					! interval between GUI outputs (timesteps) 
read(nfinp,*) SPECIES						! animal species
close(nfinp)
!call logger('Finished reading input file')

DELTA_X = MONOCYTE_DIAMETER
NX = X_SIZE/DELTA_X									! convert um -> grids
NY = Y_SIZE/DELTA_X	+NBY							! convert um -> grids
if (mod(NX,2) /= 0) NX = NX+1						! ensure that NX is even (why?)
NZ = NX
CAPILLARY_DIAMETER = CAPILLARY_DIAMETER/DELTA_X		! convert um -> grids
!chemo_radius = chemo_radius/DELTA_X					! convert um -> grids
STEM_CYCLETIME = 60*STEM_CYCLETIME					! convert hours -> minutes
CLAST_LIFETIME = CLAST_LIFETIME*24*60				! convert days -> minutes
!MAX_RESORPTION_RATE = MAX_RESORPTION_RATE/DELTA_X	! convert um/min -> grids/min
MAX_RESORPTION_D = MAX_RESORPTION_D/DELTA_X			! convert um -> grids
SIGNAL_RADIUS = SIGNAL_RADIUS/DELTA_X				! convert um -> grids
S1P1_BASERATE = 1./(60.*S1P1_RISETIME)				! convert time (hours) to rate (/min)

!write(*,*) 'DC_RADIUS, chemo_radius: ',DC_RADIUS,chemo_radius
!chemo_N = max(3,int(chemo_radius + 0.5))	! convert from um to lattice grids
!write(*,*) 'chemo_N: ',chemo_N
! Note that currently exit chemotaxis and DC chemotaxis are treated in the same way - same decay
!chemo_exp = log(1/CHEMO_MIN)/log(chemo_radius)

!sigma = log(TC_AVIDITY_SHAPE)
!TC_AVIDITY_MEDIAN = TC_AVIDITY_MEAN/exp(sigma*sigma/2)
!sigma = log(DC_ANTIGEN_SHAPE)
!DC_ANTIGEN_MEDIAN = DC_ANTIGEN_MEAN/exp(sigma*sigma/2)
!sigma = log(DC_LIFETIME_SHAPE)
!DC_LIFETIME_MEDIAN = DC_LIFETIME_MEAN/exp(sigma*sigma/2)

!sigma = log(divide_shape1)
!divide_dist1%p1 = log(60*divide_mean1/exp(sigma*sigma/2))
!divide_dist1%p2 = sigma
!sigma = log(divide_shape2)
!divide_dist2%p1 = log(60*divide_mean2/exp(sigma*sigma/2))
!divide_dist2%p2 = sigma
!write(*,*) 'divide_dist2: ',divide_dist2

!if (chemo_K_exit == 0.0) then
!    ep_factor = 2.4
!elseif (chemo_K_exit <= 0.2) then
!    ep_factor = 2.3
!elseif (chemo_K_exit <= 0.4) then
!    ep_factor = 2.2
!elseif (chemo_K_exit <= 0.6) then
!    ep_factor = 2.15
!elseif (chemo_K_exit <= 0.8) then
!    ep_factor = 2.1
!else
!    ep_factor = 2.1
!endif

!call setup_dists

Nsteps = days*60*24/DELTA_T

!open(nfout,file=outputfile,status='replace')
!if (save_input) then
!    call save_inputfile(inputfile)
!    call save_parameters
!	call save_inputfile(fixedfile)
!endif
!write(*,*) 'open resultfile: ',resultfile

ok = .true.

write(logmsg,*) 'FUSING_TIME: ',FUSING_TIME					! time taken by monocytes fusing into an osteoclast	(120) (mins)
call logger(logmsg)
write(logmsg,*) 'CLAST_LIFETIME: ',CLAST_LIFETIME				! lifetime of an osteoclast (4) days -> mins
call logger(logmsg)
write(logmsg,*) 'CLAST_DWELL_TIME: ',CLAST_DWELL_TIME				! time an osteoclast spends in one spot (180) (mins)
call logger(logmsg)
write(logmsg,*) 'MAX_RESORPTION_RATE: ',MAX_RESORPTION_RATE			! maximum bone removal rate (/grid cell) (0.02) (um/min)
call logger(logmsg)
write(logmsg,*) 'MAX_RESORPTION_D: ',MAX_RESORPTION_D				! maximum pit depth (for scaling rate) (5) (um)
call logger(logmsg)
write(logmsg,*) 'MAX_RESORPTION_N: ',MAX_RESORPTION_N				! number of monos in osteoclast corresponding to MAX_RESORPTION_RATE (30)
call logger(logmsg)

end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine connection(awp,port,error)
TYPE(winsockport) :: awp
integer :: port, error
integer :: address = 0
!!!character*(64) :: ip_address = "127.0.0.1"C      ! need a portable way to make a null-terminated C string
character*(64) :: host_name = "localhost"

if (.not.winsock_init(1)) then
    call logger("winsock_init failed")
    stop
endif
!write(nftemp,*) 'did winsock_init'

awp%handle = 0
awp%host_name = host_name
awp%ip_port = port
awp%protocol = IPPROTO_TCP
call Set_Winsock_Port (awp,error)

if (.not.awp%is_open) then
    write(nflog,*) 'Error: connection: awp not open: ',port
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
call connection(awp_0,TCP_PORT_0,error)
if (awp_0%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_0 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_0%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_0'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_0  '
call logger(logmsg)

if (use_CPORT1) then
	call connection(awp_1,TCP_PORT_1,error)
	if (awp_1%handle < 0 .or. error /= 0) then
		write(logmsg,'(a)') 'TCP connection to TCP_PORT_1 failed'
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (.not.awp_1%is_open) then
		write(logmsg,'(a)') 'No connection to TCP_PORT_1'
		call logger(logmsg)
		ok = .false.
		return
	endif
	write(logmsg,'(a)') 'Connected to TCP_PORT_1  '
	call logger(logmsg)
endif
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine snapshot
character*(128) :: msg
integer :: error

if (use_TCP) then
    if (.not.awp_1%is_open) then
        call logger("snapshot: awp_1 is not open")
        return
    endif
	write(msg,'(4i6)') istep,mono_cnt,nborn,nleft
!	call logger(msg)
    call winsock_send(awp_1,msg,len_trim(msg),error)
endif
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step
use, intrinsic :: iso_c_binding
integer(c_int) :: res
logical :: ok
integer :: error
real :: totsig, tnow, hour
integer, parameter :: NT_EVOLVE = 20

res = 0
!call logger("simulate_step")
ok = .true.
istep = istep + 1
tnow = istep*DELTA_T
if (.not. initiated .and. tnow > STARTUP_TIME*24*60) then
	call Initiation
endif
if (mod(istep,1000) == 0) then
	hour = istep/240.
	write(logmsg,'(a,i8,f8.2,6i6)') 'istep: ',istep,hour,mono_cnt,nleft	!,mono(22)%status,mono(22)%site
	call logger(logmsg)
!	call CheckMono
!	call CheckBlast
endif
!write(nflog,*) 'call updater'
call updater
!write(nflog,*) 'call MonoMover: ',nmono
call MonoMover
if (initiated .and. mod(istep,NT_EVOLVE) == 0) then
	call evolveCXCL12(NT_EVOLVE,totsig)
!	if (totsig == 0 .and. nclump > 0) then
!		call break_clumps
!	endif
!	nliveOC = 0
!	do ic = 1,nclast
!		if (clast(ic)%status /= DEAD) then
!			nliveOC = nliveOC + 1
!		endif
!	enddo
!	if (nliveOC == 0) then
!		call logger("All osteoclasts are dead")
!	endif
endif
if (nclast > 0 .and. nliveclast == 0) then
	res = -1
endif
end subroutine


!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine report
integer :: i, k, indx, status

write(logmsg,*) 'nclump: ',nclump
call logger(logmsg)
do i = 1,nclump
	if (clump(i)%status < 0) cycle
	write(logmsg,*) 'clump: ',i,'  status: ',clump(i)%status,'  ncells: ',clump(i)%ncells
	call logger(logmsg)
	do k = 1,clump(i)%ncells
		indx = clump(i)%list(k)
		write(logmsg,*) indx,mono(indx)%site
		call logger(logmsg)
	enddo
enddo
!write(logmsg,*) 'Clumped monocytes  '
!call logger(logmsg)
!do i = 1,nmono
!	status = mono(i)%status
!	if (status == LEFT .or. status == DEAD) cycle
!	if (mono(i)%iclump > 0) then
!		write(logmsg,*) i,mono(i)%iclump,mono(i)%site
!		call logger(logmsg)
!	endif
!enddo
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine wrapup
call par_zigfree
if (allocated(occupancy)) deallocate(occupancy)
if (allocated(surface)) deallocate(surface)
if (allocated(mono)) deallocate(mono)
if (allocated(clast)) deallocate(clast)
if (allocated(blast)) deallocate(blast)
if (allocated(stem)) deallocate(stem)
if (allocated(capillary)) deallocate(capillary)
if (allocated(entrysite)) deallocate(entrysite)
if (allocated(CXCL12_conc)) deallocate(CXCL12_conc)
if (allocated(CXCL12_grad)) deallocate(CXCL12_grad)
if (allocated(CXCL12_influx)) deallocate(CXCL12_influx)
if (allocated(S1P_conc)) deallocate(S1P_conc)
if (allocated(S1P_grad)) deallocate(S1P_grad)

end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run
use, intrinsic :: iso_c_binding
integer(c_int) :: res
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

!call report
call wrapup

if (res == 0) then
	call logger(' Execution successful')
else
	call logger('  === Execution failed ===')
	call sleeper(1)
endif
close(nflog)

if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    if (use_CPORT1) call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
!	    call logger("closed PORT_0")
		if (use_CPORT1) then
		    call winsock_send(awp_1,quit,8,error)
		    call winsock_close(awp_1)
!			call logger("closed PORT_1")
		endif
	endif
endif

end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine execute(infile_array,buflen) BIND(C)
!!DEC$ ATTRIBUTES DLLEXPORT :: EXECUTE
!!DEC$ ATTRIBUTES C, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"execute" :: execute
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128)
integer(c_int) :: buflen
character*(128) :: infile
logical :: ok
integer :: i, res

use_CPORT1 = .false.	! TESTING DIRECT CALLING FROM C++
infile = ''
do i = 1,buflen
	infile(i:i) = infile_array(i)
enddo

open(nflog,file='bone.log',status='replace')
awp_0%is_open = .false.
awp_1%is_open = .false.

#ifdef GFORTRAN
    write(logmsg,'(a)') 'Built with GFORTRAN'
	call logger(logmsg)
#endif

logmsg = 'OS??'
#ifdef LINUX
    write(logmsg,'(a)') 'OS is Linux'
#endif
#ifdef OSX
    write(logmsg,'(a)') 'OS is OS-X'
#endif
#ifdef _WIN32
    write(logmsg,'(a)') 'OS is Windows'
#endif
#ifdef WINDOWS
    write(logmsg,'(a)') 'OS is Windows'
#endif
call logger(logmsg)

!#ifdef OPENMP
#if defined(OPENMP) || defined(_OPENMP)
    write(logmsg,'(a)') 'Executing with OpenMP'
	call logger(logmsg)
#endif

write(logmsg,*) 'inputfile:  ', infile
call logger(logmsg)
if (use_tcp) then
	call logger('call connecter')
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
	call logger('did connecter')
endif
call setup(infile,ok)
call logger('did setup')
if (ok) then
	call logger('returned OK')
	clear_to_send = .true.
	simulation_start = .true.
	istep = 0
	
	return
	
	call simulate(ok)
	call logger('Ended simulation')
	call logger('Execution successful')
else
	call logger('setup failed')
endif
if (ok) then
	res = 0
else
	res = 1
endif
call terminate_run(res)
end subroutine

end module

