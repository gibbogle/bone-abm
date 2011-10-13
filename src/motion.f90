module motion
use global
use fields
use clumping
use rkf45

implicit none
save

real, parameter :: Kchemo = 1.0
real, parameter :: K_MUTUAL_ATTRACT = 0.6	!0.2 
real, parameter :: K_MUTUAL_REPELL = 0.3
real, parameter :: K_SIGNAL = 0.001	!0.001
real, parameter :: K_JOINING = 100
real, parameter :: SEAL_REMOVAL_RATE = 0.0002	!0.0001
real, parameter :: JOIN_DISP_THRESHOLD = 0.0005
real, parameter :: JOIN_TIME_LIMIT = 100
real, parameter :: K_OC_DRAG = 10.0
real, parameter :: OC_MASS = 0.1
real, parameter :: SPEEDUP = 1
logical, parameter :: SEAL_SWITCH = .true.

! For testing OC dynamics, need lookup tables to translate between
! OC index and state index (actually (index-1)/4 + 1)
type(osteoclast_type), allocatable, target :: OClist(:)
real(8), allocatable :: OCstate(:), OCstatep(:)
integer :: nOClist
integer :: OC_NV
integer :: OC_to_index(100)
integer :: index_to_OC(100)
integer, parameter :: NF=100
real(8) :: F(NF,2), FSIG(NF,2)
integer :: istep_OC
logical :: first_OC
logical :: dbug, WOK, ts_dbug

contains

!---------------------------------------------------------------------
! Redesign:
! Taking into account Parfitt 1994, 1996.  The key idea is that the
! lining cells control the initial placement of an OC.  New OCs are
! located near the apex of the "hemicutting cone", and move laterally
! very little, if at all, after finding their position.
! The rate of arrival of pre-OCs and their clumping and fusion must
! be synchronized with the rate of removel of lining cells, exposing
! the bare bone surface.
! Parfitt's diagrams showing pre-OC ==> OC ==> OB are 2D (X-Z), while
! we need to think about the location of cells on the Y axis as well.
! Assuming that it is correct that the leading face of the trench has
! roughly the shape of a half-cone, we need to put in place the rules
! and mechanisms that will generate such a shape.
!
! Parfitt 1996 data: (based on cortical BMU data: Jaworski1981)
! Total lifetime of an OC nucleus is about 16 days: 3 1/2 days in
! the pre-OC stage, and 12 1/2 days in an OC.  Total BMU team has
! about 80 nuclei (60-100), with a turnover of 8%/day, i.e. 
! 7 new nuclei/day.  BMU lasts 2-8 months.  Team has ~9 OCs (4-16), 
! with ~9 nuclei/OC (3-19). 
! Zallone1984 observed 2-30 nuclei/OC, with most in the range 10-20.
! Roodman1996 says 2-100, most 10-20, and size up to 100 um diameter.
! Trench is ~100-200 um across.  In cortical bone, the BMU travels
! at about 20-40 um/day for a distance of 2-6 mm.
!
! How big is an OC?  Assume we know the size of a monocyte: Rm = 5
! If an OC is made of N monocytes the volume is (4/3)N.Rm^3 = 167.N
! If we assume the shape is hemispherical then 167.N = (2/3).Rc^3
! Therefore Rc = (250N)^(1/3)
! For N = 9, Rc = 2250^(1/3) = 13.1, diameter Dc = 26.2
! In fact the diameter is more like 50, because the shape is more
! flattened than a hemisphere.
!
! Assume that there are 4 OCs across the trench, and that the speed of 
! excavation is the same as for cortical bone, i.e. 30 um/d.  If the
! OCs simply excavate where they sit, and do not move laterally, this
! implies that about 4 new OCs are needed each day, i.e. 36 pre-OCs
! to make new OCs.  In addition the rest of the team will need about 
! 7 new nuclei/d, giving a total pre-OC requirement of 43/d.
!
! Letter to Michael Parfitt
! I've been trying to put the various numbers together to see how rates 
! of monocyte recruitment, nucleus lifetime and trenching rate are related.
! In the absence of good data for cancellous bone, I'm currently using the 
! cortical numbers as a rough first guess.  Let's assume that the BMU team 
! comprises N monocyte nuclei (where we think N is of the order of 80).  
! The average count of nuclei/OC is about 9, according to Jaworski 1981,
! although some others (e.g. Roodman 1996) talk about 10-20 nuclei. 
! It seems reasonable to assume that the number of nuclei in a new OC 
! should be greater than the team average, so for the sake of argument 
! take this number to be 15.  Then if one new OC is created every M days,
! and if the rate of recruitment of monocytes by the existing BMU team 
! members is R (/OC/day), and the death rate of nuclei is 8%/day, then
! under steady-state conditions:
!
! (0.08 - R)N = 15/M
!
! This implies that the maximum rate of creation of new OCs (the minimum value for M), 
! which occurs when R = 0, is one OC created every 15/(0.08N) days.  If N = 80,
! this gives one OC every 2.3 days.  If we take R = 0.04, the interval between 
! the arrival of new OCs becomes 4.6 days.  The rate of OC creation can be 
! related to the average speed of the OCs in the direction of trench advance. 
! I have to make some wild guesses for OC size, trench width and advance rate.
! If we say the new OCs have a diameter of 50 um, and there are three across 
! the active face, and the rate of advance is 25 um/day, then three new OCs
! every two days would seem to be the approximate requirement if the OCs do
! not move forward.  On the other hand, if R=0 and one new OC arrives every
! 2.3 days, i.e. very roughly pushing the leading edge of the BMU forward at
! 50/3 um every 2.3 days, or 50/6.9 = 7 um/day, the team would have to move
! forward at about 18 um/day.  These are very hand-waving numbers, but the
! conclusion seems clear that the OCs in the BMU team do move at not much less
! than the trench advance rate (as implied in Parfitt 1994).  Monocyte recruitment
! by the team makes the required rate of OC forward motion even closer to the
! trench advance rate.
!
! In the limit of large M, no more OCs are created, and R = 0.08, and all
! OCs move forward at 25 um/day.
!
! Assume that the monocyte-attracting signal from a lining cell is suppressed
! when a new OC is created in the vicinity?
!
! Do the lining cells need to disappear (die) as the trench advances?
!
! The addition of a new OC to the team requires a major repositioning,
! to integrate the OC into the MBU.  The forces acting on OCs are two-fold:
! (1) attraction to bone replacement signal (osteocyte)
! (2) mutual attraction - OCs like to be close to each other
! The combination of these two effects determines OC movement.
!
! When a new OC is first created (JOINING), it needs to squeeze itself between the
! front row of the team, i.e. two OCs need to move apart to make space for
! the new OC.  This is a challenge to model.  There probably needs to be
! an active aversion to having the cell surface "exposed", i.e. not 
! within the BMU cluster of OCs.  For cells at the leading edge, this
! is balanced by the greater signal strength there.
! probably need another force:
! (3) force towards the BMU that depends on the degree of exposure of the OC.
! The idea is that the force does not just depend on the distance of the OC
! from the centre of the BMU, it also depends on the proximity of other OCs
! on all sides.  This could be implemented as an amplifying factor on the 
! force of attraction (2).
! Say there are no OCs closer than a threshold distance in the range of
! angles A1 to A2. Then there is an effective force in the direction
! (A1+A2)/2 + Pi, i.e. a vector that bisects A1,A2, but with opposite sense.
! The strength is proportional to the amount by which A2-A1 exceeds Pi.
! In general, the attractive-repulsive force between two OCs with
! radii R1 and R2 needs to be increasingly repulsive for decreasing D,
! D < R1+R2, aapproximately zero near D = R1+R2, increasingly attractive
! as D increases, with a maximum at D = Dmax, then decreasing asymptotically
! to zero as D increases further.
! Note: (3) may not be necessary if bone signal is not turned on beyond
! the leading edge, i.e. lining cells have not lifted off.
!
! With forces (1) and (2), and possibly (3) if needed, we have the makings
! of a procedure for simulating "steady-state" conditions.  What isn't
! clear is how to start the BMU.  Perhaps the first few OCs just clump together.
! How is the direction of advance determined?  As far as we know two BMUs do
! not start moving in opposite directions along the same osteocyte signal 
! fault line.  Perhaps the direction of the capillary is somehow (arbitrarily?)
! determined, then the BMU follows.  In any case, we could ensure that the
! turning on of bone signal (lifting of lining cells) occurs only on the
! forward side of the neighbourhood of an OC, after some time.  The turning
! on could be a continuous variation, up to full turn-on (= 1).
! In the early stages we have to allow for JOINING OCs that inevitably are
! located over bone that has not yet been prepared for them.  In this
! situation the turning-on process may need to be accelerated, as the new
! OC somehow forces the lining cell to lift off.  This should occur only
! after the JOINING OC has been incorporated into the BMU, i.e. -> RESORBING.
! At this point the sites under the OC are turned on.  We just need a criterion 
! for the completion of incorporation into the team (JOINING), possibly related
! to the amount of motion of the OC, or alternately based on elapsed time.
! 
! OC motion could be simulated in two ways: crude and less crude.
! Crudely, all the forces on an OC could be computed, then a
! displacement made that is proportional to the nett force.  In 
! this case OCs would be processed sequentially.
! Less crudely, the full equations of motion for all OCs could be
! solved using Runge-Kutta.  Let's explore this approach.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Need to:
!	create some OCs
!	give them positions
!	compute forces (initially just attraction-repulsion)
!	solve for motion
! Now need to test the arrival of new OCs to join the BMU team.
! For this we can use pre-specified arrival times
!---------------------------------------------------------------------
subroutine test_OCdynamics
!integer :: it

WOK = .not.use_TCP
DELTA_T_OC = DELTA_T
call OCsetup

call PrepareSurface

OC_NV = 4*nclast
if (allocated(OCstate)) deallocate(OCstate)
if (allocated(OCstatep)) deallocate(OCstatep)
allocate(OCstate(OC_NV))
allocate(OCstatep(OC_NV))

! Initialize state, statep
call InitState(OCstate,OCstatep)

istep_OC = 0
first_OC = .true.
dbug = .false.

!do it = 1,5000
!	call simulate_OC_step
!enddo

!call OCforces1(F,NF)
!call OCforces(state,F,NF)
!write(*,*) 'F: ',F(4,1),F(4,2)
!write(*,*) 'OCstate: '
!write(*,'(4f10.6)') OCstate(1:OC_NV)
!write(*,*) 'OCstatep: '
!write(*,'(4f10.6)') OCstatep(1:OC_NV)
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine simulate_OC_step(ierr)
integer :: i, site(3), ierr
type(osteoclast_type), pointer :: pclast
logical :: changed	! this would be set when an OC is either added or removed (or both)
real :: t, dt, sig10, sig01

istep_OC = istep_OC + 1
dt = DELTA_T_OC
t = (istep_OC-1)*DELTA_T_OC
dbug = .false.
ts_dbug = .false.
call OCfsignal(OCstate) 
call OC_mover(OC_NV,t,dt,OCstate,OCstatep,first_OC,ierr)
if (ierr /= 0) return
first_OC = .false.
call UpdateSurface
call MovePits
call Resorb
call LiftSeal
if (mod(istep_OC,10) == 0 .or. istep_OC >= 1000000) then	
!	dbug = .true.
	if (WOK) then
		write(*,'(2i8,6(2x,2f6.2))') istep_OC,nclast,((OCstate(4*i-3),OCstate(4*i-1)),i=1,nclast)
!		call LiftSeal
		call OCfsignal(OCstate)
!		ts_dbug = .true.
		pclast => clast(1)
		site = pclast%cm
		site(1) = site(1) + 1
		sig10 = TotalSignal(pclast,site)
		site(1) = site(1) - 1
		site(3) = site(3) - 1
		sig01 = TotalSignal(pclast,site)
		write(*,*) 'sig10,sig01: ',sig10,sig01
		ts_dbug = .false.
	endif
	dbug = .false.
endif
call CheckJoining(t,changed)
if (changed) then
	write(logmsg,*) 'New OC joined the team: ',istep_OC
	call logger(logmsg)
	call ReinitState(OCstate,OCstatep)
	first_OC = .true.
	changed = .false.
!		pclast =>clast(2)
!		call show(pclast)
!		pclast =>clast(4)
!		call show(pclast)
endif
end subroutine

!--------------------------------------------------------------------------
! Set up a list of OCs, all except the first with entry times in the future.
!--------------------------------------------------------------------------
subroutine OCsetup
type(osteoclast_type), pointer :: pclast
integer :: iclast
real, parameter :: interval = 5	! days

nOClist = 6
allocate(OClist(nOClist))
do iclast = 1,nOClist
	pclast => OClist(iclast)
	pclast%ID = iclast
	pclast%entrytime = (iclast-1)*interval*24*60
	pclast%cm(1) = 3.0
	pclast%cm(2) = NBY + 1
	pclast%cm(3) = NZ/2
	pclast%count = 6
	pclast%radius = (25/DELTA_X)*sqrt(pclast%count/10.)
	pclast%status = QUEUED
	pclast%prevcm = -999
	call MakePits(pclast)
enddo
OClist(1)%status = RESORBING
nclast = 1
clast(nclast) = OClist(1)

end subroutine

!--------------------------------------------------------------------------
subroutine OCsetup1
type(osteoclast_type), pointer :: pclast

nclast = nclast+1
pclast => clast(nclast)
pclast%ID = nclast
pclast%entrytime = 0
pclast%cm(1) = 3.0
pclast%cm(2) = NBY + 1
pclast%cm(3) = NZ/2	+ 2.5
pclast%count = 10
pclast%radius = (25/DELTA_X)*sqrt(pclast%count/10.)
pclast%status = RESORBING
pclast%prevcm = -999
call MakePits(pclast)

nclast = nclast+1
pclast => clast(nclast)
pclast%ID = nclast
pclast%entrytime = 0
pclast%cm(1) = 3.0
pclast%cm(2) = NBY + 1
pclast%cm(3) = NZ/2 - 3
pclast%count = 7
pclast%radius = (25/DELTA_X)*sqrt(pclast%count/10.)
pclast%status = RESORBING
pclast%prevcm = -999
call MakePits(pclast)

nclast = nclast+1
pclast => clast(nclast)
pclast%ID = nclast
pclast%entrytime = 0
pclast%cm(1) = 6.2
pclast%cm(2) = NBY + 1
pclast%cm(3) = NZ/2 + 3.1
pclast%count = 9
pclast%radius = (25/DELTA_X)*sqrt(pclast%count/10.)
pclast%status = RESORBING
pclast%prevcm = -999
call MakePits(pclast)

nclast = nclast+1
pclast => clast(nclast)
pclast%ID = nclast
pclast%entrytime = 0
pclast%cm(1) = 7.5
pclast%cm(2) = NBY + 1
pclast%cm(3) = NZ/2 - 3.8
pclast%count = 12
pclast%radius = (25/DELTA_X)*sqrt(pclast%count/10.)
pclast%status = RESORBING
pclast%prevcm = -999
call MakePits(pclast)
!
nclast = nclast+1
pclast => clast(nclast)
pclast%ID = nclast
pclast%entrytime = 0
pclast%cm(1) = 9.2
pclast%cm(2) = NBY + 1
pclast%cm(3) = NZ/2 + 3 
pclast%count = 15
pclast%radius = (25/DELTA_X)*sqrt(pclast%count/10.)
pclast%status = JOINING
pclast%prevcm = -999
call MakePits(pclast)

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine show(pclast)
type(osteoclast_type), pointer :: pclast
type(surface_type), pointer :: ps
integer :: ipit, x0, z0, x, z
real :: rate

write(*,*) 'clast: ',pclast%ID
x0 = pclast%cm(1) + 0.5
z0 = pclast%cm(3) + 0.5
do ipit = 1,pclast%npit
	x = x0 + pclast%pit(ipit)%delta(1)
	z = z0 + pclast%pit(ipit)%delta(3)
	if (x < 1 .or. x > NX .or. z < 1 .or. z > NZ) cycle
!	ps => surface(x,z)
	write(*,'(4i4,f8.3)') ipit,pclast%pit(ipit)%delta,pclast%pit(ipit)%cover
enddo
end subroutine

!--------------------------------------------------------------------------
! A new OC makes the transition from JOINING to RESORBING when one of two
! conditions is met:
! (1) The displacement in a timestep drop below a threshold
! (2) A specified time limit since the OC's creation is reached.
!--------------------------------------------------------------------------
subroutine CheckJoining(t,joined)
real :: t
logical :: joined
integer :: iclast
real :: dcm(3), d
type(osteoclast_type), pointer :: pclast

joined = .false.
do iclast = 1,nOClist
	pclast => OClist(iclast)
	if (pclast%status /= QUEUED) cycle
	if (pclast%entrytime < t) then
		pclast%status = JOINING
		joined = .true.
		nclast = nclast + 1
		clast(nclast) = OClist(iclast)
		pclast => clast(nclast)
		call SetJoiningLocation(pclast)
		exit
	endif
enddo
do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status /= JOINING) cycle
	dcm = pclast%dcm
	d = sqrt(dcm(1)*dcm(1)+dcm(3)*dcm(3))
!	write(*,'(a,i3,3f10.6)') 'd: ',iclast,pclast%dcm(1),pclast%dcm(3),d 
	if (d < JOIN_DISP_THRESHOLD .or. t-pclast%entrytime > JOIN_TIME_LIMIT) then
		pclast%status = RESORBING
	endif
enddo
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine SetJoiningLocation(pclast0)
type(osteoclast_type), pointer :: pclast0
type(osteoclast_type), pointer :: pclast
real :: xmax, cm(3)
integer :: iclast, imax
integer :: kpar=0

xmax = 0
do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status == DEAD) cycle
	if (pclast%ID == pclast0%ID) cycle
	if (pclast%cm(1) > xmax) then
		xmax = pclast%cm(1)
		imax = iclast
	endif
enddo
cm = clast(imax)%cm
if (cm(3) > NZ/2) then
	pclast0%cm(3) = NZ/2 - 2*(1 + par_uni(kpar))
else
	pclast0%cm(3) = NZ/2 + 2*(1 + par_uni(kpar))
endif
pclast0%cm(1) = cm(1) + 1.5*pclast0%radius
end subroutine

!--------------------------------------------------------------------------
! Set up a test surface, with a target line.
! For initial testing, remove seal wherever there is signal.
!--------------------------------------------------------------------------
subroutine PrepareSurface
integer :: x, z, n0, iclast, dx, dz
real :: c, x0, z0, xx, zz, r0, d, v(2)
type(osteoclast_type), pointer :: pclast

patch%b = LACUNA_B
patch%z0 = NZ/2
do x = 1,NX
	do z = 1,NZ
		surface(x,z)%signal = 0
		surface(x,z)%depth = 0
		surface(x,z)%target_depth = 0
		surface(x,z)%seal = 1
		c = ((z-patch%z0)/patch%b)**2
		if (c < 1) then
			surface(x,z)%target_depth = (1 - c)*MAX_PIT_DEPTH
			surface(x,z)%signal = GetSignal(x,z)
			surface(x,z)%seal = 1
!			if (x <= 2) then	! TEMPORARY measure for startup
!				surface(x,z)%signal = 0	!0.1*surface(x,z)%signal
!				surface(x,z)%seal = 0
!			endif
		endif
	enddo
enddo
do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status /= RESORBING) cycle
	x0 = pclast%cm(1)
	z0 = pclast%cm(3)
	r0 = pclast%radius
	n0 = r0 + 1
	do dx = -n0,n0
		xx = x0 + dx
		if (xx < 1 .or. xx > NX) cycle
		do dz = -n0,n0
			zz = z0 + dz
			if (zz < 1 .or. zz > NZ) cycle
			v(1) = xx - x0
			v(2) = zz - z0
			d = sqrt(v(1)**2 + v(2)**2)
			if (d < r0) then
				x = xx + 0.5
				z = zz + 0.5
				surface(x,z)%seal = 0
			endif
		enddo
	enddo
enddo
	
end subroutine

!--------------------------------------------------------------------------
! %prevcm(:) is set every time site coverage by pits is recomputed.
! The coverage is recomputed when a minimum OC movement has occurred.
!--------------------------------------------------------------------------
subroutine MovePits
type(osteoclast_type), pointer :: pclast
integer :: iclast
real :: dx, dz, d

do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status /= RESORBING) cycle
	dx = pclast%cm(1) - pclast%prevcm(1)
	dz = pclast%cm(3) - pclast%prevcm(3)
	d = sqrt(dx*dx + dz*dz)
	if (d > OC_MOVE_THRESHOLD) then
		pclast%prevcm = pclast%cm
		call SetPitCover(pclast)
	endif
enddo
end subroutine

!--------------------------------------------------------------------------
! Based on the radius of the OC, the list of pits is defined.  Each pit
! has an integer offset from the centre %cm(:), and a value %cover indicating
! what fraction of the site at that offset is currently covered.
! When the OC moves the pit%cover values need to be updated.
!--------------------------------------------------------------------------
subroutine MakePits(pclast)
type(osteoclast_type), pointer :: pclast
type(pit_type) :: temp(200)
integer :: n, dx, dz, k
real :: d, x, z

n = pclast%radius + 3
k = 0
do dx = -n,n
	do dz = -n,n
		d = sqrt(real(dx*dx + dz*dz))
		if (d < pclast%radius + 2.0) then
			k = k+1
			temp(k)%delta(1) = dx
			temp(k)%delta(2) = 0
			temp(k)%delta(3) = dz
			temp(k)%cover = 0
		endif
	enddo
enddo
if (allocated(pclast%pit)) then
	deallocate(pclast%pit)
endif
pclast%npit = k
allocate(pclast%pit(pclast%npit))
do k = 1,pclast%npit
	pclast%pit(k)%delta = temp(k)%delta
	pclast%pit(k)%fraction = 1.0/pclast%npit	! crude first approx.
!	write(*,*) pclast%ID,k,pclast%pit(k)%delta(1),pclast%pit(k)%delta(3)
!	x = pclast%cm(1) + pclast%pit(k)%delta(1)
!	z = pclast%cm(3) + pclast%pit(k)%delta(3)
enddo
call SetPitCover(pclast)
end subroutine

!--------------------------------------------------------------------------
! Whenever the OC moves, or when it is resized, the cover must be recomputed.
! For each pit offset, the fraction within a circle pclast%radius centred at
! pclast%cm must be determined.
! OC movement needs to be more than some threshold distance since the last move.
! OC_MOVE_THRESHOLD
!--------------------------------------------------------------------------
subroutine SetPitCover(pclast)
type(osteoclast_type), pointer :: pclast
integer :: k, site0(3), site(3)
real :: x0, z0, dx, dz, d

!write(*,*) 'SetPitCover: ',pclast%ID
x0 = pclast%cm(1)
z0 = pclast%cm(3)
site0(1) = x0 + 0.5
site0(3) = z0 + 0.5
do k = 1,pclast%npit
	site = site0 + pclast%pit(k)%delta
	if (site(1) < 1 .or. site(1) > NX .or. site(3) < 1 .or. site(1) > NZ) then
		pclast%pit(k)%cover = 0.0
		cycle
	endif	
	dx = site(1) - x0
	dz = site(3) - z0
	d = sqrt(dx*dx + dz*dz)
	if (d < pclast%radius - 0.5) then
		pclast%pit(k)%cover = 1.0
	elseif (d < pclast%radius + 0.5) then
		pclast%pit(k)%cover = pclast%radius + 0.5 - d
	else
		pclast%pit(k)%cover = 0.0
	endif
enddo
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine Resorb
type(osteoclast_type), pointer :: pclast
integer :: iclast, ipit, x0, z0, x, z
real :: rate

do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status /= RESORBING) cycle
	x0 = pclast%cm(1) + 0.5
	z0 = pclast%cm(3) + 0.5
	do ipit = 1,pclast%npit
		x = x0 + pclast%pit(ipit)%delta(1)
		z = z0 + pclast%pit(ipit)%delta(3)
!		if (pclast%pit(ipit)%cover > 0 .and. surface(x,z)%signal > 0 .and. surface(x,z)%seal == 0) then
		if (pclast%pit(ipit)%cover > 0 .and. surface(x,z)%seal == 0) then
			rate = pclast%pit(ipit)%fraction*resorptionRate(pclast%count,pclast%npit)*pclast%pit(ipit)%cover
			surface(x,z)%depth =  surface(x,z)%depth + rate*DELTA_T_OC
			if (SIGNAL_POSITIVE) then
				surface(x,z)%signal = max(0.0,GetSignal(x,z))	
			else
				surface(x,z)%signal = GetSignal(x,z)	
			endif
!			if (iclast == 1 .and. ipit == 23) then
!				write(*,'(2i3,4e12.3)') iclast,ipit,pclast%pit(ipit)%fraction,resorptionRate(pclast%count,pclast%npit),pclast%pit(ipit)%cover,rate
!			endif
		endif
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------
! The surface seal on sites adjacent to an OC pit (with cover > 0) is
! gradually removed.
! This method relies on the set of possible pits,pclast%pit(k)%delta,
! spreading a bit beyond the radius of the OC.
! (see: 	if (d < pclast%radius + 1.5)  in MakePits)
!--------------------------------------------------------------------------
subroutine LiftSeal
integer :: iclast, k
integer :: site0(3), site(3)
type(osteoclast_type), pointer :: pclast
real :: r2, d2
real, parameter :: EDGE = 1.5	! TESTING!!!!!!!!!!!!!!

do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status /= RESORBING) cycle
	site0(1) = pclast%cm(1) + 0.5
	site0(3) = pclast%cm(3) + 0.5
	r2 = (pclast%radius + EDGE)**2
!	if (dbug) write(*,*) 'LiftSeal: ',iclast,pclast%radius,r2
	do k = 1,pclast%npit
		site = site0 + pclast%pit(k)%delta
		if (site(1) < 1 .or. site(1) > NX) cycle
		if (site(3) < 1 .or. site(3) > NZ) cycle
		if (surface(site(1),site(3))%seal == 0) cycle
		d2 = pclast%pit(k)%delta(1)**2 + pclast%pit(k)%delta(3)**2
!		if (dbug) write(*,*) k,site(1),site(3),pclast%pit(k)%delta(1),pclast%pit(k)%delta(3),d2
		if (d2 > r2) cycle
		if (surface(site(1),site(3))%target_depth > 0) then
			!NOTE:  This hard-wires the direction of BMU travel to be +x.  TESTING ONLY  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if ( pclast%pit(k)%delta(1) >= 0) then
!			if (pclast%pit(k)%cover > 0 .or. pclast%pit(k)%delta(1) > 0) then
!			if (pclast%pit(k)%cover > 0) then
				surface(site(1),site(3))%seal = max(0.0,surface(site(1),site(3))%seal - SEAL_REMOVAL_RATE*DELTA_T_OC)
			endif
		endif
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------
! Set up state(:), statep(:) for the first time.
!--------------------------------------------------------------------------
subroutine InitState(state,statep)
real(8) :: state(:), statep(:)
integer :: iclast, indx, k
type(osteoclast_type), pointer :: pclast
indx = 0
do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status == DEAD) cycle
	indx = indx + 1
	OC_to_index(iclast) = indx
	index_to_OC(indx) = iclast
	k = (indx-1)*4 + 1
	state(k) = pclast%cm(1)
	k = k+1
	state(k) = 0
	k = k+1
	state(k) = pclast%cm(3)
	k = k+1
	state(k) = 0
enddo
statep = 0
end subroutine

!--------------------------------------------------------------------------
! When the OC team changes, we need to set up state(:), statep(:) again.
! First deallocate the arrays, then determine the new dimension, then
! allocate again.
! In this simple version, OC velocities and accelerations are all reset
! to zero (statep = 0).  This is not an issue.
!--------------------------------------------------------------------------
subroutine ReinitState(state,statep)
real(8), allocatable :: state(:), statep(:)
integer :: iclast, indx, k
type(osteoclast_type), pointer :: pclast

deallocate(state)
deallocate(statep)

indx = 0
do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status == DEAD) cycle
	indx = indx + 1
enddo
OC_NV = 4*indx
allocate(state(OC_NV))
allocate(statep(OC_NV))

indx = 0
do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status == DEAD) cycle
	indx = indx + 1
	OC_to_index(iclast) = indx
	index_to_OC(indx) = iclast
	k = (indx-1)*4 + 1
	state(k) = pclast%cm(1)
	k = k+1
	state(k) = 0
	k = k+1
	state(k) = pclast%cm(3)
	k = k+1
	state(k) = 0
enddo
statep = 0
end subroutine

!--------------------------------------------------------------------------
! For now assume that the number of active OCs is fixed - none die, none created.
!--------------------------------------------------------------------------
subroutine OC_mover(NV,t,dt,state,statep,first,ierr)
integer :: NV
real :: t, dt
real(8) :: state(:), statep(:)
logical :: first
integer :: ierr
real(8) ::  tstart, tend, relerr, abserr
integer :: flag, iclast, k, indx, x
type(osteoclast_type), pointer :: pclast

ierr = 0
tstart = t
tend = tstart + dt
if (first) then
    call OC_deriv(tstart,state,statep)
	first = .false.
    flag = 1
else
    flag = 2
endif

abserr = sqrt ( epsilon ( abserr ) )
relerr = sqrt ( epsilon ( relerr ) )

call r8_rkf45 ( OC_deriv, OC_NV, state, statep, tstart, tend, relerr, abserr, flag )

do indx = 1,OC_NV/4
	iclast = index_to_OC(indx)
	pclast => clast(iclast)
	k = (indx-1)*4 + 1
	pclast%dcm(1) = state(k) - pclast%cm(1)
	pclast%cm(1) = state(k)
	k = (indx-1)*4 + 3
	pclast%dcm(3) = state(k) - pclast%cm(3)
	pclast%cm(3) = state(k)
	x = pclast%cm(1) + pclast%radius
	if (x >= NX) then	! limit of grid
		call logger('OC has reached the grid limit')
		ierr = 1
		return
	endif
enddo
end subroutine

!--------------------------------------------------------------------------
! Each OC has associated with it 4 variables: x,x',z,z'
! therefore if there are N active OCs, there are 4N variables.
!--------------------------------------------------------------------------
subroutine OC_deriv(t,y,yp)
real(8) :: t, y(*), yp(*)
integer :: indx, nindx, k

call OCforces(y,F,NF)

!do iclast = 1,nclast
!	write(*,'(i4,2f8.4)') iclast,F(iclast,:)
!enddo
nindx = OC_NV/4
do indx = 1,nindx
	k = (indx-1)*4 + 1
	yp(k) = y(k+1)
	k = k+1
	yp(k) = (F(indx,1) - K_OC_DRAG*y(k))/OC_MASS
	k = k+1
	yp(k) = y(k+1)
	k = k+1
	yp(k) = (F(indx,2) - K_OC_DRAG*y(k))/OC_MASS
enddo
!write(*,'(a,8f8.4)') 'y: ',y(1:OC_NV)
!write(*,'(a,8f8.4)') 'yp: ',yp(1:OC_NV)
end subroutine

!---------------------------------------------------------------------
! First, compute forces on OCs from other OCs.  This reflects both the
! tendency of cells to cluster together, and their occupation of space.
! Next need to include attraction towards bone signal.
! Pass with n = dimension of array F(:,2)
!---------------------------------------------------------------------
subroutine OCforces(y,F,n)
real(8) :: y(*)
integer :: n
real(8) :: F(n,2)
integer :: iclast0, iclast1, indx0, indx1, nindx, nm0
integer :: ix, iz, n0, idx, idz, ixx, izz, site0(3), site(3), idxmax, idzmax
real(8) :: x0, z0, x1, z1, r0, r1, v(2), vn(2), d, df, s, fs(2), x, z, ax, az, dfac, dx, dz
real(8) :: ss(0:1,0:1), tot(-1:1,-1:1), tot0, sigmax
type(osteoclast_type), pointer :: pclast0, pclast1
real, parameter :: BS_REACH_FACTOR = 1.5
integer, parameter :: idea = 0

nindx = OC_NV/4
do indx0 = 1,nindx
	iclast0 = index_to_OC(indx0)
!	if (dbug) write(*,*) 'iclast0: ',iclast0
	F(indx0,:) = 0
	pclast0 => clast(iclast0)
	if (pclast0%status == DEAD) then
		write(*,*) 'ERROR: indx0, iclast0: DEAD: ',indx0,iclast0
		stop
	endif
	x0 = y((indx0-1)*4+1)	! pclast0%cm(1)
	z0 = y((indx0-1)*4+3)	! pclast0%cm(3)
	r0 = pclast0%radius
	nm0 = pclast0%count
	
	! Forces from other OCs
	do indx1 = 1,nindx
		if (indx1 == indx0) cycle
		iclast1 = index_to_OC(indx1)
		pclast1 => clast(iclast1)
		if (pclast1%status == DEAD) then
			write(*,*) 'ERROR: indx1, iclast1: DEAD: ',indx1,iclast1
			stop
		endif
		x1 = y((indx1-1)*4+1)
		z1 = y((indx1-1)*4+3)
		r1 = pclast1%radius
		v(1) = x1 - x0
		v(2) = z1 - z0
		d = sqrt(v(1)**2 + v(2)**2)
		if (d > 3*r0) cycle		! otherwise remote OCs exert attraction - probably not a good idea
		vn = v/d
		df = OCattraction(d,r0,r1)
!		F(indx0,:) = F(indx0,:) + K_MUTUAL*df*vn
		if (df > 0) then
			F(indx0,:) = F(indx0,:) + K_MUTUAL_ATTRACT*df*vn
		else
			F(indx0,:) = F(indx0,:) + K_MUTUAL_REPELL*df*vn
		endif
	enddo
	
	! Bone signal forces
	! Need to distinguish between an OC that is JOINING and one that is RESORBING.
	! A JOINING OC must be responsive to signals that are quite remote, in order
	! to be able to squeeze into the BMU team.  This is not an issue for a RESORBING OC.
	! 
	! First attempt:
	! Assume that an OC can sample sites within some radius (e.g. 1.5r) and
	! the attractiveness of each is determined and added vectorially.
	!
	! Second idea:
	! RESORBING OC: two steps
	! (a) Measure totsignal0 under OC
	!     If totsignal0 > threshold
	!         FS = 0
	!     else
	!         Choose best direction of movement.  Consider OC located at a set of
	!         sites nearby, evaluate totsignal at each.  Select maximum = sigmax
	!         if (sigmax < totsignal)
	!             FS = 0
	!         else
	!             FS magnitude ~ (sigmax-totsignal0)
	!             FS direction is site-cm
	!         endif
	!     endif
	!
	! This module can be removed from the OCforces subroutine, and FSIG(:,2),
	! which will be a global array, can be computed less frequently, e.g. every 10 timesteps.
	
	if (idea == 1) then
	fs = 0
	n0 = BS_REACH_FACTOR*r0 + 1
	do idx = -n0,n0
		x = x0 + idx
		if (x < 1 .or. x > NX) cycle
		
		do idz = -n0,n0
			z = z0 + idz
			if (z < 1 .or. z > NZ) cycle
!			if (dx == 0 .and. dz == 0) cycle
			v(1) = x - x0
			v(2) = z - z0
			d = sqrt(v(1)**2 + v(2)**2)
			if (d > n0) cycle
			if (d < 1) then
				dfac = 1
			else
				dfac = 1/d
			endif
			ax = x - int(x)
			az = z - int(z)
			vn = v/d
			do ixx = 0,1
				do izz = 0,1
!					if (surface(ix+ixx,iz+izz)%seal == 1) then
!						ss(ixx,izz) = 0
!					else
!						ss(ixx,izz) = surface(ix+ixx,iz+izz)%signal*(1-surface(ix+ixx,iz+izz)%seal)
!					endif
					if (surface(ix+ixx,iz+izz)%seal /= 0) then
						ss(ixx,izz) = 0
					else
						ss(ixx,izz) = surface(ix+ixx,iz+izz)%signal
					endif
				enddo
			enddo
			ix = x
			iz = z
!			s = (1-ax)*(1-az)*surface(ix,iz)%signal*(1 - surface(ix,iz)%seal) &
!			  + ax*(1-az)*surface(ix+1,iz)%signal*(1 - surface(ix+1,iz)%seal) &
!			  + (1-ax)*az*surface(ix,iz+1)%signal*(1 - surface(ix,iz+1)%seal) &
!			  + ax*az*surface(ix+1,iz+1)%signal*(1 - surface(ix+1,iz+1)%seal)
			s = (1-ax)*(1-az)*ss(0,0) + ax*(1-az)*ss(1,0) + (1-ax)*az*ss(0,1) + ax*az*ss(1,1)
			fs = fs + K_SIGNAL*nm0*s*vn/dfac
!			if (dbug .and. iclast0 == 2) then
!				write(*,*) 
!			if (dbug .and. dz == 2) write(*,'(4f12.6)') x,z,s,s*vn(1)		!surface(x,z)%signal,surface(x,z)%seal
		enddo
	enddo
	elseif (idea == 2) then
		if (x0 < 1 .or. z0 < 1) then
			write(*,*) 'ERROR: OCforces: x0 or z0 < 0.5'
			stop
		endif
		if (x0 > NX .or. z0 > NZ) then
			write(*,*) 'ERROR: OCforces: x0 or z0 < 0.5'
			stop
		endif
		tot = -1
		site0(1) = x0 + 0.5
		site0(2) = NBY + 1
		site0(3) = z0 + 0.5
		tot(0,0) = TotalSignal(pclast0,site0)
		site = site0
		dx = x0 - site0(1)
		dz = z0 - site0(3)
		if (dx >= 0) then
			idx = 1
			ax = dx
		else
			idx = -1
			ax = -dx
		endif
		site(1) = site(1) + idx
		tot(idx,0) = TotalSignal(pclast0,site)
		site(1) = site0(1)
		if (dz >= 0) then
			idz = 1
			az = dz
		else
			idz = -1
			az = -dz
		endif
		site(3) = site(3) + idz
		tot(0,idz) = TotalSignal(pclast0,site)
		site(1) = site0(1) + idx
		site(3) = site0(3) + idz
		tot(idx,idz) = TotalSignal(pclast0,site)
		tot0 = (1-ax)*(1-az)*tot(0,0) + ax*(1-az)*tot(idx,0) &
		     + (1-ax)*az*tot(0,idz)   + ax*az*tot(idx,idz)
		if (tot0 > 0.4*pclast0%npit) then	! arbitrary threshold
!			if (dbug) then
!				write(*,*) 'ax,az: ',idx,idz,ax,az
!				write(*,*) 'tot0: ',tot0,0.4*pclast0%npit,tot(idx,idz)
!			endif
			fs = 0
		else
			sigmax = 0
			do idx = -1,1
				do idz = -1,1
					if (tot(idx,idz) < 0) then
						site(1) = site0(1) + idx
						site(2) = NBY + 1
						site(3) = site0(3) + idz
						tot(idx,idz) = TotalSignal(pclast0,site)
					endif
					if (tot(idx,idz) > sigmax) then
						idxmax = idx
						idzmax = idz
						sigmax = tot(idx,idz)
					endif
				enddo
			enddo
			if (sigmax > tot0) then
				v(1) = site0(1) + idxmax - x0
				v(2) = site0(3) + idzmax - z0
				d = sqrt(v(1)**2 + v(2)**2)
				if (d == 0) then
					fs = 0
				else
					vn = v/d
					fs = K_SIGNAL*nm0*(sigmax-tot0)*vn
				endif
			else
				fs = 0
			endif
		endif
!		if (dbug) then
!			write(*,*) 'fs: ',fs
!			stop
!		endif
	endif
!	if (dbug .and. iclast0 == 1) then
!		write(*,'(a,i3,4f8.4)') 'iclast0, F_S, F_OC: ',iclast0,fs,F(indx0,:)
!	endif
	pclast0%foc = F(indx0,:)
	if (idea == 0) then
		F(indx0,:) = F(indx0,:) + FSIG(indx0,:)
	else
		F(indx0,:) = F(indx0,:) + fs
	endif
	pclast0%ftot = F(indx0,:)
enddo
F = SPEEDUP*F
end subroutine

!---------------------------------------------------------------------
! Bone signal forces
! Need to distinguish between an OC that is JOINING and one that is RESORBING.
! A JOINING OC must be responsive to signals that are quite remote, in order
! to be able to squeeze into the BMU team.  This is not an issue for a RESORBING OC.
!
! RESORBING OC: two steps
! (a) Measure totsignal0 under OC
!     If totsignal0 > threshold
!         FS = 0
!     else
!         Choose best direction of movement.  Consider OC located at a set of
!         sites nearby, evaluate totsignal at each.  Select maximum = sigmax
!         if (sigmax < totsignal)
!             FS = 0
!         else
!             FS magnitude ~ (sigmax-totsignal0)
!             FS direction is site-cm
!         endif
!     endif
!
! JOINING OC:
! Needs to push its way in.
!
! This module can be removed from the OCforces subroutine, and FSIG(:,2),
! which will be a global array, can be computed less frequently, e.g. every 10 timesteps.
!---------------------------------------------------------------------
subroutine OCfsignal(y)
real(8) :: y(*)
integer :: iclast0, indx0, nindx, nm0
integer :: ix, iz, n0, idx, idz, ixx, izz, site0(3), site(3), idxmax, idzmax
real(8) :: x0, z0, r0, v(2), vn(2), d, fs(2), x, z, ax, az, dx, dz
integer, parameter :: Q=2
real(8) :: tot(-Q:Q,-Q:Q), tot0, sigmax
type(osteoclast_type), pointer :: pclast0

nindx = OC_NV/4
do indx0 = 1,nindx
	iclast0 = index_to_OC(indx0)
	pclast0 => clast(iclast0)
	if (pclast0%status == DEAD) then
		write(logmsg,*) 'ERROR: indx0, iclast0: DEAD: ',indx0,iclast0
		call logger(logmsg)
		stop
	endif
	x0 = y((indx0-1)*4+1)	! pclast0%cm(1)
	z0 = y((indx0-1)*4+3)	! pclast0%cm(3)
	r0 = pclast0%radius
	nm0 = pclast0%count
	if (x0 < 1 .or. z0 < 1) then
		write(logmsg,*) 'ERROR: OCfsignal: x0 or z0 < 0.5'
		call logger(logmsg)
!		stop
	endif
	if (x0 > NX .or. z0 > NZ) then
		write(logmsg,*) 'ERROR: OCforces: x0 or z0 < 0.5'
		call logger(logmsg)
!		stop
	endif
	
	fs = 0
	if (pclast0%status == RESORBING .or. pclast0%status == JOINING) then
		tot = -1
		site0(1) = x0 + 0.5
		site0(2) = NBY + 1
		site0(3) = z0 + 0.5
		tot(0,0) = TotalSignal(pclast0,site0)
		site = site0
		dx = x0 - site0(1)
		dz = z0 - site0(3)
		if (dx >= 0) then
			idx = 1
			ax = dx
		else
			idx = -1
			ax = -dx
		endif
		site(1) = site(1) + idx
		tot(idx,0) = TotalSignal(pclast0,site)
		site(1) = site0(1)
		if (dz >= 0) then
			idz = 1
			az = dz
		else
			idz = -1
			az = -dz
		endif
		site(3) = site(3) + idz
		tot(0,idz) = TotalSignal(pclast0,site)
		site(1) = site0(1) + idx
		site(3) = site0(3) + idz
		tot(idx,idz) = TotalSignal(pclast0,site)
		tot0 = (1-ax)*(1-az)*tot(0,0) + ax*(1-az)*tot(idx,0) &
			 + (1-ax)*az*tot(0,idz)   + ax*az*tot(idx,idz)
		if (pclast0%status == RESORBING .and. tot0 > 0.4*pclast0%npit) then	! arbitrary threshold
!			if (dbug) then
!				write(*,*) 'ax,az: ',idx,idz,ax,az
!				write(*,*) 'tot0: ',tot0,0.4*pclast0%npit,tot(idx,idz)
!			endif
			fs = 0
		else
			sigmax = 0
			do idx = -Q,Q
				do idz = -Q,Q
					if (tot(idx,idz) < 0) then
						site(1) = site0(1) + idx
						site(2) = NBY + 1
						site(3) = site0(3) + idz
						tot(idx,idz) = TotalSignal(pclast0,site)
					endif
					if (tot(idx,idz) > sigmax) then
						idxmax = idx
						idzmax = idz
						sigmax = tot(idx,idz)
					endif
				enddo
			enddo
			if (sigmax > tot0) then
				v(1) = site0(1) + idxmax - x0
				v(2) = site0(3) + idzmax - z0
				d = sqrt(v(1)**2 + v(2)**2)
				if (d == 0) then
					fs = 0
				else
					vn = v/d
!					if (dbug) write(*,*) 'vn: ',vn
					fs = K_SIGNAL*nm0*(sigmax-tot0)*vn
				endif
			else
				fs = 0
			endif
			if (WOK .and. dbug) then
				write(*,'(a,4f10.6)') 'sigmax, tot0, fs: ',sigmax,tot0,fs
				do idz = -Q,Q
					write(*,'(5f10.6)') (tot(idx,idz),idx=-Q,Q)
				enddo
			endif
			if (pclast0%status == JOINING) then
				fs = K_JOINING*fs
			endif
		endif
		pclast0%fsig = fs

	endif
!	if (dbug) then
!		write(*,*) 'fs: ',fs
!	endif
	FSIG(indx0,:) = fs
enddo
!write(*,'(a,4(2x,2e12.3))') 'FS: ',(FSIG(iclast0,:),iclast0=1,nclast)	
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine OCforces1(F,n)
integer :: n
real(8) :: F(n,2)
integer :: iclast0, iclast1
real(8) :: x0, z0, x1, z1, r0, r1, v(2), vn(2), d, df
type(osteoclast_type), pointer :: pclast0, pclast1

do iclast0 = 1,nclast
	F(iclast0,:) = 0
	pclast0 => clast(iclast0)
	if (pclast0%status == DEAD) cycle
	x0 = pclast0%cm(1)
	z0 = pclast0%cm(3)
	r0 = pclast0%radius
	
	! Forces from other OCs
	do iclast1 = 1,nclast
		if (iclast1 == iclast0) cycle
		pclast1 => clast(iclast1)
		if (pclast1%status == DEAD) cycle
		x1 = pclast1%cm(1)
		z1 = pclast1%cm(3)
		r1 = pclast1%radius
		v(1) = x1 - x0
		v(2) = z1 - z0
		d = sqrt(v(1)**2 + v(2)**2)
		vn = v/d
		write(*,*) 'd,r0,r1,vn: ',d,r0,r1,vn
		df = OCattraction(d,r0,r1)
		write(*,*) 'df*vn: ',df*vn
		F(iclast0,:) = F(iclast0,:) + df*vn
		write(*,*) 'F: ',F(iclast0,:)
	enddo
enddo
	
end subroutine


!---------------------------------------------------------------------
! The desired shape of the attraction-repulsion function is not obvious.
! Clearly f -> 0 as d -> inf., and f -> -inf. as d -> K1(r0+r1) (K1 < 1)
! Near d = r0+r1, say at d = K2(r0+r1), we want f = 0, and the slope
! near the axis crossing should be small.
! The attractive force will rise to a maximum, at say d = K3(r0+r1),
! then decrease (in an inverse square way?) with increasing d.
! As in attraction-repulsion.xlsx:
! Using d = r0+r1 as the axis-crossing point (simple case to start)
! and setting x = d/(r0+r1),
! For x <= 1
!	f(x) = -b/(x-1+a) + b/a
!   Note: x < 1-a => ERROR (attraction)
! for x >= 1
!	f(x) = h.exp(-k(x-1)).(exp(g(x-1))-1)/(exp(g(x-1))+1)
! The slopes are matched (= s) at x=1 when:
!	b = s.a^2
!	g = (2s)/h
! Reasonable curve is obtained with:
!	s = 0.1
!	a = 0.5
!	h = 20
!	k = 3
! Note that the slope is still a bit steep near x=1
!---------------------------------------------------------------------
real(8) function OCattraction(d,r0,r1)
real(8) :: d, r0, r1, x
real, parameter :: s = 0.1
real, parameter :: a = 0.5
real, parameter :: h = 20
real, parameter :: k = 3
real(8) :: b, g

b = s*a**2
g = (2*s)/h
x = d/(r0+r1)
!if (dbug) write(*,*) 'x: ',d,r0+r1,x
if (x <= 1) then
	if (x < 1-a) then
		write(logmsg,*) 'ERROR: in OCattraction: x < 1-a: ',x,1-a
		call logger(logmsg)
		x = 1.01*(1-a)
!		stop
	endif
	OCattraction = -b/(x-1+a) + b/a
!	if (dbug) write(*,*) -b/(x-1+a),b/a
else
	OCattraction = h*exp(-k*(x-1))*(exp(g*(x-1))-1)/(exp(g*(x-1))+1)
endif
end function

!---------------------------------------------------------------------
! We now take account of the surface seal, and of the presence of
! another OC.
! This method computes total signal based on OC centred on a site.
! Note that %cover is not used, and that the total extent of %delta(:)
! is summed.
!---------------------------------------------------------------------
real function TotalSignal(pclast,site)
type(osteoclast_type), pointer :: pclast
integer :: site(3)
integer :: k, x, z, ic, dx, dz
real :: sig, r2

if (ts_dbug) write(*,*) 'TotalSignal: site: ',site
r2 = pclast%radius**2
TotalSignal = 0
do k = 1,pclast%npit
	dx = pclast%pit(k)%delta(1)
	dz = pclast%pit(k)%delta(3)
	x = site(1) + dx
	z = site(3) + dz
	if (x < 1 .or. x > NX) cycle
	if (z < 1 .or. z > NZ) cycle
	if (dx*dx+dz*dz > r2) cycle
	ic = surface(x,z)%iclast
	if (pclast%status == RESORBING .and. ic /= 0 .and. ic /= pclast%ID) cycle
	if (SEAL_SWITCH) then
		if (surface(x,z)%seal /= 0) then
			sig = 0
		else
			sig = surface(x,z)%signal
		endif
	else
		sig = surface(x,z)%signal*(1 - surface(x,z)%seal)
	endif
	if (ts_dbug) write(*,'(3i4,3f8.4)') k,x,z,surface(x,z)%signal,surface(x,z)%seal,sig
	TotalSignal = TotalSignal + sig
!	write(*,*) k,x,z,surface(x,z)%signal,surface(x,z)%seal
enddo
TotalSignal = max(0.0,TotalSignal)
!if (dbug) write(*,*) 'TotalSignal: ',TotalSignal
end function

!---------------------------------------------------------------------
! The signal in the OC foorprint is summed, but set to -1 if it
! overlaps another OC footprint.
!---------------------------------------------------------------------
real function ClearSignal(pclast,site)
type(osteoclast_type), pointer :: pclast
integer :: site(3)
integer :: k, x, z, iclast, k1, x1, z1, site1(3)

ClearSignal = 0
do k = 1,pclast%npit
	x = site(1) + pclast%pit(k)%delta(1)
	z = site(3) + pclast%pit(k)%delta(3)
	do iclast = 1,nclast
		if (clast(iclast)%ID == pclast%ID) cycle
		if (clast(iclast)%status /= RESORBING) cycle
		site1 = clast(iclast)%site		! site of another resorbing OC
		do k1 = 1,clast(iclast)%npit
			x1 = site1(1) + clast(iclast)%pit(k1)%delta(1)
			z1 = site1(3) + clast(iclast)%pit(k1)%delta(3)
			if (x == x1 .and. z == z1) then	! footprints overlap
				ClearSignal = -1
				return
			endif
		enddo
	enddo
	ClearSignal = ClearSignal + surface(x,z)%signal
enddo
end function

!---------------------------------------------------------------------
! An OC moves when the total signal falls below a threshold level, and
! it moves in the direction that maximizes total signal.  The only
! constraint is that it must avoid other OCs.
! (Note that when an OC is first created it moves to the best signal
! location immediately, without beginning to resorb bone.  This tends
! to ensure that a new OC moves to the mine face.)
! The osteoclast is moved one lattice jump, and the locations of active
! pit sites are recomputed, together with their resorption rates.
! Chemotaxis
! When osteoclast motion is influenced by the bone signal, the jump 
! probabilities depend on the relative signal strengths of the 
! neighbour sites.
! There are two cases:
! (1) OC is on or near sites needing excavation (signal > 0)
! (2) OC is not near sites needing excavation.  In this case the OC 
!     must first move towards sites with signal.  It's probably safe to
!     assume that the OC has a range for detection of signal.
!
! Result values (res):
!     0 = prefer to stay, continue resorbing
!     1 = prefer or need to move, move possible
!     2 = blocked, will die if blocked for too long
!---------------------------------------------------------------------
subroutine MoveClast(pclast,res)
type(osteoclast_type), pointer :: pclast
integer :: res
integer :: ipit, iy, x, y, z, site(3), i, imono, kdir, dx, dz, dirmax, iclast
integer :: ddir, nsig, ocsite(3), kmax, targetsite(3), kmin
real :: prob(0:4) = (/ 0.5, 0.15, 0.07, 0.025, 0.01 /)
integer :: jump(3), lastjump(3), kpar=0
real :: d, bf, v(3), proj, size0, size1, djump, dsum, depth
real :: totsig(0:8), siglim, sig, patchtot, amp(8), cosa, amax, tnow, dmin, sigmax
real(8) :: psum, pmax, R, dp, p(8), ptemp(8)
!real :: p(3) = (/0.25,0.5,0.25/)	! arbitrary, interim
logical :: covered, bdryhit, nearclast, freshbone, possible(8)
type(osteoclast_type), pointer :: pclast1
integer, save :: count = 0
logical :: dbug
real, parameter :: SIGNAL_EXCESS = 1.05

if (pclast%ID == -1) then
	dbug = .true.
else
	dbug = .false.
endif

tnow = istep*DELTA_T
res = 0
! First choose a direction.  Quantify the new bone available in each direction
! at a distance approx equal to the long axis dimension of the osteoclast, i.e.
! sqrt(count).

ocsite = pclast%site
if (pclast%status == MOVING) then	! special case, OC moving fast, possibly over other cells
	targetsite = pclast%targetsite
	if (ocsite(1) == targetsite(1) .and. ocsite(3) == targetsite(3)) then	! OC has reached the target site
		pclast%status = RESORBING
		res = 0
		return
	endif
	dmin = 1.0e10
	do kdir = 1,8
		jump = dir2D(:,kdir)
		site = ocsite + jump
		v = targetsite - site
		v(2) = 0
		d = sqrt(dot_product(v,v))
		if (d < dmin) then
			dmin = d
			kmin = kdir
		endif
	enddo
	call ocmove(pclast,kmin)
	if (dmin == 0) then
		pclast%status = RESORBING
	endif
	res = 1
	return
endif

totsig(0) = TotalSignal(pclast,ocsite)
!write(*,*) 'totsig(0): ',istep,totsig(0)/pclast%npit
if (totsig(0)/pclast%npit > OC_SIGNAL_THRESHOLD) then
	! no need to move
	res = 0
	return
endif
!write(*,*) 'totsig(0): ',totsig(0)/pclast%npit
! Explore possible directions of motion.
possible = .true.
nsig = 0
size0 = pclast%radius
kmax = 0
sigmax = 0
!if (dbug) write(*,'(a,i4,4f8.1)') 'move_clast: ',pclast%ID,pclast%cm,size0
do kdir = 1,8
	jump = dir2D(:,kdir)
	site = ocsite + jump
	totsig(kdir) = TotalSignal(pclast,site)
!	write(*,*) kdir,totsig(kdir)/pclast%npit
!	if (dbug) write(*,*) 'kdir,i,site: ',kdir,i,site
	bdryhit = .false.
	if ((site(1) <= 1 .or. site(1) >= NX) .or. (site(3) <= 1 .or. site(3) >= NZ)) then
		bdryhit = .true.
		possible(kdir) = .false.
		cycle
	endif
	! Check for nearby osteoclasts.
	! Treat osteoclast footprint as a circle with radius = clast_size
	nearclast = .false.
	do iclast = 1,nclast
		pclast1 => clast(iclast)
		if (pclast%ID == pclast1%ID) cycle
		if (pclast1%status == DEAD) cycle
		size1 = pclast1%radius
		v = pclast1%site - pclast%site
		d = sqrt(dot_product(v,v))
		if (d < (size0 + size1)) then	! we don't want to move in the direction of v
			proj = dot_product(v,real(jump))
			if (proj > 0) then
!				write(*,*) 'Too close: ', iclast
				nearclast = .true.
				exit
!				possible(kdir) = .false.
			endif
		endif
	enddo	
	if (nearclast) then
		if (dbug) write(*,*) 'nearclast: ',kdir,iclast,size0,size1
		possible(kdir) = .false.
		cycle
	endif
	if (totsig(kdir) > sigmax) then
		kmax = kdir
		sigmax = totsig(kdir)
	endif
enddo

if (kmax == 0) then	
	! OC is not near sites needing excavation.
	! Blocked
	res = 2
else
	! make a move if pushed by an OC, or if there is a better site and dwell time is exceeded
	if (nearclast .or. (tnow > pclast%movetime .and. totsig(kmax) > totsig(0))) then
		call ocmove(pclast,kmax)
		res = 1
	else
		res = 0
	endif
endif
end subroutine

!---------------------------------------------------------------------
! An OC wants to move but cannot.  Need to find the best OC move(s) to 
! improve the OC distribution. 
!---------------------------------------------------------------------
subroutine unblocker
integer :: iclast, icmax, kdir, kmax
type(osteoclast_type), pointer :: pclast
real :: dvmax, tnow
real :: dV(MAX_CLAST,8)		! change in value associated with the configuration
							! that results from a move of clast(iclast) in direction kdir
tnow = istep*DELTA_T
dvmax = -1.0e10
icmax = 0
kmax = 0
!write(*,*) 'unblocker: ',nclast
do iclast = 1,nclast
	if (clast(iclast)%status == DEAD) cycle
	do kdir = 1,8
		dV(iclast,kdir) = GetValueChange(iclast,kdir)
!		write(logmsg,*) 'unblocker: ',iclast,kdir,dV(iclast,kdir)
!		call logger(logmsg)
		if (dV(iclast,kdir) > dvmax) then
			dvmax = dV(iclast,kdir)
			icmax = iclast
			kmax = kdir
		endif
	enddo
enddo
!write(logmsg,*) 'best move: ',icmax,kmax,dir2D(:,kmax)
!call logger(logmsg)
pclast => clast(icmax)
if (tnow > pclast%movetime) then
	call ocmove(pclast,kmax)
endif
end subroutine

!---------------------------------------------------------------------
! The change in value associated with an OC move is the increase in
! attractiveness - the increase in proximity cost.  A deficit of signal
! has negative attractiveness.  A positive value change is a good thing.
!---------------------------------------------------------------------
real function GetValueChange(icmove,kdir)
integer :: icmove, kdir
integer :: iclast, site1(3), site2(3)
real :: d1, d2, v(3), r0, r1, val1, val2, dval_repulsion, dval_attraction
real, parameter :: Krepulsion = 1.0

site1 = clast(icmove)%site
site2 = site1 + dir2D(:,kdir)
r0 = clast(icmove)%radius
dval_repulsion = 0
do iclast = 1,nclast
	if (clast(iclast)%status == DEAD) cycle
	if (iclast == icmove) cycle
	r1 = clast(iclast)%radius
	v = clast(iclast)%site - site1
	d1 = sqrt(dot_product(v,v))
	v = clast(iclast)%site - site2
	d2 = sqrt(dot_product(v,v))
	if (d1 < (r0 + r1)) then
		val1 = (r0 + r1 - d1)**2
	else
		val1 = 0
	endif	
	if (d2 < (r0 + r1)) then
		val2 = (r0 + r1 - d2)**2
	else
		val2 = 0
	endif
	dval_repulsion = dval_repulsion + Krepulsion*(val1 - val2)
enddo
dval_attraction = attractiveness(icmove,site2) - attractiveness(icmove,site1)
!write(*,*) 'attract, repulse: ',dval_attraction,dval_repulsion
GetvalueChange = dval_attraction + dval_repulsion
end function

!---------------------------------------------------------------------
! Signal above threshold makes a site attractive, below is penalized.
!---------------------------------------------------------------------
real function attractiveness(iclast,site)
integer :: iclast, site(3)
type(osteoclast_type), pointer :: pclast
integer :: k, x, z, ic
real :: sig, total
real, parameter :: Kpenalty = 5

pclast => clast(iclast)
total = 0
do k = 1,pclast%npit
	x = site(1) + pclast%pit(k)%delta(1)
	z = site(3) + pclast%pit(k)%delta(3)
	ic = surface(x,z)%iclast
	if (ic /= 0 .and. ic /= iclast) cycle
	sig = GetSignal(x,z) 
!	write(*,*) 'attractiveness: ',k,x,z,sig
	sig = sig - OC_SIGNAL_THRESHOLD
	if (sig > 0) then
		total = total + sig*(1 - surface(x,z)%seal)
	else
		total = total + Kpenalty*sig
	endif
enddo
attractiveness = total
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine ChooseTarget(pclast,targetsite)
type(osteoclast_type), pointer :: pclast
integer :: targetsite(3)
integer :: ocsite0(3), ocsite1(3), x, z
real :: v(3), d, sig, amp, amax

ocsite0 = pclast%site
targetsite = 0
amax = 0
do x = 1,NX
	do z = 1,NZ
		if (surface(x,z)%target_depth > 0) then	! consider this as a possible OC target site
			ocsite1 = (/ x, NBY+1, z /)
			v = ocsite1 - ocsite0
			v(2) = 0
			sig = ClearSignal(pclast,ocsite1)
			if (sig <= 0) cycle
			d = sqrt(dot_product(v,v))
			if (d <= OC_SIGNAL_SENSING_RANGE) then
				amp = sig/d
				if (amp > amax) then
					amax = amp
					targetsite = ocsite1
				endif
			endif
		endif
	enddo
enddo

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine ocmove(pclast,kdir)
type(osteoclast_type), pointer :: pclast
integer :: kdir
integer :: jump(3), site(3), i, imono
real :: tnow

tnow = istep*DELTA_T
jump = dir2D(:,kdir)
if (pclast%ID == -1) then
	write(logmsg,'(a,5i4)') 'clast moves: kdir: ',pclast%ID,kdir,jump
	call logger(logmsg)
endif
pclast%site = pclast%site + jump
pclast%lastdir = kdir
pclast%movetime = tnow + CLAST_DWELL_TIME
!do i = 1,pclast%count
!	imono = pclast%mono(i)
!	site = mono(imono)%site
!	occupancy(site(1),site(2),site(3))%indx = 0
!	site = site + jump
!	mono(imono)%site = site
!	occupancy(site(1),site(2),site(3))%indx = imono
!enddo
call UpdateSurface
end subroutine

!------------------------------------------------------------------------------------------------
! Refresh the locations of OCs on the surface. 
!------------------------------------------------------------------------------------------------
subroutine UpdateSurface
type(osteoclast_type), pointer :: pclast
integer :: iclast, i, site(3), x, z
real :: r2, d2

surface%iclast = 0
do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status == DEAD) cycle 
	r2 = pclast%radius**2
	pclast%site = pclast%cm + 0.5
	do i = 1,pclast%npit
		d2 = pclast%pit(i)%delta(1)**2 + pclast%pit(i)%delta(3)**2
		if (d2 > r2) cycle
		site = pclast%site + pclast%pit(i)%delta
		x = site(1)
		z = site(3)
		if (surface(x,z)%iclast == 0) then
			surface(x,z)%iclast = iclast
		elseif (surface(x,z)%iclast < 100) then
			surface(x,z)%iclast = 100*surface(x,z)%iclast + iclast
		endif
	enddo
enddo
end subroutine

!------------------------------------------------------------------------------------------------
! An OB will try to stay away from other OBs, and from OCs.
! An OB will try to locate to a high-signal region.
! An OB may have a clump attached, which must move with it.
! An attached clump slows down the OB motility.
! When an OB is unattached, it may travel over the top of OCs, to get to a better bone site.  
! In this state (CROSSING) a clump cannot begin to form on it (nearOB() will return 0).
! To avoid clump interference, a hard limit is placed on OB separation: the distance between
! OBs must not be less than 7 grids.  (Note that clumps attached to different OBs cannot merge.)
!------------------------------------------------------------------------------------------------
subroutine MoveBlast(pblast,res)
type(osteoblast_type), pointer :: pblast
integer :: res

end subroutine

!------------------------------------------------------------------------------------------------
! We need a way to determine whether a monocyte has moved close enough to the bone surface to
! make contact with an OB - needed to initiate clustering.
! Base the criterion on nearness to a signalling site.
!------------------------------------------------------------------------------------------------
integer function NearOB(pmono)
type(monocyte_type), pointer :: pmono
real :: v(3), d2, d2lim, d2min
integer :: msite(3), bsite(3), iblast, ibmin

d2lim = (OB_REACH/DELTA_X)**2
msite = pmono%site
d2min = 1.0e10
ibmin = 0
do iblast = 1,nblast
	if (blast(iblast)%status /= ALIVE) cycle
	bsite = blast(iblast)%site
	v = msite - bsite
	d2 = dot_product(v,v)
	if (d2 < d2min) then
		d2min = d2
		ibmin = iblast
	endif
enddo
if (d2min < d2lim) then
	NearOB = ibmin
else
	NearOB = 0
endif
end function

!---------------------------------------------------------------------
! Currently monocytes move around while in the marrow, and may pass
! into the blood, effectively removed from further consideration.
!---------------------------------------------------------------------
subroutine MonoMover
type (monocyte_type), pointer :: pmono
integer :: kcell, status, iblast
logical :: go
integer :: kpar = 0
integer :: kdbug = -300
integer :: site(3)
real :: tnow

tnow = istep*DELTA_T
do kcell = 1,nmono
	pmono => mono(kcell)
	status = pmono%status
	if (status == LEFT .or. status == DEAD) cycle
!	if (pmono%site(2) == NBY+1 .and. status < FUSED) then
!		write(*,*) 'y = NBY+1: ',istep,pmono%ID,pmono%site,status
!		stop
!	endif
	if (status == CROSSING) then
		if (tnow >= pmono%exittime) then
			pmono%status = LEFT
			pmono%region = BLOOD
			site = pmono%site
			occupancy(site(1),site(2),site(3))%indx = 0
			mono_cnt = mono_cnt - 1
			nleft = nleft + 1
!			write(*,'(a,2i6,2f6.3,i6)') 'monocyte leaves: ',istep,kcell,mono(kcell)%S1P1,mono_cnt
		endif
	elseif (status >= MOTILE .and. pmono%iclump == 0) then	! interim criterion
		call MonoJumper(kcell,go,kpar)
	endif
	iblast = NearOB(pmono)
	if (kcell == kdbug) then
		write(nflog,'(2i6,6i4)') istep,kcell,pmono%site,iblast,status
	endif
	if (iblast == 0) cycle
	if (blast(iblast)%status == CROSSING) cycle
!	if (pmono%stickiness > 0 .and. pmono%status < FUSED .and. NearOB(pmono)) then
	if (pmono%stickiness > 0 .and. pmono%status < FUSED) then
		call sticker(kcell,iblast)
	endif
enddo
	
end subroutine

!--------------------------------------------------------------------------------
! When a monocyte is within the field of influence of a signal site, the
! jump probabilities are modified to increase the prob of jumps towards
! the site.  When close enough to the signal site (intensity above a threshold)
! the jump probabilities are determined solely by the relative intensities.
! This is to ensure that monocytes cluster.
! Note: now CXCL12 chemotaxis does not depend on S
!--------------------------------------------------------------------------------
subroutine MonoJumper(kcell,go,kpar)
integer :: kpar,kcell
logical :: go
type (monocyte_type), pointer :: cell
integer :: site1(3),site2(3)
integer :: region, kcell2
integer :: irel,dir1,lastdir1,status
integer :: savesite2(3,26), jmpdir(26)
real(8) :: psum, p(26), R, wS1P(27), g(3), gamp, S1Pfactor, wCXCL12(27), CXCL12factor
real :: f0, f, motility_factor
logical :: free, cross, field

!write(*,*) 'istep,kcell: ',istep,kcell
cell => mono(kcell)
status = cell%status
if (status == MOTILE) then
	motility_factor = 1
elseif (status == CHEMOTACTIC) then
	motility_factor = 0.6
elseif (status == STICKY) then
	motility_factor = 0.4
endif
!if (istep - cell%lastmovestep > 1050) then
!	write(logmsg,*) 'No move: ',cell%ID,istep,kcell,istep - cell%lastmovestep,cell%site
!	call logger(logmsg)
!	stop
!elseif (istep - cell%lastmovestep > 1000) then
!	write(logmsg,*) 'No move: ',cell%ID,istep,kcell,istep - cell%lastmovestep,cell%site
!	call logger(logmsg)
!endif
site1 = cell%site
go = .false.
field = .false.
f0 = 0
!if (cell%status >= CHEMOTACTIC .and. occupancy(site1(1),site1(2),site1(3))%signal /= 0) then
!	field = .true.
!	f0 = occupancy(site1(1),site1(2),site1(3))%intensity
!endif

! TESTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
field = .false.

R = par_uni(kpar)
!call random_number(R)
if (R <= dirprob(0)) then    ! case of no jump
	return
endif

! Now we must jump (if possible)

! Set up weights in the 6 principal axis directions corresponding to S1P_grad(:,:,:,:)
if (S1P_chemotaxis) then
	g = S1P_grad(:,site1(1),site1(2),site1(3))
	gamp = sqrt(dot_product(g,g))	! Magnitude of S1P gradient
	g = g/gamp
	call chemo_weights(g,wS1P)		! w(:) now holds the normalized gradient vector components
	S1Pfactor = S1P_CHEMOLEVEL*min(1.0,gamp/S1P_GRADLIM)*cell%S1P1
!	write(*,*) 'g,gamp: ',g,gamp
!	write(*,*) 'wS1P: ',wS1P
!	write(*,*) 'S1P1: ',cell%S1P1
endif
! Set up weights in the 6 principal axis directions corresponding to CXCL12_grad(:,:,:,:)
if (CXCL12_chemotaxis .and. initiated) then
	g = CXCL12_grad(:,site1(1),site1(2),site1(3))
	gamp = dot_product(g,g)	! Magnitude of CXCL12 gradient
	if (gamp < 1.0e-6) then		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! check threshold !!!!!!!!!
		CXCL12factor = 0
		wCXCL12 = 0
	else
		g = g/gamp
		call chemo_weights(g,wCXCL12)		! w(:) now holds the normalized gradient vector components
	!	CXCL12factor = CXCL12_CHEMOLEVEL*min(1.0,gamp/CXCL12_GRADLIM)*cell%RANKSIGNAL
		CXCL12factor = CXCL12_CHEMOLEVEL*min(1.0,gamp/CXCL12_GRADLIM)
		if (istep >= 5760 .and. kcell == 300) then
			write(nflog,*) 'g,gamp: ',g,gamp,CXCL12_GRADLIM
			write(nflog,*) 'wCXCL12: ',wCXCL12
			write(nflog,*) 'CXCL12factor: ',CXCL12factor
		endif
	endif
else
	CXCL12factor = 0
	wCXCL12 = 0
endif
!write(*,*) 'S1Pfactor, CXCL12factor: ',S1Pfactor, CXCL12factor
!stop
lastdir1 = cell%lastdir
p = 0
psum = 0
do irel = 1,nreldir
    p(irel) = 0
	dir1 = reldir(lastdir1,irel)
	site2 = site1 + jumpvec(:,dir1)
	if (site2(1) < 1 .or. site2(1) > NX) cycle
	if (site2(2) < NBY+2 .or. site2(2) > NY) cycle
	if (site2(3) < 1 .or. site2(3) > NZ) cycle
	if (site2(3) > NZ) then
		write(*,*) 'Error: MonoJumper: ',kcell,site1,irel,dir1,jumpvec(:,dir1),site2
		stop
	endif
	! With the call to free_site() returning region and kcell, decisions can be made
	! about transition to BLOOD, for example.
	free = free_site(site2,region,kcell2)
	if (free) then
		p(irel) = motility_factor*dirprob(irel)
		if (field) then
			! The probability of a jump is modified by the relative signal intensities of the two sites
			! This is a very crude interim treatment
!			f = occupancy(site2(1),site2(2),site2(3))%intensity
			if (f >= SIGNAL_THRESHOLD) then
				p(irel) = max(0.0,f-f0)
			else
				p(irel) = max(0.0,p(irel) + SIGNAL_AFACTOR*(f-f0))
			endif
		endif
		if (S1P_chemotaxis) then
			! the increment to the probability depends on S1Pfactor and wS1P(dir1)
			p(irel) = p(irel) + S1Pfactor*wS1P(dir1)
		endif
		if (CXCL12_chemotaxis) then
			! the increment to the probability depends on CXCL12factor and wCXCL12(dir1)
			p(irel) = p(irel) + CXCL12factor*wCXCL12(dir1)
		endif
		jmpdir(irel) = dir1
		psum = psum + p(irel)
		savesite2(:,irel) = site2
	elseif (region == BLOOD .and. vn_adjacent(dir1)) then
		if (CrossToBlood(kcell,site1)) then
			return
		endif
!		write(*,*) 'Not free_site: ',site2
	endif
enddo

if (psum == 0) then
!	cell%lastdir = random_int(1,6,kpar)
!	if (field) then
!		write(*,*) 'clustering: ',kcell
!	endif
	write(logmsg,*) 'psum = 0: ',cell%ID
	call logger(logmsg)
	return
else
    go = .true.
endif

! Now choose a direction on the basis of these probs p()
R = par_uni(kpar)
!call random_number(R)
R = psum*R
psum = 0
do irel = 1,nreldir
   	psum = psum + p(irel)
   	if (R <= psum) then
   		exit
   	endif
enddo
if (irel > nreldir) then
	irel = nreldir
endif
site2 = savesite2(:,irel)
if (irel <= nreldir) then
	dir1 = reldir(lastdir1,irel)
	if (dir1 == 0) then
		write(logmsg,*) 'dir1 = 0: ',kcell,lastdir1,irel,dir1
		call logger(logmsg)
	endif
endif
if (diagonal_jumps) then
	dir1 = fix_lastdir(dir1,kpar)
elseif (dir1 == 0) then
	dir1 = random_int(1,6,kpar)
endif

cell%site = site2
cell%lastdir = dir1
cell%lastmovestep = istep
occupancy(site2(1),site2(2),site2(3))%indx = kcell
occupancy(site1(1),site1(2),site1(3))%indx = 0
end subroutine

!---------------------------------------------------------------------
! Use for both S1P and CXCL12 chemotaxis
!---------------------------------------------------------------------
subroutine chemo_weights(g,w)
real(8) :: g(3), w(27)

w = 0
if (g(1) < 0) then		! -x
	w(5) = -g(1)		! reldir26(1,1)
else					! +x
	w(23) = g(1)		! reldir26(2,1)
endif
if (g(2) < 0) then		! -y
	w(11) = -g(2)		! reldir26(3,1)
else					! +y
	w(17) = g(2)		! reldir26(4,1)
endif
if (g(3) < 0) then		! -z
	w(13) = -g(3)		! reldir26(5,1)
else					! +z
	w(15) = g(3)		! reldir26(6,1)
endif
end subroutine



!---------------------------------------------------------------------
! For a candidate destination site, returns .true. if the site is free,
! else returns .false.  Wrapping at the boundaries is handled.
! A site may be unavailable because region /= MARROW, 
! or because kcell /= 0
!---------------------------------------------------------------------
logical function free_site(site,region,kcell)
integer :: site(3), region, kcell

if (site(1) < 1) then
	free_site = wrapx1(site,region,kcell)
	return
endif
if (site(1) > NX) then
	free_site = wrapx2(site,region,kcell)
	return
endif
if (site(2) > NY) then
	free_site = wrapy2(site,region,kcell)
	return
endif
if (site(3) < 1) then
	free_site = wrapz1(site,region,kcell)
	return
endif
if (site(3) > NZ) then
	free_site = wrapz2(site,region,kcell)
	return
endif
!if (site(2) <= NBY+1) then	!???????????  what's wrong with this?
!	kcell = 0
!	region = LAYER
!	free_site = .false.
!	return
!endif
region = occupancy(site(1),site(2),site(3))%region
if (region /= MARROW) then
	kcell = 0
	free_site = .false.
!	write(*,*) 'bad region'
	return
endif
kcell = occupancy(site(1),site(2),site(3))%indx
if (kcell == 0) then
	free_site = .true.
else
	free_site = .false.
!	write(*,*) 'occupied'
endif
end function

!---------------------------------------------------------------------
! site(1) = 0, and the location is wrapped to site(1) = NX
! The wrapped site may be an inaccessible region (non-MARROW), 
! or it may be occupied by a cell.
! Returns .true. for a vacant site, .false. otherwise
! site is replaced by the wrapped location, and the region and 
! cell index (kcell) are also returned as arguments.
!---------------------------------------------------------------------
logical function wrapx1(site,region,kcell)
integer :: site(3), region, kcell

site(1) = NX
site(2) = min(site(2),NY)
site(3) = min(site(3),NZ)
site(3) = max(site(3),1)
region = occupancy(site(1),site(2),site(3))%region
if (region /= MARROW) then
	kcell = 0
	wrapx1 = .false.
	return
endif
kcell = occupancy(site(1),site(2),site(3))%indx
if (kcell == 0) then
	wrapx1 = .true.
else
	wrapx1 = .false.
endif
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
logical function wrapx2(site,region,kcell)
integer :: site(3), region, kcell

site(1) = 1
site(2) = min(site(2),NY)
site(3) = min(site(3),NZ)
site(3) = max(site(3),1)
region = occupancy(site(1),site(2),site(3))%region
if (region /= MARROW) then
	kcell = 0
	wrapx2 = .false.
	return
endif
kcell = occupancy(site(1),site(2),site(3))%indx
if (kcell == 0) then
	wrapx2 = .true.
else
	wrapx2 = .false.
endif
end function

!---------------------------------------------------------------------
! This case is different because instead of wrapping to the opposite 
! side of the cube, the cell enters the top face (y = NY) at another
! location.  The selection is a bit arbitrary.  We choose to map 
! half of the surface, divided at x = NX/2, into the other half.
! z is preserved.
!---------------------------------------------------------------------
logical function wrapy2(site,region,kcell)
integer :: site(3), region, kcell

site(2) = NY
site(1) = min(site(1),NX)
site(1) = max(site(1),1)
site(3) = min(site(3),NZ)
site(3) = max(site(3),1)
if (site(1) > NX/2) then
	site(1) = site(1) - NX/2
else
	site(1) = site(1) + NX/2
endif	
region = occupancy(site(1),site(2),site(3))%region
if (region /= MARROW) then
	kcell = 0
	wrapy2 = .false.
	return
endif
kcell = occupancy(site(1),site(2),site(3))%indx
if (kcell == 0) then
	wrapy2 = .true.
else
	wrapy2 = .false.
endif
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
logical function wrapz1(site,region,kcell)
integer :: site(3), region, kcell

site(3) = NZ
site(2) = min(site(2),NY)
site(1) = min(site(1),NX)
site(1) = max(site(1),1)
region = occupancy(site(1),site(2),site(3))%region
if (region /= MARROW) then
	kcell = 0
	wrapz1 = .false.
	return
endif
kcell = occupancy(site(1),site(2),site(3))%indx
if (kcell == 0) then
	wrapz1 = .true.
else
	wrapz1 = .false.
endif
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
logical function wrapz2(site,region,kcell)
integer :: site(3), region, kcell

site(3) = 1
site(2) = min(site(2),NY)
site(1) = min(site(1),NX)
site(1) = max(site(1),1)
region = occupancy(site(1),site(2),site(3))%region
if (region /= MARROW) then
	kcell = 0
	wrapz2 = .false.
	return
endif
kcell = occupancy(site(1),site(2),site(3))%indx
if (kcell == 0) then
	wrapz2 = .true.
else
	wrapz2 = .false.
endif
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
logical function crossToBlood(kcell,site)
real :: R, prob
integer :: kcell, site(3), kpar=0
type(monocyte_type), pointer :: cell
real :: tnow

cell => mono(kcell)
crossToBlood = .false.
if (cell%S1P1 > S1P1_THRESHOLD) then
	prob = CROSS_PROB*(cell%S1P1 - S1P1_THRESHOLD)/(1 - S1P1_THRESHOLD)
	R = par_uni(kpar)
	if (R < prob) then
		tnow = istep*DELTA_T
		crossToBlood = .true.
		cell%status = CROSSING
		cell%exittime = tnow + CROSSING_TIME
!		write(logmsg,*) 'S1P1, prob, R: ',cell%S1P1,prob,R,cell%exittime
!		call logger(logmsg)
	endif
endif
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine chemo(S,site,p,jmpdir,pstar)
integer :: site(3), jmpdir(MAXRELDIR)
real :: S, p(0:MAXRELDIR), pstar(0:MAXRELDIR)
integer :: u(3), i
real :: g(3), gdotu, psum

g = S1P_grad(:,site(1),site(2),site(3))
do i = 1,MAXRELDIR
	if (p(i) == 0) cycle
	u = jumpvec(:,jmpdir(i))
	gdotu = g(1)*u(1) + g(2)*u(2) + g(3)*u(3)
	pstar(i) = p(i) + Kchemo*S*gdotu
	pstar(i) = max(pstar(i),0.0)
	pstar(i) = min(pstar(i),1.0)
	psum = psum + pstar(i)
enddo
if (psum > 1) then
	pstar(1:MAXRELDIR) = pstar(1:MAXRELDIR)/psum
	pstar(0) = 0
else
	pstar(0) = 1 - psum
endif
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
real function clast_size(count)
integer :: count
clast_size = CLAST_RADIUS_FACTOR*sqrt(real(count))
end function

!---------------------------------------------------------------------
! Choose one of the 6 axis directions to allocate to the direction of
! previous step given by jump, which takes values from 0 - 27.  
! This could be already on an axis, or could be a D2 or D3 diagonal, 
! in which case the choice is random.
!---------------------------------------------------------------------
integer function fix_lastdir(jump,kpar)
integer :: jump, kpar
integer :: k,nax
!                            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
integer :: naxes(0:27) = (/  0, 3, 2, 3, 2, 1, 2, 3, 2, 3, 2, 1, 2, 1, 0, 1, 2, 1, 2, 3, 2, 3, 2, 1, 2, 3, 2, 3 /)
integer :: axes(3,27) = reshape( (/ &
!   1      2      3      4      5      6      7      8      9     10     11     12     13     14
1,3,5, 1,3,0, 1,3,6, 1,5,0, 1,0,0, 1,6,0, 1,4,5, 1,4,0, 1,4,6, 3,5,0, 3,0,0, 3,6,0, 5,0,0, 0,0,0, &
!  15     16     17     18     19     20     21     22     23     24     25     26     27 
6,0,0, 4,5,0, 4,0,0, 4,6,0, 2,3,5, 2,3,0, 2,3,6, 2,5,0, 2,0,0, 2,6,0, 2,4,5, 2,4,0, 2,4,6 /), (/3,27/) )

if (diagonal_jumps) then
	nax = naxes(jump)
	if (nax == 0) then
	    write(*,*) 'Should not get here: fix_lastdir: nax=0: jump: ',jump
	    stop
		fix_lastdir = random_int(1,6,kpar)
	elseif (nax == 1) then
		fix_lastdir = axes(1,jump)
	else
		k = random_int(1,nax,kpar)
		fix_lastdir = axes(k,jump)
	endif
else
	write(*,*) 'fix_lastdir: Not MOORE model: ',jump
	stop
endif
end function

!---------------------------------------------------------------------
! Need to precompute array reldir(:,:)
! FOR NRELDIR = 6 (diagonal_jumps = .false.)
! The directions 1 - 6 are defined according to the neumann array,
! i.e. 1 is -x, 2 is +x, 3 is -y, 4 is +y, 5 is -z, 6 is +z
! The relative directions are numbered as follows:
! irel = 1 gives the same direction as lastdir
! irel = 6 gives the opposite direction to lastdir
! irel = 2,3,4,5 cover the remaining directions - order not important
! at the moment since all directions normal to lastdir are equally likely
! FOR NRELDIR = 17 (diagonal_jumps = .true.)
! The reldir(:,:) array uses the same set of 6 previous jump directions,
! restricted to the axes.  Diagonal previous jumps are accommodated by
! making a random selection of an axis direction from either 2 or 3
! possibilities (depending on whether it was a D2 or D3 jump).
! The possible jumps, now including diagonal moves, are stored in jumpvec(:)
! The jumpvec(:) entries are ordered with z varying fastest, then y, then x,
! each ranging -1, 0, +1.
! For dirprob(:) calculation:
! There are 3 groups of jumps, with sets of probs
! Group 1: Q(1) = P(1), D1 jump in the same as last direction
! Group 2: Q(2) = P(2)+P(3)+P(4)+P(5) (D2 jumps) + P(6)+P(7)+P(8)+P(9) (D3 jumps)
! Group 3: Q(3) = P(10)+P(11)+P(12)+P(13) (D1 jumps) + P(14)+P(15)+P(16)+P(17) (D2 jumps)
! i.e. Q(1) = prob of no direction change, Q(2) = prob of direction change < 90 deg,
! Q(3) = prob of direction change = 90 deg.
! Setting D1 jumps (von Neumann) to 1, the D2 jumps have length L2 = sqrt(2)
! and the D3 jumps have length L3 = sqrt(3)
! Therefore jump probs must be scaled appropriately to give symmetry within a set
! In other words:
!		P(i)/P(j) = L2/L3 for j=6,7,8,9 j=2,3,4,5
! and	P(i)/P(j) = 1/L2  for i=14,15,16,17 j=10,11,12,13
! Then choosing P(2) to represent D2 jump prob in Group 2, and P(10) to
! represent D1 jump prob in Group 3:
!		Q(2) = 4(1 + L2/L3).P(2)
!		Q(3) = 4(1 + 1/L2).P(10)
! It is always the case that Q(1) + Q(2) + Q(3) = BETA
!
! When RHO = 1, Q(1) = BETA = p1, Q(2) + Q(3) = 0
!
! When RHO = 0, P(10) = P(1), and P(2) = P(1)/L2
! therefore:	Q(3) = 4(1 + 1/L2).P(1)
! and			Q(2) = 4(1 + L2/L3).P(1)/L2
! giving:		P(1).[1 + 4(1 + L2/L3)/L2 + 4(1 + 1/L2)] = BETA
!				P(1) = BETA/[1 + 4(1 + L2/L3)/L2 + 4(1 + 1/L2)] = p0
! and:			Q(3)/Q(2) = (1 + L2)/(1 + L2/L3) = q320
!
! Now for 0 < RHO < 1
! first set Q(1) = P(1) to vary linearly with RHO from 
! the value at 0 (p0) to the value at 1 (p1)
!	P(1) = p0 + RHO*(p1-p0)
! Now we know that Q(2) + Q(3) = BETA - Q(1)
! and the third parameter ALPHA is used to determine how Q(3)/Q(2)
! varies with RHO, by making it change linearly from q23 at RHO=0
! to ALPHA*q23 at RHO = 1
! Now we have all the info to solve for Q(1), Q(2), Q(3) and all P(:)
!
! Revised Moore model.
! Now allow reverse directions, but disallow D3 diagonal jumps, giving
! 18 possible directions.
! In the preferred direction the surface is a prolate spheroid, with
! a = 1, b = 1-RHO.
! In the reverse direction the surface is either a sphere (a = b) or
! an oblate spheroid (a = b^2).
!---------------------------------------------------------------------
subroutine make_reldir
integer :: lastdir,irel,k,ix,iy,iz,i,site(3)
integer :: reldir18(6,18),reldir26(6,26)
real :: Q(3),L2,L3,p0,p1,q320,q321,q32,qsum,a,b,c,d(MAXRELDIR),R(MAXRELDIR),E, theta, p(8), psum, u2, u(3)

if (MODEL == NEUMANN_MODEL) then
	diagonal_jumps = .false.
	njumpdirs = 6
	nreldir = 6
	reldir = 0
	do lastdir = 1,6
		reldir(lastdir,1) = lastdir
		if (mod(lastdir,2) == 0) then
			reldir(lastdir,6) = lastdir - 1
		else
			reldir(lastdir,6) = lastdir + 1
		endif
		irel = 2
		do k = 1,6
			if (k /= reldir(lastdir,1) .and. k /= reldir(lastdir,6)) then
				reldir(lastdir,irel) = k
				irel = irel + 1
			endif
		enddo		
	enddo
!   	write(*,*) 'reldir'
!    do lastdir = 1,6
!	    write(*,'(6i6)') (reldir(lastdir,k),k=1,6)
!    enddo
	jumpvec(:,1:6) = neumann(:,1:6)

else
	if (MODEL == MOORE18_MODEL) then
    	nreldir = 18
	elseif (MODEL == MOORE26_MODEL) then
    	nreldir = 26
    endif
    njumpdirs = 27
	diagonal_jumps = .true.
	k = 0
	do ix = -1,1
		do iy = -1,1
			do iz = -1,1
				k = k+1
				jumpvec(:,k) = (/ix,iy,iz/)
				u = jumpvec(:,k)
				u2 = dot_product(u,u)
				if (u2 > 0) then
					unitjump(:,k) = u/sqrt(u2)
				else
					unitjump(:,k) = 0
				endif
			enddo
		enddo
	enddo

! Data for revised Moore18 model.  D3 jumps are excluded, but reverse directions are allowed.
	reldir18(1,:) = (/  5,  2, 4, 6, 8, 11,13,15,17, 10,12,16,18, 22,24,26,20, 23 /)	! -x
	reldir18(2,:) = (/ 23, 20,22,24,26, 11,13,15,17, 10,12,16,18,  2, 4, 6, 8,  5 /)	! +x
	reldir18(3,:) = (/ 11,  2,10,12,20,  5,13,15,23,  4, 6,22,24,  8,16,18,26, 17 /)	! -y
	reldir18(4,:) = (/ 17,  8,16,18,26,  5,13,15,23,  4, 6,22,24,  2,10,12,20, 11 /)	! +y
	reldir18(5,:) = (/ 13,  4,10,16,22,  5,11,17,23,  2, 8,20,26,  6,12,18,24, 15 /)	! -z
	reldir18(6,:) = (/ 15,  6,12,18,24,  5,11,17,23,  2, 8,20,26,  4,10,16,22, 13 /)	! +z

! Data for revised Moore26 model.  D3 jumps are included, but reverse directions are allowed.
	reldir26(1,:) = (/  5,  2, 4, 6, 8,  1, 3, 7, 9, 11,13,15,17, 10,12,16,18, 19,21,27,25, 20,22,24,26, 23 /)	! -x
	reldir26(2,:) = (/ 23, 20,22,24,26, 19,21,27,25, 11,13,15,17, 10,12,16,18,  1, 3, 7, 9,  2, 4, 6, 8,  5 /)	! +x
	reldir26(3,:) = (/ 11,  2,10,12,20,  1, 3,19,21,  5,13,15,23,  4, 6,22,24,  7, 9,25,27,  8,16,18,26, 17 /)	! -y
	reldir26(4,:) = (/ 17,  8,16,18,26,  7, 9,25,27,  5,13,15,23,  4, 6,22,24,  1, 3,19,21,  2,10,12,20, 11 /)	! +y
	reldir26(5,:) = (/ 13,  4,10,16,22,  1, 7,19,25,  5,11,17,23,  2, 8,20,26,  3, 9,21,27,  6,12,18,24, 15 /)	! -z
	reldir26(6,:) = (/ 15,  6,12,18,24,  3, 9,21,27,  5,11,17,23,  2, 8,20,26,  1, 7,19,25,  4,10,16,22, 13 /)	! +z

    if (MODEL == MOORE18_MODEL) then
        reldir(1:6,1:18) = reldir18
    elseif (MODEL == MOORE26_MODEL) then
        reldir = reldir26
    endif

endif
call compute_dirprobs
!write(*,*)  'dirprob: '
!write(*,'(7f6.3)') dirprob(0:nreldir)
qsum = 0
do i = 0,nreldir
    qsum = qsum + dirprob(i)
enddo
call make_vn_adjacent
!write(*,*) 'sum: ',qsum
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine make_vn_adjacent
integer :: dir, v(3)

do dir = 1,njumpdirs
	v = jumpvec(:,dir)
	if ((v(1) == 0 .and. v(2) == 0 .and. v(3) /= 0) .or. &
	    (v(1) == 0 .and. v(2) /= 0 .and. v(3) == 0) .or. &
	    (v(1) /= 0 .and. v(2) == 0 .and. v(3) == 0)) then
	   vn_adjacent(dir) = .true.
	else
		vn_adjacent(dir) = .false.
	endif
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine compute_dirprobs
integer :: i
real :: L2,L3,b,c,e,d(MAXRELDIR),R(MAXRELDIR)

if (MODEL == NEUMANN_MODEL) then
!	diagonal_jumps = .false.
	nreldir = 6
	! New Neumann model (with reverses)
	b = 1-RHO
	b = b*b
	e = BETA/(1 + 4*b + b**2)
	dirprob(0) = 1 - BETA	! this is the prob of no jump
	dirprob(1) = e
	dirprob(2:5) = e*b
	dirprob(6) = e*b**2
else
    ! New prolate spheroid approach
!	diagonal_jumps = .true.
    if (MODEL == MOORE18_MODEL) then
	    nreldir = 18
	    L2 = sqrt(2.0)
	    b = 1 - RHO
	    b = b*b		! revised version
	    c = b*b		! oblate spheroid on reverse side
	    d(1) = 1
	    R(1) = 1
	    d(2:5) = L2
	    R(2:5) = L2*b/sqrt(b*b+1)
	    d(6:9) = 1
	    R(6:9) = b
	    d(10:13) = L2
	    R(10:13) = b
	    d(14:17) = L2
	    R(14:17) = L2*b*c/sqrt(b*b+c*c)
	    d(18) = 1
	    R(18) = c
	    e = BETA/(1 + 4*b/sqrt(b*b+1) + 4*b + 4*b/L2 + 4*b*c/sqrt(b*b+c*c) + c)
	    dirprob(0) = 1 - BETA
	    dirprob(1:nreldir) = e*R(1:nreldir)/d(1:nreldir)
    elseif (MODEL == MOORE26_MODEL) then
    	nreldir = 26
	    L2 = sqrt(2.0)
	    L3 = sqrt(3.0)
	    b = 1 - RHO
	    b = b*b		! revised version
	    c = b*b		! oblate spheroid on reverse side
	    d(1) = 1
	    R(1) = 1
	    d(2:5) = L2
	    R(2:5) = L2*b/sqrt(b*b+1)
	    d(6:9) = L3
	    R(6:9) = L2*b/sqrt(b*b+1)
	    d(10:13) = 1
	    R(10:13) = b
	    d(14:17) = L2
	    R(14:17) = b
	    d(18:21) = L3
	    R(18:21) = L2*b*c/sqrt(b*b+c*c)
	    d(22:25) = L2
	    R(22:25) = L2*b*c/sqrt(b*b+c*c)
	    d(26) = 1
	    R(26) = c
	    e = BETA/(1 + 4*(1+L2/L3)*b/sqrt(b*b+1) + 4*b + 4*b/L2 + 4*(1+L2/L3)*b*c/sqrt(b*b+c*c) + c)
	    dirprob(0) = 1 - BETA
	    dirprob(1:nreldir) = e*R(1:nreldir)/d(1:nreldir)
    endif
endif
end subroutine

end module

