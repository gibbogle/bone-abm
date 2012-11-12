module global
use omp_lib
use defines
use par_zig_mod
use winsock

implicit none
save

type(occupancy_type), allocatable, target :: occupancy(:,:,:)
type(surface_type), allocatable, target :: surface(:,:)
type(monocyte_type), allocatable, target :: mono(:)
type(osteoclast_type), allocatable, target :: clast(:)
type(osteoblast_type), allocatable, target :: blast(:)
type(stem_type), allocatable :: stem(:)

type(capillary_type), allocatable :: capillary(:)
type(signal_type) :: signal(MAX_SIGNAL)

logical :: diagonal_jumps, clear_to_send, simulation_start, stopped

integer :: istep

character*(128) :: inputfile
character*(128) :: outputfile
character*(128) :: resultfile
character*(128) :: runningfile
character*(256) :: logmsg
TYPE(winsockport) :: awp_0, awp_1, awp_2, awp_3
logical :: use_TCP = .true.
logical :: use_CPORT1 = .true.
!DEC$ ATTRIBUTES DLLEXPORT :: use_TCP

! Parameters read from inputfile
real :: MONOCYTE_DIAMETER = 10			! um
real :: BETA							! speed: 0 < beta < 1 
real :: RHO								! persistence: 0 < rho < 1
real :: S1P_CHEMOLEVEL
real :: S1P_KDIFFUSION
real :: S1P_KDECAY
real :: S1P_GRADLIM
real :: S1P1_THRESHOLD
real :: S1P1_BASERATE

real :: X_SIZE
real :: Y_SIZE
real :: CAPILLARY_DIAMETER = 3
integer :: MONO_PER_MM3 = 2000
integer :: STEM_PER_MM3
real :: STEM_CYCLETIME = 6*60	! 6 hours
real :: CROSSING_TIME = 2*60

! Osteoclast parameters
real :: FUSING_TIME = 2*60				! mins
real :: CLAST_LIFETIME = 96*60			! days -> mins
!real :: CLAST_DWELL_TIME0 = 4*60		! mins
real :: CLAST_DWELL_TIME = 60			! mins
real :: MAX_RESORPTION_RATE = 0.2		! um/min
real :: MAX_RESORPTION_D = 10			! um
integer :: MAX_RESORPTION_N = 10

! Signal parameters (NOT USED)
real :: SIGNAL_RADIUS					! radius of influence of bone signal (um -> grids) (10)
real :: SIGNAL_THRESHOLD				! defines the high-signal region, near the source (0.14)
real :: SIGNAL_AFACTOR					! field amplification factor (0.4)
integer :: MTHRESHOLD					! number of monocytes in the high-signal region that triggers fusing (25)

! Motion variables
integer :: MODEL = MOORE26_MODEL
integer :: reldir(6,26)
integer :: njumpdirs, nreldir
integer :: jumpvec(3,MAXRELDIR+1)
real :: unitjump(3,MAXRELDIR+1)
real :: dirprob(0:MAXRELDIR)
logical :: vn_adjacent(MAXRELDIR+1)
integer :: dir2D(3,8) = reshape((/ -1,0,-1, -1,0,0, -1,0,1, 0,0,1, 1,0,1, 1,0,0, 1,0,-1, 0,0,-1/),(/3,8/))
integer :: clumpoffset(3,MAX_NCLUMP)
logical :: CXCL12_chemotaxis = .false.

integer :: in_per_hour					! rate of influx of OP monocytes from the blood (cells/hour)
integer :: exit_rule					! 1 = no chemotaxis, 2 = chemotaxis
integer :: exit_region					! region for cell exits 1 = capillary, 2 = sinusoid
real :: cross_prob						! probability (/timestep) of monocyte egress to capillary

real :: days							! number of days to simulate
integer :: seed(2)						! seed vector(1) for the RNGs
integer :: ncpu							! number of threads, not used currently
integer :: Mnodes
integer :: NT_GUI_OUT					! interval between GUI outputs (timesteps)
integer :: SPECIES						! animal species (0=mouse, 1=human)

! Misc parameters
real :: DELTA_X, PI, DELTA_T_OC
integer :: Nsteps
integer :: NX,NY,NZ						! size of region
type(patch_type) :: patch
integer :: NMONO_INITIAL, NSTEM, NBLAST_INITIAL
integer :: nmono, mono_cnt, nsignal, nclast, nblast, nborn, nleft, ncap, nentrysites, nclump, nliveclast
integer, allocatable :: entrysite(:,:)
type(clump_type), target :: clump(MAX_NCLUMP)
real :: RANKSIGNAL_decayrate			! from RANKSIGNAL_halflife
real :: CXCL12_KDECAY					! from CXCL12_HALFLIFE
real :: CXCL12_GRADLIM = 5.0e-4			! was equal to the initial max CXCL12 gradient at NBY+3
										! Now need to guess the value to use!!!!!!!!
logical :: stuck
logical :: initiated
logical :: use_capillary

logical, parameter :: TESTING_OC = .true.

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine logger(msg)
character*(*) :: msg
integer :: error
logical :: isopen
character*(1) :: LF = char(94)

error = 0
if (use_TCP) then
    if (awp_0%is_open) then
        call winsock_send(awp_0,trim(msg)//LF,len_trim(msg)+1,error)
    endif
else
	write(*,*) trim(msg)
endif
inquire(unit=nflog,OPENED=isopen)
if (isopen) then
	write(nflog,*) 'msg: ',trim(msg)
	if (error /= 0) then
	    write(nflog,'(a,i4)') 'winsock_send error: ',error
	    close(nflog)
	endif
endif
!if (error /= 0) stop
end subroutine

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
!use IFPORT
integer :: n1,n2,kpar
integer :: k,R

R = par_shr3(kpar)
k = abs(R)
random_int = n1 + mod(k,(n2-n1+1))

end function

!---------------------------------------------------------------------
! Randomly selects from n possibilities with probabilities p(:)
!---------------------------------------------------------------------
integer function random_selection(p,n)
real(8) :: p(*)
integer :: n
integer :: i, kpar = 0
real(8) :: psum, R

R = par_uni(kpar)
psum = 0
do i = 1,n
	psum = psum + p(i)
	if (psum >= R) exit
enddo
random_selection = min(i,n)
end function

!---------------------------------------------------------------------
! Magnitude of real vector v(:)
!---------------------------------------------------------------------
real function rnorm(v)
real :: v(3)
rnorm = sqrt(dot_product(v,v))
end function

!---------------------------------------------------------------------
! Magnitude of integer vector v(:)
!---------------------------------------------------------------------
real function inorm(v)
integer :: v(3)
inorm = sqrt(real(dot_product(v,v)))
end function

!------------------------------------------------------------------------------------------------
! Changed to make signal proportional to fraction of target depth left to remove.
! This still is not satisfactory!!!
! Signal is now proportional to the remaining bone depth to be removed, normalized by MAX_PIT_DEPTH
! When the attractiveness of a site to an OC is being computed, the surface seal is taken into
! account - the effective signal is signal*(1-seal).
!------------------------------------------------------------------------------------------------
real function GetSignal(x,z)
integer :: x, z

!GetSignal = (surface(x,z)%target_depth - surface(x,z)%depth)/MAX_PIT_DEPTH
if (surface(x,z)%target_depth > 0) then
!	GetSignal = 2*(surface(x,z)%target_depth - surface(x,z)%depth)/(surface(x,z)%target_depth + MAX_PIT_DEPTH)
	GetSignal = (surface(x,z)%target_depth - surface(x,z)%depth)/(MAX_PIT_DEPTH)
else
	GetSignal = 0
endif
!write(*,*) 'GetSignal: ',x,z,surface(x,z)%target_depth,surface(x,z)%depth
end function

!------------------------------------------------------------------------------------------------
! The bone resorption rate at a given pit site (x,z) depends on:
!	Nm = the number of monocytes that fused to make the osteoclast
!   Np = number of pit sites that the osteoclast covers
!	(d = the distance of the target bone site from the osteoclast centre
!   The depth factor df decreases linearly to zero as d goes from 0 to MAX_RESORPTION_D)
! The volume rate of resorption per pit site is:
!   (MAX_RESORPTION_RATE/MAX_RESORPTION)*(Nm/Np) um^3/min
! which is converted to grids/min, /(DELTA_X^3)
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
resorptionRate = (MAX_RESORPTION_RATE*Nm)/(MAX_RESORPTION_N*Np*DELTA_X**3)
!write(*,*) 'resorptionRate: ',Nm,Np,resorptionRate
!write(*,*) MAX_RESORPTION_RATE,Nm,MAX_RESORPTION_N,Np,24*60,DELTA_X 
end function


!------------------------------------------------------------------------------------------------
! Compute the signal that an OB is receiving from neighbouring bone sites.
! The OB should be kept out of the pit.  How?  It should always keep its distance from all OCs,
! and avoid sites with depth > 0, while trying to be near high signal.
! Try making blast signal on/off.
!------------------------------------------------------------------------------------------------
real function BlastSignal(pblast)
type(osteoblast_type), pointer :: pblast
integer :: bsite(3), dx, dz, x, z
real :: d2, r2, sum, factor

!if (pblast%status /= ALIVE) then
!	BlastSignal = 0
!else
!	BlastSignal = 30
!endif
!return

BlastSignal = 0
if (pblast%status /= ALIVE) return
bsite = pblast%site
if (surface(bsite(1),bsite(3))%seal == 1) return
r2 = OB_SIGNAL_RADIUS**2
factor = 1.0
sum = 0
do dx = -OB_SIGNAL_RADIUS,OB_SIGNAL_RADIUS
	x = bsite(1) + dx
	do dz = -OB_SIGNAL_RADIUS,OB_SIGNAL_RADIUS
		z = bsite(3) + dz
		d2 = dx**2 + dz**2
		if (d2 > r2) cycle
		if (surface(x,z)%iclast > 0) then
			factor = 0.5
		endif
!		sum = sum + surface(x,z)%signal*(1-surface(x,z)%seal)
		sum = sum + factor*surface(x,z)%signal
	enddo
enddo
BlastSignal = OB_SIGNAL_FACTOR*sum
end function

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
real function OCradius(n)
integer :: n

OCradius = (25/DELTA_X)*sqrt(n/20.)
end function

end module
