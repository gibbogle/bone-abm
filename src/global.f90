module global
use defines
use par_zig_mod
use winsock

implicit none
save

type(occupancy_type), allocatable, target :: occupancy(:,:,:)
type(monocyte_type), allocatable, target :: mono(:)
type(osteoclast_type), allocatable, target :: clast(:)
type(osteoblast_type), allocatable :: blast(:)
type(stem_type), allocatable :: stem(:)

type(capillary_type), allocatable :: capillary(:)
type(signal_type) :: signal(MAX_SIGNAL)

logical :: diagonal_jumps, clear_to_send, simulation_start, stopped
integer :: Mnodes = 1

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
real :: CLAST_DWELL_TIME = 3*60			! mins
real :: MAX_RESORPTION_RATE = 0.02		! um/min
real :: MAX_RESORPTION_D = 10			! um
integer :: MAX_RESORPTION_N = 30

! Signal parameters (NOT USED)
real :: SIGNAL_RADIUS					! radius of influence of bone signal (um -> grids) (10)
real :: SIGNAL_THRESHOLD				! defines the high-signal region, near the source (0.14)
real :: SIGNAL_AFACTOR					! field amplification factor (0.4)
integer :: MTHRESHOLD					! number of monocytes in the high-signal region that triggers fusing (25)

integer :: in_per_hour					! rate of influx of OP monocytes from the blood (cells/hour)
integer :: exit_rule					! 1 = no chemotaxis, 2 = chemotaxis
integer :: exit_region					! region for cell exits 1 = capillary, 2 = sinusoid
real :: cross_prob						! probability (/timestep) of monocyte egress to capillary
real :: days							! number of days to simulate
integer :: seed(2)						! seed vector(1) for the RNGs
integer :: ncpu							! number of threads, not used currently
integer :: NT_GUI_OUT					! interval between GUI outputs (timesteps)
integer :: SPECIES						! animal species (0=mouse, 1=human)

! Misc parameters
real :: DELTA_X, PI
integer :: Nsteps
integer :: NX,NY,NZ						! size of region
integer :: NMONO_INITIAL, NSTEM
integer :: nmono, mono_cnt, nsignal, nclast, nborn, nleft, ncap, nentrysites, nclump
integer, allocatable :: entrysite(:,:)
type(clump_type), target :: clump(MAX_NCLUMP)
real :: RANKSIGNAL_decayrate			! from RANKSIGNAL_halflife

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

end module
