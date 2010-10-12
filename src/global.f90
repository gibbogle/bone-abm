module global
use defines
use par_zig_mod
use winsock

implicit none

type(occupancy_type), allocatable, target :: occupancy(:,:,:)
type(monocyte_type), allocatable, target :: mono(:)
type(osteoclast_type), allocatable, target :: clast(:)
type(osteoblast_type), allocatable :: blast(:)
type(stem_type), allocatable :: stem(:)

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
!DEC$ ATTRIBUTES DLLEXPORT :: use_TCP

! Parameters read from inputfile
real :: BETA							! speed: 0 < beta < 1
real :: RHO								! persistence: 0 < rho < 1
integer :: NX,NY,NZ						! size of cubical region
integer :: exit_rule					! 1 = no chemotaxis, 2 = chemotaxis
integer :: exit_region					! region for cell exits 1 = capillary, 2 = sinusoid
real :: cross_prob						! probability (/timestep) of monocyte egress to capillary
real :: days							! number of days to simulate
integer :: seed(2)						! seed vector(1) for the RNGs
integer :: ncpu							! number of threads, not used currently
integer :: NT_GUI_OUT					! interval between GUI outputs (timesteps)
integer :: SPECIES						! animal species (0=mouse, 1=human)

! Misc parameters
real :: DELTA_X
integer :: Nsteps
integer :: NMONO_INITIAL
integer :: nmono, mono_cnt, nsignal, nclast, nleft

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine logger(msg)
character*(*) :: msg
integer :: error
logical :: isopen

error = 0
if (use_TCP) then
    if (awp_0%is_open) then
        call winsock_send(awp_0,trim(msg),len_trim(msg),error)
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
