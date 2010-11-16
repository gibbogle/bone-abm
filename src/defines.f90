module defines
implicit none

integer, parameter :: nfinp=10, nflog=11, nfpos=12
integer, parameter :: TCP_PORT_0 = 5000		! main communication port (logging)
integer, parameter :: TCP_PORT_1 = 5001		! data transfer port (plotting)

integer, parameter :: MARROW = 0
integer, parameter :: PIT = -1
integer, parameter :: BLOOD = -2
integer, parameter :: BONE = -3

integer, parameter :: MONOCYTE = 10
integer, parameter :: OSTEOBLAST = 11
integer, parameter :: STROMAL = 12
integer, parameter :: OSTEOCLAST = 13
integer, parameter :: STEMCELL = 14
integer, parameter :: OFF = 0
integer, parameter :: ON = 1

integer, parameter :: DEAD = -1
integer, parameter :: MOTILE = 1
integer, parameter :: ALIVE = 1
integer, parameter :: CHEMOTACTIC = 2
integer, parameter :: STICKY = 3
integer, parameter :: FUSING = 4
integer, parameter :: FUSED = 5
integer, parameter :: CLOSED = 6
integer, parameter :: MATURE = 6
integer, parameter :: CROSSING = 7
integer, parameter :: LEFT = 8

integer, parameter :: NEUMANN_MODEL = 1
integer, parameter :: MOORE18_MODEL = 2
integer, parameter :: MOORE26_MODEL = 3
integer, parameter :: MAXRELDIR = 26
integer, parameter :: neumann(3,6) = reshape((/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /), (/3,6/))


integer, parameter :: NBY = 8
integer, parameter :: MAX_MONO = 50000
integer, parameter :: MAX_CLAST = 100
integer, parameter :: MAX_CAP = 100
integer, parameter :: MAX_SIGNAL = MAX_CLAST
integer, parameter :: MAX_NCLUMP = 20
integer, parameter :: MAX_CLUMP_CELLS = 40
real, parameter :: DELTA_T = 0.25		! minutes
real, parameter :: BIGTIME = 1.0e10
logical, parameter :: FAST_DISPLAY = .false.

! GUI parameters
character*(12), parameter :: stopfile = 'stop_dll'
character*(13), parameter :: pausefile = 'pause_dll'

!real, parameter :: capR = 1.5
!real, parameter :: MONOCYTE_DIAMETER = 10	! um
!real, parameter :: MONO_PER_MM3 = 2000
!integer, parameter :: NSTEM = 20
!real, parameter :: STEM_CYCLETIME = 6*60	! 6 hours

!real, parameter :: SIGNAL_RADIUS = 10	! radius of influence (in lattice sites) of bone signal
!real, parameter :: SIGNAL_THRESHOLD = 0.14
!real, parameter :: SIGNAL_AFACTOR = 0.4	! field amplification factor
!integer, parameter :: MTHRESHOLD = 25

! Osteoclast parameters
!real, parameter :: MAX_RESORPTION_RATE = 0.002
!real, parameter :: MAX_RESORPTION_D = 10
!integer, parameter :: MAX_RESORPTION_N = 30
!real, parameter :: CLAST_LIFETIME = 96*60
!real, parameter :: CROSSING_TIME = 2*60
!real, parameter :: FUSING_TIME = 2*60
!real, parameter :: CLAST_DWELL_TIME0 = 4*60
!real, parameter :: CLAST_DWELL_TIME = 3*60


! S1P1 parameters
logical, parameter :: S1P_chemotaxis = .true.
real, parameter :: S1P1_THRESHOLD = 0.5
real, parameter :: S1P1_BASERATE = 1.0/(6*60)	! 6 hours

! RANKL parameters
logical, parameter :: use_RANK = .true.
real, parameter :: RANKSIGNAL_rateconstant = 4.0
real, parameter :: RANKSIGNAL_halflife = 12		! hours
real, parameter :: ST1 = 0.3	! -> CHEMOTACTIC
real, parameter :: ST2 = 0.5	! -> STICKY

! Clump parameters
integer, parameter :: CLUMP_THRESHOLD = 25

type pit_type
!	logical :: active
	integer :: site(3)
	real :: rate
!	real :: fraction
end type

type monocyte_type
    integer :: ID
    integer(2) :: site(3)
	integer(1) :: region
	integer(1) :: iclast
    integer(1) :: status
	integer(1) :: lastdir
	integer(1) :: iclump
	real :: S1P1				! level of S1P1 expression
	real :: RANKSIGNAL			! integrated RANK signal
	real :: stickiness
    real :: entrytime			! time that the cell entered the marrow (from blood or stem cell division)
    real :: exittime			! time that the cell left the marrow (for the blood, or to form an osteoclast)
    real :: dietime             ! time cell will die
!    type(clump_type), pointer :: clump
!    type(O_type),    pointer :: optr    ! because NULL is used by winsock (from ifwinty).  NULLIFY() instead.
end type

type osteoblast_type
    integer :: ID
    integer :: site(3)
    integer :: step
    integer(2) :: status
	integer(2) :: lastdir
    real :: entrytime			! time that the cell entered the marrow (from blood or stem cell division)
    real :: dietime             ! time cell will die
end type

type osteoclast_type
    integer :: ID
    integer :: site(3)
	real :: normal(3)
    integer(2) :: status
    real :: fusetime			! time that the monocytes began to fuse to form the osteoclast
    real :: entrytime			! time that the cell became a mature osteoclast
    real :: movetime			! time that the cell will move on
    real :: dietime             ! time cell will die
	integer :: count
	integer :: mono(100)
	integer :: npit
	type(pit_type), allocatable :: pit(:)
end type

type stromal_type
    integer :: ID
    integer :: site(3)
    integer :: step
    integer(2) :: status
	integer(2) :: lastdir
    real :: entrytime			! time that the cell entered the marrow (from blood or stem cell division)
    real :: dietime             ! time cell will die
end type

type stem_type
    integer :: ID
    integer :: site(3)
!    integer :: step
    integer(2) :: status
!	integer(2) :: lastdir
!	real :: entrytime			! time that the cell entered the marrow (from blood or stem cell division)
	real :: dividetime          ! time cell will divide
end type

type occupancy_type
	integer :: region
	integer :: species
    integer :: indx
    integer :: signal
    real :: intensity
    real :: bone_fraction
end type

type signal_type
	logical :: active
	integer :: site(3)
	real :: normal(3)
end type

type capillary_type
	real :: pos1(3)
	real :: pos2(3)
	real :: radius
	real :: length
	real :: surface_area
end type

type clump_type
	integer :: ncells
	integer :: list(MAX_CLUMP_CELLS)
	integer :: status
	real :: starttime
end type

end module