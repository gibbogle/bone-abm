module defines
implicit none
save

integer, parameter :: nfinp=10, nflog=11, nfpos=12
integer, parameter :: TCP_PORT_0 = 5000		! main communication port (logging)
integer, parameter :: TCP_PORT_1 = 5001		! data transfer port (plotting)

integer, parameter :: MARROW = 0
integer, parameter :: PIT = -1
integer, parameter :: BLOOD = -2
integer, parameter :: BONE = -3
integer, parameter :: LAYER = -4

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
integer, parameter :: CLUMPED = 4
integer, parameter :: FUSING = 5
integer, parameter :: FUSED = 6
integer, parameter :: OSTEO = 7
integer, parameter :: CROSSING = 8
integer, parameter :: LEFT = 9

integer, parameter :: RESORBING = 2
integer, parameter :: MOVING = 3

integer, parameter :: NEUMANN_MODEL = 1
integer, parameter :: MOORE18_MODEL = 2
integer, parameter :: MOORE26_MODEL = 3
integer, parameter :: MAXRELDIR = 26
integer, parameter :: neumann(3,6) = reshape((/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /), (/3,6/))


integer, parameter :: NBY = 8
integer, parameter :: MAX_MONO = 50000
integer, parameter :: MAX_CLAST = 100
integer, parameter :: MAX_BLAST = 200
integer, parameter :: MAX_CAP = 100
integer, parameter :: MAX_SIGNAL = 10000
integer, parameter :: MAX_NCLUMP = 50
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

!real, parameter :: S1P_CHEMOLEVEL = 0.1		! 0 -> 1
!real, parameter :: S1P_GRADLIM = 0.02
!real, parameter :: S1P1_THRESHOLD = 0.5
!real, parameter :: S1P1_BASERATE = 1.0/(6*60)	! 6 hours 

! RANKL parameters
logical, parameter :: use_RANK = .true.
real, parameter :: RANKL_KDIFFUSION = 2		! 1/10 of approx 1000 um^2/min
! http://www.math.ubc.ca/~ais/website/status/diffuse.html
!real, parameter :: RANKL_KDECAY = 0.00001	! <=== compute from RANKL_HALFLIFE
real, parameter :: RANKL_HALFLIFE = 12*60
real, parameter :: RANKSIGNAL_rateconstant = 1.0
real, parameter :: RANKSIGNAL_halflife = 6*60		! mins
real, parameter :: RANK_BONE_RATIO = 0.3			! ratio of RANKL secretion to bone signal strength.
real, parameter :: ST1 = 0.3	! -> CHEMOTACTIC
real, parameter :: ST2 = 0.5	! -> STICKY

logical, parameter :: RANKL_chemotaxis = .true.
real, parameter :: RANKL_CHEMOLEVEL = 0.3
!real, parameter :: RANKL_GRADLIM = 0.0005

! Clump parameters
integer, parameter :: CLUMP_THRESHOLD = 25
real, parameter :: CLUMP_SEPARATION = 6
real, parameter :: CLUMP_FALL_PROB = 0.001	! arbitrary

! Osteoclast parameters
real, parameter :: CLAST_STOP_TIME = 12*60	! max time for an OC to be blocked
real, parameter :: DT_FAST_MOVE = 10.0
real, parameter :: CLAST_RADIUS_FACTOR = 0.5	! relates OC radius to the sqrt of the number of monocytes

! Pit parameters
real, parameter :: LACUNA_A = 20		! Parameters of the elliptical region to be excavated
real, parameter :: LACUNA_B = 10
real, parameter :: MAX_PIT_DEPTH = 3	! grids (should be input parameter)
real, parameter :: OC_SIGNAL_SENSING_RANGE = 20	! grids
real, parameter :: Kattraction = 4

! Osteoblast parameters
real, parameter :: BLAST_PER_UM3 = 7.0e-6
real, parameter :: OB_SIGNAL_RADIUS = 5		! radius of disk over which OB integrates signal
real, parameter :: OB_REACH = 40			! to make contact with a monocyte (microns)
real, parameter :: BLAST_DWELL_TIME = 2*60

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
    integer :: lastmovestep
end type

type osteoblast_type
    integer :: ID
    integer :: site(3)
    integer :: iclump
    integer :: step
    integer(2) :: status
	integer(2) :: lastdir
    real :: movetime			! time that the cell will be checked for move
    real :: entrytime			! time that the cell entered the marrow (from blood or stem cell division)
    real :: dietime             ! time cell will die
end type

type pit_type
	integer :: delta(3)
	real :: rate
	real :: fraction
end type

type osteoclast_type
    integer :: ID
    integer :: site(3)
    real :: cm(3), radius
	real :: normal(3)
    integer(2) :: status
    integer(2) :: lastdir
    real :: fusetime			! time that the monocytes began to fuse to form the osteoclast
    real :: entrytime			! time that the cell became a mature osteoclast
    real :: movetime			! time that the cell will be checked for move
    real :: blocktime			! time that the cell became blocked
    real :: dietime             ! time cell will die
	integer :: count
!	integer :: mono(100)
	integer :: targetsite(3)
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
!    integer :: signal
!    real :: intensity
    real :: bone_fraction
end type

type surface_type
	real :: target_depth	! desired excavation depth (units of grids)
	real :: signal			! current bone signal strength
	real :: depth			! current excavation depth (units of grids)
	integer :: iclast
end type

type signal_type
	logical :: active
	integer :: site(3)
	real :: intensity
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
	integer :: ID
	integer :: iblast
	integer :: ncells
	integer :: list(MAX_CLUMP_CELLS)
	integer :: status
	real :: starttime
	real :: fusetime
	real :: cm(3)
end type

type patch_type
	integer :: x0, z0
	real :: a, b
	real :: volume
end type

end module