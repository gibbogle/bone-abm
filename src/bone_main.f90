!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
program main
use bone_mod
integer :: nsteps = 100000 

use_tcp = .false.
call execute(nsteps)
end program
