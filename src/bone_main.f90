!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
program main
use bone_mod
!integer :: nsteps = 100000 
character*(64) :: infile = 'basecase.inp'

use_tcp = .false.
call execute(infile)
end program
