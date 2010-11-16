!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
program main
use bone_mod
integer :: nt = 10000
integer :: ndays = 2
character*(64) :: infile = 'basecase.inp'
integer :: res

use_tcp = .false.
call execute(infile)

nt = ndays*24*60*4
do it = 1,nt
	call simulate_step(res)
enddo
end program
