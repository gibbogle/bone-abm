!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine get_ndays(infile,ndays)
character*(*) :: infile
integer :: ndays
real :: val
character*(32) :: str
integer :: nfinp = 10
character*(64) :: line

open(nfinp,file=infile,status='old')
do
	read(nfinp,'(a)') line
	read(line,*) val,str
	if (trim(str) == 'NDAYS') then
		ndays = val
		close(nfinp)
		return
	endif
enddo
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
program main
use bone_mod
integer :: nt = 10000
integer :: ndays = 2
character*(64) :: infile = 'basecase.inp'
integer :: res

use_tcp = .false.
call get_ndays(infile,ndays)
write(*,*) 'NDAYS: ',ndays
call execute(infile)

nt = ndays*24*60*4
do it = 1,nt
	call simulate_step(res)
enddo
call terminate_run(res)
end program
