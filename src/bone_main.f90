!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine get_days(infile,days)
character*(*) :: infile
real :: days
real :: val
character*(32) :: str
integer :: nfinp = 10
character*(64) :: line

open(nfinp,file=infile,status='old')
do
	read(nfinp,'(a)') line
	read(line,*) val,str
	if (trim(str) == 'NDAYS') then
		days = val
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
real :: simdays = 2
character*(64) :: infile = 'basecase.inp'
integer :: res, nlen

use_tcp = .false.
call get_days(infile,simdays)
write(*,*) 'DAYS: ',simdays
nlen = len(infile)
call execute(infile,nlen)

nt = simdays*24*60*4
do it = 1,nt
	call simulate_step(res)
enddo
call terminate_run(res)
end program
