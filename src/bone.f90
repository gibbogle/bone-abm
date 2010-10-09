
module bone_mod
use global
use fields
use motion

implicit none

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine rng_initialisation()
integer, allocatable :: zig_seed(:) 
integer :: i
integer :: npar, grainsize = 32

npar = Mnodes
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
deallocate(zig_seed)
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine setup(ok)
logical :: ok
integer :: x, y, z, del, ix, iy, i, site(3)
real :: xc, yc, fac, xx, yy
integer :: kpar = 0
real :: R

ok = .true.
call rng_initialisation
call make_reldir

allocate(occupancy(NX,NY,NZ))
occupancy%region = MARROW
occupancy%indx = 0
occupancy%signal = 0
occupancy%intensity = 0
do y = 1,NBY
	occupancy(:,y,:)%region = BONE
enddo
!write(*,*) 'z=1:  x,y: ',cap(1,1),cap(1,2)
!write(*,*) 'z=NZ: x,y: ',cap(2,1),cap(2,2)
del = capR+1
do z = 1,NZ
	fac = (z-1.)/(NZ-1.)
	xc = cap(1,1)*(1-fac) + cap(2,1)*fac
	yc = cap(1,2)*(1-fac) + cap(2,2)*fac
	do ix = -del,del
		do iy = -del,del
			x = xc + ix + 0.5
			y = yc + iy + 0.5
			if ((x-xc)**2 + (y-yc)**2 <= capR**2) then
				occupancy(x,y,z)%region = BLOOD
!				write(*,'(3i4)') x,y,z
			endif
		enddo
	enddo
enddo

nclast = 0
nmono = 0
allocate(mono(MAX_MONO))
allocate(clast(MAX_CLAST))
do while (nmono < NMONO_INITIAL)
	x = random_int(1,NX,kpar)
	y = random_int(NBY+1,NY,kpar)
	z = random_int(1,NZ,kpar)
	if (occupancy(x,y,z)%region /= MARROW) cycle
	site = (/x,y,z/)
	call addMono(site)
enddo

!nstem = 500
allocate(stem(NSTEM))
i = 0
do while (i < NSTEM)
	x = random_int(1,NX,kpar)
	y = random_int(NBY+1,NY,kpar)
	z = random_int(1,NZ,kpar)
	if (occupancy(x,y,z)%region /= MARROW) cycle
	i = i + 1
	stem(i)%ID = i
	stem(i)%site = (/x,y,z/)
	call random_number(R)
	stem(i)%dividetime = R*STEM_CYCLETIME
	occupancy(x,y,z)%species = STEMCELL
enddo

nsignal = 1
signal(1)%site = (/ NX/2,NBY,NZ/2 /)
call setSignal(1,ON,ok)
if (ok == .false.) return

if (S1P_chemotaxis) then
	call init_fields
endif
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine updater
real :: S
integer :: i, k, irel, dir, region, kcell, site(3)
real :: tnow, R

tnow = istep*DELTA_T
do i = 1,nmono
	if (mono(i)%region /= MARROW) cycle
	S = mono(i)%S1P1
	S = S + DELTA_T*rate_S1P1(S)
	mono(i)%S1P1 = min(S,1.0)
enddo
do i = 1,NSTEM
	if (tnow > stem(i)%dividetime) then
		call random_number(R)
		irel = nreldir*R
		do k = 1,nreldir
			irel = irel + 1
			if (irel > nreldir) irel = 1
			dir = reldir(1,irel)
			site = stem(i)%site + jumpvec(:,dir)
			if (free_site(site,region,kcell)) then
				call addMono(site)
!				write(*,*) 'added monocyte: ',i,site
				stem(i)%dividetime = tnow + STEM_CYCLETIME
				exit
			endif
		enddo
	endif
enddo
do i = 1,nclast
	if (clast(i)%status /= ALIVE) cycle
	if (tnow > clast(i)%dietime) then
		call clastDeath(i)
		cycle
	endif
	do k = 1,clast(i)%npit
		if (clast(i)%pit(k)%fraction > 0) then
			clast(i)%pit(k)%fraction = clast(i)%pit(k)%fraction - clast(i)%pit(k)%rate*DELTA_T
			if (clast(i)%pit(k)%fraction <= 0) then
				site = clast(i)%pit(k)%site
				occupancy(site(1),site(2),site(3))%region = PIT
				occupancy(site(1),site(2),site(3))%indx = 0
				if (site(2) > 1) then
					write(logmsg,*) 'Created pit site: ',site
					call logger(logmsg)
					clast(i)%pit(k)%site(2) = site(2) - 1
					clast(i)%pit(k)%fraction = 1
				else
					clast(i)%pit(k)%fraction = 0
				endif
			endif
		endif
	enddo
enddo
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine addMono(site)
integer :: site(3)
integer :: kpar = 0

if (nmono == MAX_MONO) then
	call logger('Error: addMono: Exceeded monocyte array, increase MAX_MONO')
	stop
endif
nmono = nmono + 1
mono_cnt = mono_cnt + 1
mono(nmono)%ID = nmono
mono(nmono)%site = site
mono(nmono)%region = MARROW
mono(nmono)%status = ALIVE
mono(nmono)%lastdir = random_int(1,6,kpar)
mono(nmono)%S1P1 = 0
occupancy(site(1),site(2),site(3))%species = MONOCYTE
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
real function rate_S1P1(S)
real :: S
rate_S1P1 = S1P1_BASERATE
end function

!------------------------------------------------------------------------------------------------
! To avoid problems with wrapping, no signal site is allowed to occur near the boundary.
! The bone normal is used to give higher preference to sites near the bone.
!------------------------------------------------------------------------------------------------
subroutine setSignal(isig,onoff,ok)
logical :: ok
integer :: isig, onoff
integer :: site(3)
integer :: x1, x2, y1, y2, z1, z2
integer :: x, y, z
real :: R2, d2, vn(3), v(3)

ok = .true.
site = signal(isig)%site
x1 = site(1) - (SIG_RADIUS+1)
x2 = site(1) + (SIG_RADIUS+1)
y1 = site(2) - (SIG_RADIUS+1)
y2 = site(2) + (SIG_RADIUS+1)
z1 = site(3) - (SIG_RADIUS+1)
z2 = site(3) + (SIG_RADIUS+1)
R2 = SIG_RADIUS**2
if (onoff == ON) then
	if (x1 < 1 .or. x2 > NX .or. y1 < 1 .or. y2 > NY .or. z1 < 1 .or. z2 > NZ) then
		write(logmsg,*) 'Error: setSignal: signal site too close to boundary: ',site
		call logger(logmsg)
		ok = .false.
		return
	endif
	signal(isig)%active = .true.
	signal(isig)%site = site
	call boneNormal(site,vn)
	signal(isig)%normal = vn
else
	signal(isig)%active = .false.
	signal(isig)%site = (/0,0,0/)
endif
do x = x1,x2
	do y = y1, y2
		do z = z1,z2
			d2 = (x-site(1))**2 + (y-site(2))**2 + (z-site(3))**2
			if (d2 <= R2) then
				if (occupancy(x,y,z)%region == MARROW) then
					if (onoff == ON) then
						v = (/x,y,z/) - site
						occupancy(x,y,z)%signal = isig
						occupancy(x,y,z)%intensity = 1/(d2 + dot_product(v,vn))
					else
						occupancy(x,y,z)%signal = 0
						occupancy(x,y,z)%intensity = 0
					endif
				endif
			endif
		enddo
	enddo
enddo
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine checkSignals(ok)
logical :: ok
integer :: isig
integer :: site(3)
integer :: x1, x2, y1, y2, z1, z2
integer :: x, y, z, n, nt
real :: R2, d2
type(occupancy_type), pointer :: p

ok = .true.
R2 = 0.33*(SIG_RADIUS)**2
do isig = 1,nsignal
	if (.not.signal(isig)%active) cycle
	site = signal(isig)%site
	x1 = site(1) - (SIG_RADIUS+1)
	x2 = site(1) + (SIG_RADIUS+1)
	y1 = site(2) - (SIG_RADIUS+1)
	y2 = site(2) + (SIG_RADIUS+1)
	z1 = site(3) - (SIG_RADIUS+1)
	z2 = site(3) + (SIG_RADIUS+1)
	n = 0
	nt = 0
	do x = x1,x2
		do y = y1, y2
			do z = z1,z2
				p => occupancy(x,y,z)
!				d2 = (x-site(1))**2 + (y-site(2))**2 + (z-site(3))**2
!				if (d2 <= R2) then
				if (p%signal == isig .and. p%intensity >= FTHRESHOLD) then
					if (p%region == MARROW) then
						nt = nt + 1
						if (p%indx /= 0) then
							n = n+1
						endif
					endif
				endif
			enddo
		enddo
	enddo
!	write(*,*) 'Total near sites: ',nt
	if (n >= 0.8*nt) then
		write(logmsg,*) 'fusing monocytes: ',n,nt
		call logger(logmsg)
		call fuser(isig,n,ok)
		if (ok == .false.) return
		call setSignal(isig,OFF,ok)
		if (ok == .false.) return
		signal(isig)%active = .false.
	endif
enddo
end subroutine

!------------------------------------------------------------------------------------------------
! For now it is assumed, for simplicity, that the bone surface lies in an X-Z plane,
! i.e. that it is normal to the Y axis.  This makes it easy to determine which
! bone sites are subject to resorption.
!------------------------------------------------------------------------------------------------
subroutine fuser(isig,n,ok)
integer :: isig,n
logical :: ok
integer :: site(3)
integer :: x1, x2, y1, y2, z1, z2
integer :: x, y, z, cnt, imin, npit, i, yb
real :: tnow, d
integer :: bonesite(3,50)
logical :: inlist
type(occupancy_type), pointer :: p

ok = .true.
tnow = istep*DELTA_T
nclast = nclast + 1
site = signal(isig)%site
x1 = site(1) - (SIG_RADIUS+1)
x2 = site(1) + (SIG_RADIUS+1)
y1 = site(2) - (SIG_RADIUS+1)
y2 = site(2) + (SIG_RADIUS+1)
z1 = site(3) - (SIG_RADIUS+1)
z2 = site(3) + (SIG_RADIUS+1)
cnt = 0
npit = 0
do x = x1,x2
	do y = y1, y2
		do z = z1,z2
			p => occupancy(x,y,z)
			if (p%signal == isig .and. p%intensity >= FTHRESHOLD &
				.and. p%region == MARROW .and. p%indx /= 0) then
				cnt = cnt+1
				clast(nclast)%mono(cnt) = p%indx
				mono(p%indx)%status = FUSED
				yb = y
				do
					yb = yb - 1
					if (yb < 1) then
						write(logmsg,*) 'Error: fuser: yb < 1'
						call logger(logmsg)
						ok = .false.
						return
					endif
					if (occupancy(x,yb,z)%region == BONE) then
						inlist = .false.
						do i = 1,npit
							if (bonesite(1,i) == x .and. bonesite(3,i) == z) then
								inlist = .true.
								exit
							endif
						enddo
						if (.not.inlist) then
							npit = npit + 1
							bonesite(:,npit) = (/x,yb,z/)
						endif
						exit
					endif
				enddo
				write(logmsg,'(5i6)') cnt,p%indx,x,y,z
				call logger(logmsg)
			endif
		enddo
	enddo
enddo
clast(nclast)%ID = nclast
clast(nclast)%site = site
clast(nclast)%normal = signal(isig)%normal
clast(nclast)%status = ALIVE
clast(nclast)%entrytime = tnow
clast(nclast)%dietime = tnow + clastLifetime()
clast(nclast)%count = cnt
clast(nclast)%npit = npit
allocate(clast(nclast)%pit(npit))
do i = 1,npit
	clast(nclast)%pit(i)%site = bonesite(:,i)
	clast(nclast)%pit(i)%fraction = 1.0
	d = sqrt(real((bonesite(1,i)-site(1))**2 + (bonesite(2,i)-site(2))**2 + (bonesite(3,i)-site(3))**2))
	clast(nclast)%pit(i)%rate = resorptionRate(cnt,d)
enddo
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
real function resorptionRate(n,d)
integer :: n
real :: d
real :: f

if (d >= MAX_RESORPTION_D) then
	f = 0
else
	f = 1 - d/MAX_RESORPTION_D
endif
resorptionRate = MAX_RESORPTION_RATE*(real(n)/MAX_RESORPTION_N)*f
end function

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine clastDeath(i)
integer :: i
integer :: k, site(3), kcell

do k = 1,clast(i)%npit
	if (clast(i)%pit(k)%fraction > 0 .and. clast(i)%pit(k)%fraction < 0.5) then
		site = clast(i)%pit(k)%site
		occupancy(site(1),site(2),site(3))%region = PIT
		occupancy(site(1),site(2),site(3))%indx = 0
	endif
enddo
clast(i)%status = DEAD
deallocate(clast(i)%pit)
do k = 1,clast(i)%count
	kcell = clast(i)%mono(k)
	site = mono(kcell)%site
	mono(kcell)%status = DEAD
	occupancy(site(1),site(2),site(3))%indx = 0
enddo
write(logmsg,*) 'clast death: ',i
call logger(logmsg)
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
real function clastLifetime()

clastLifetime = MEAN_CLAST_LIFETIME
end function

!------------------------------------------------------------------------------------------------
! Determine approximate normal to the bone surface at site().
! Method: inspect all neighbours of site, average the (x,y,z) of locations in marrow.
! The normal vector is from site to this average point.
!------------------------------------------------------------------------------------------------
subroutine boneNormal(site,v)
integer :: site(3)
real :: v(3)
integer :: k, irel, cnt, site1(3)

v = 0
cnt = 0
do irel = 1,nreldir
	k = reldir(1,irel)
	site1 = site + jumpvec(:,k)
	if (occupancy(site1(1),site1(2),site1(3))%region == MARROW) then
		cnt = cnt + 1
		v = v + site1
	endif
enddo
v = v/cnt - site
v = v/rnorm(v)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine save_pos(ok)
logical :: ok
integer :: error
character*(64) :: msg

ok = .true.
if (use_TCP) then
	call save_cell_positions
	msg = 'VTK'
	clear_to_send = .false.
    call winsock_send(awp_1,msg,len_trim(msg),error)
    if (error /= 0) ok = .false.
endif
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine save_cell_positions
!!!use ifport
integer :: k, kcell, site(3), j, nfused, x, y, z
real :: mono_state
!integer :: itcstate, stype, ctype
real :: clast_diam = 0.9
real :: mono_diam = 0.5
logical :: ex
character*(12) :: fname = 'cell_pos.dat'
character*(9) :: removefile = 'TO_REMOVE'

if (simulation_start) then
	inquire(file=fname,exist=ex)
	if (ex) then
		call unlink(fname)
	endif
	inquire(file=removefile,exist=ex)
	if (ex) then
		call unlink(removefile)
	endif
endif
simulation_start = .false.

if (.not.clear_to_send) then
	! wait until the file called removefile exists, then remove it
	inquire(file=removefile,exist=ex)
	if (.not.ex) then
!		call logger('wait')
		do
			inquire(file=removefile,exist=ex)
			if (.not.ex) then
!				call millisleep(10) ! no good at all
			else
				exit
			endif
		enddo
	endif
	call unlink(removefile)
	clear_to_send = .true.
endif

open(nfpos,file=fname,status='new')
! Monocyte section
if (nmono > 0) then
	nfused = 0
    do kcell = 1,nmono
		if (mono(kcell)%status == DEAD) cycle
		if (mono(kcell)%region /= MARROW) cycle
		if (mono(kcell)%status == FUSED) then
			mono_state = 1.0
			nfused = nfused + 1
		else
            mono_state = 0.0
		endif
        site = mono(kcell)%site
        write(nfpos,'(a2,i4,3i4,f4.1,f5.2)') 'M ',kcell-1, site, mono_diam, mono_state
    enddo
    if (nfused > 0) then
	    write(logmsg,*) 'nfused: ',nfused
		call logger(logmsg)
	endif
endif

! Capillary section
write(nfpos,'(a,6i4,f5.2)') 'C ',cap(1,1),cap(1,2),1,cap(2,1),cap(2,2),NZ,capR

! Pit section
do x = 1,NX
	do z = 1,NZ
		do y = 1,NBY 
			if (occupancy(x,y,z)%region == PIT) then
				write(nfpos,'(a,3i4)') 'P ',x,y,z
				exit
			endif
		enddo
	enddo
enddo

! T cell section
!do kcell = 1,nclast
!	if (cellist(kcell)%ID == 0) cycle  ! gap
!	site = cellist(kcell)%site
!	bnd = cellist(kcell)%DCbound
!	ctype = cellist(kcell)%ctype
!	stype = struct_type(ctype)
!	if (stype == NONCOG_TYPE_TAG) then
!		itcstate = -1
!	else
!		gen = get_generation(cellist(kcell)%cptr)
!		if (get_stage(cellist(kcell)%cptr) == NAIVE) then
!			itcstate = 0
!		else
!			if (bnd(1) == 0 .and. bnd(2) == 0) then
!				itcstate = gen
!			else
!				itcstate = 99
!			endif
!		endif
!	endif
!	! Need tcstate to convey non-activated status, i.e. 0 = non-activated
!	write(nfpos,'(a2,i6,3i4,f4.1,i3)') 'T ',kcell-1, site, Tcell_diam, itcstate
!enddo
! Bond section
!do kcell = 1,nlist
!	if (cellist(kcell)%ID == 0) cycle  ! gap
!	site = cellist(kcell)%site
!	do j = 1,2
!		idc = cellist(kcell)%DCbound(j)
!		if (idc /= 0) then
!			if (DClist(idc)%capable) then
!				dcsite = DClist(idc)%site
!				dcstate = 1
!				write(nfpos,'(a2,2i5)') 'B ',kcell-1,idc-1
!			endif
!		endif
!	enddo
!enddo
write(nfpos,'(a2,i6)') 'E ',istep
close(nfpos)

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_pause
logical :: paused

inquire(file=pausefile,exist=paused)
if (paused) then
	call logger('Pause order received')
	do while (paused)
        call sleeper(1)   ! Too coarse!
		inquire(file=pausefile,exist=paused)
	enddo
	call logger('Resuming ...')
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine connection(awp,port,error)
TYPE(winsockport) :: awp
integer :: port, error
integer :: address = 0
!!!character*(64) :: ip_address = "127.0.0.1"C      ! need a portable way to make a null-terminated C string
character*(64) :: host_name = "localhost"

if (.not.winsock_init(1)) then
    call logger("winsock_init failed")
    stop
endif
!write(nftemp,*) 'did winsock_init'

awp%handle = 0
awp%host_name = host_name
awp%ip_port = port
awp%protocol = IPPROTO_TCP
call Set_Winsock_Port (awp,error)

if (.not.awp%is_open) then
    write(nflog,*) 'Error: connection: awp not open: ',port
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
call connection(awp_0,TCP_PORT_0,error)
if (awp_0%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_0 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_0%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_0'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_0  '
call logger(logmsg)

call connection(awp_1,TCP_PORT_1,error)
if (awp_1%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_1 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_1%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_1'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_1  '
call logger(logmsg)
end subroutine


!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine simulate(nsteps,ok)
integer :: nsteps
logical :: ok
integer :: error

clear_to_send = .true.
simulation_start = .true.

ok = .true.
do istep = 1,nsteps
	inquire(file=stopfile,exist=stopped)
	if (stopped) then
		write(logmsg,'(a)') 'Stop order received'
		call logger(logmsg)
		return
	endif
	call check_pause	
	if (mod(istep,100) == 0) then
		write(logmsg,*) 'istep: ',istep,mono_cnt
		call logger(logmsg)
	endif
	if (mod(istep,100) == 0) then
		call save_pos(ok)
		if (.not.ok) return
	endif
	call updater
	call mover
	if (nsignal > 0) then
		call checkSignals(ok)
		if (.not.ok) return
	endif
enddo
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine terminate
call par_zigfree
if (allocated(occupancy)) deallocate(occupancy)
if (allocated(mono)) deallocate(mono)
if (allocated(clast)) deallocate(clast)
if (allocated(stem)) deallocate(stem)
close(nflog)
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine execute(nsteps)
!DEC$ ATTRIBUTES DLLEXPORT :: EXECUTE
!DEC$ ATTRIBUTES C, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"EXECUTE" :: execute
integer :: nsteps
logical :: ok

open(nflog,file='bone.log',status='replace')
if (use_tcp) then
	call logger('call connecter')
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
	call logger('did connecter')
endif
call setup(ok)
if (ok) then
	call simulate(nsteps,ok)
	call logger('Ended simulation')
	call logger('Execution successful')
else
	call logger('setup failed')
endif
call terminate
end subroutine

end module

