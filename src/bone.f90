
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
subroutine setup(infile,ok)
logical :: ok
character*(*) :: infile
integer :: x, y, z, del, dx, dy, dz, ix, i, site(3), ndiv, icap, k
real :: xc, yc, zc, fac, x2, y2, z2
integer :: kpar = 0
real :: R, alfa

ok = .true.
inputfile = infile
call read_inputfile(ok)
if (.not.ok) return
call rng_initialisation
call make_reldir

write(logmsg,*) 'NX,NY,NZ: ',NX,NY,NZ
call logger(logmsg)
allocate(occupancy(NX,NY,NZ))
occupancy%region = MARROW
occupancy%indx = 0
occupancy%signal = 0
occupancy%intensity = 0
do y = 1,NBY
	occupancy(:,y,:)%region = BONE
enddo

! Create a test capillary
ncap = 1
allocate(capillary(ncap))
capillary(1)%radius = capR
capillary(1)%pos1 = (/ 0.5, NY/2., NZ/2. /)
capillary(1)%pos2 = (/ NX+0.5, NY/2., NZ/2. /)

! Set up capillary sites - this is currently VERY CRUDE
! If the centre of a site (cube) falls inside the capillary tube,
! the site is tagged BLOOD.  
! Note that (xc,yc,zc), which lies on the capillary centreline,
! is in (x,y,z) coord values, offset from site index values by 0.5
! E.g. if site(:) = (/2, 3, 4/), corresponding (x,y,z) = (/1.5, 2.5, 3.5/)
! This has the origin of the axis system at the lower left rear corner
! of the region.  (Note that in OpenGL/VTK Z axis is out of the screen.)

ndiv = 100
do icap = 1,ncap
	del = capillary(icap)%radius + 2
	do k = 1,ndiv
		alfa = (k-1.0)/(ndiv-1.0)
		xc = (1-alfa)*capillary(icap)%pos1(1) + alfa*capillary(icap)%pos2(1)
		yc = (1-alfa)*capillary(icap)%pos1(2) + alfa*capillary(icap)%pos2(2)
		zc = (1-alfa)*capillary(icap)%pos1(3) + alfa*capillary(icap)%pos2(3)		
		do dx = -del,del
			x = xc + dx
			if (x < 1 .or. x > NX) cycle
			do dy = -del,del
				y = yc + dy
				if (y <= NBY .or. y > NY) cycle
				do dz = -del,del
					z = zc + dz
					if (z < 1 .or. z > NZ) cycle
!					x2 = (x - 0.5 - xc)*(x - 0.5 - xc)
!					y2 = (y - 0.5 - yc)*(y - 0.5 - yc)
!					z2 = (z - 0.5 - zc)*(z - 0.5 - zc)
					x2 = (x - xc)*(x - xc)
					y2 = (y - yc)*(y - yc)
					z2 = (z - zc)*(z - zc)
					! These are the squared distances, in three axis directions, from the
					! site (cube) midpoint to the capillary centreline location (xc,yc,zc)
					if (x2+y2+z2 <= (capillary(icap)%radius)**2) then
						occupancy(x,y,z)%region = BLOOD
						if (k <= 10) then
							write(nflog,'(4f6.1,3i4)') xc,yc,zc,sqrt(x2+y2+z2),x,y,z
						endif
					endif
				enddo
			enddo
		enddo
	enddo
enddo
NMONO_INITIAL = (NX*(NY-NBY)*NZ)/40
nclast = 0
nmono = 0
nborn = 0
nleft = 0
mono_cnt = 0
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
	stem(i)%dividetime = R*STEM_CYCLETIME	! stem cells are due to divide at random times
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
integer :: i, k, iclast, irel, dir, region, kcell, site(3)
real :: tnow, R
type(osteoclast_type), pointer :: pclast

tnow = istep*DELTA_T
! Monocyte S1P1 level grows
do i = 1,nmono
	if (mono(i)%region /= MARROW) cycle
	S = mono(i)%S1P1
	S = S + DELTA_T*rate_S1P1(S)
	mono(i)%S1P1 = min(S,1.0)
enddo
! Stem cells divide
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
				nborn = nborn + 1
				exit
			endif
		enddo
	endif
enddo
! Osteoclasts complete fusing, dissolve bone, or die.
do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status == DEAD) cycle
	if (pclast%status == ALIVE .and. tnow > pclast%dietime) then
		call clastDeath(iclast)
		cycle
	endif
	if (pclast%status == FUSING) then
		if (tnow >= pclast%entrytime) then
			call completeFusing(iclast)
			cycle
		endif
	endif
	do k = 1,pclast%npit
		if (pclast%pit(k)%fraction > 0) then
			pclast%pit(k)%fraction = pclast%pit(k)%fraction - pclast%pit(k)%rate*DELTA_T
			if (pclast%pit(k)%fraction <= 0) then
				site = pclast%pit(k)%site
				occupancy(site(1),site(2),site(3))%region = PIT
				occupancy(site(1),site(2),site(3))%indx = 0
				if (site(2) > 1) then
					write(logmsg,*) 'Created pit site: ',site
					call logger(logmsg)
					pclast%pit(k)%site(2) = site(2) - 1
					pclast%pit(k)%fraction = 1
				else
					pclast%pit(k)%fraction = 0
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
mono(nmono)%status = MOTILE
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
! Each osteocyte(?) signal is checked to see if the number of monocytes that have gathered is
! sufficient to initiate osteoclastogenesis.  Currently only the number of monocytes within
! a specified volume is used to determine the initiation.  It would be 
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
! Was using 0.8*nt
	if (n >= MTHRESHOLD) then
		write(logmsg,*) 'fusing monocytes: ',n,nt
		call logger(logmsg)
		call startFusing(isig,n,ok)
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
subroutine startFusing(isig,n,ok)
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
site = signal(isig)%site
nclast = nclast + 1
clast(nclast)%ID = nclast
clast(nclast)%site = site
clast(nclast)%normal = signal(isig)%normal
clast(nclast)%status = FUSING
clast(nclast)%fusetime = tnow
clast(nclast)%entrytime = tnow + FUSING_TIME
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
				mono(p%indx)%iclast = nclast
				mono(p%indx)%status = FUSING
				mono(p%indx)%exittime = clast(nclast)%entrytime
				yb = y
				do
					yb = yb - 1
					if (yb < 1) then
						write(logmsg,*) 'Error: startFusing: yb < 1'
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
subroutine completeFusing(iclast)
integer :: iclast
integer :: i
real :: tnow

tnow = istep*DELTA_T
clast(iclast)%status = ALIVE
clast(iclast)%dietime = tnow + clastLifetime()
do i = 1,clast(iclast)%count
	mono(clast(iclast)%mono(i))%status = FUSED
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
integer :: k, kcell, iclast, site(3), j, nfused, x, y, z, status, mono_state, icap
real :: tnow, t1, t2, fraction
!integer :: itcstate, stype, ctype
real :: clast_diam = 0.9
real :: mono_diam = 0.5
logical :: ex
character*(12) :: fname = 'cell_pos.dat'
character*(9) :: removefile = 'TO_REMOVE'

tnow = istep*DELTA_T
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
		status = mono(kcell)%status
		if (status == DEAD) cycle
		if (mono(kcell)%region /= MARROW) cycle
		if (status == CROSSING) then
			mono_state = 1
		elseif (status == FUSING) then
			iclast = mono(kcell)%iclast
			t1 = clast(iclast)%fusetime
			t2 = clast(iclast)%entrytime
			fraction = min(1.0,(tnow-t1)/(t2-t1))
			mono_state = 2 + fraction*99
		elseif (status == FUSED) then
			mono_state = 101
			nfused = nfused + 1
		else
            mono_state = 0
		endif
        site = mono(kcell)%site
        write(nfpos,'(a2,i6,4i4)') 'M ',kcell-1, site, mono_state
        if (status == FUSING .or. status == FUSED) then
			write(nflog,*) istep, kcell, site, mono_state
		endif
    enddo
!    if (nfused > 0) then
!	    write(logmsg,*) 'nfused: ',nfused
!		call logger(logmsg)
!	endif
endif

! Capillary section
do icap = 1,ncap
	write(nfpos,'(a,i4,6f6.1,f5.2)') 'C ',icap,capillary(icap)%pos1,capillary(icap)%pos2,capillary(icap)%radius + 0.25
enddo

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

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine read_inputfile(ok)
logical :: ok
real :: TC_AVIDITY_MEAN,TC_AVIDITY_SHAPE,TC_STIM_RATE_CONSTANT,TC_STIM_HALFLIFE,divide_mean1,divide_shape1
real :: DC_LIFETIME_MEAN,DC_LIFETIME_SHAPE
real :: IL2_THRESHOLD,ACTIVATION_THRESHOLD,FIRST_DIVISION_THRESHOLD,DIVISION_THRESHOLD,EXIT_THRESHOLD,STIMULATION_LIMIT
real :: chemo_radius,chemo_K_exit,chemo_K_DC

ok = .false.
open(nfinp,file=inputfile,status='old')

read(nfinp,*) TC_AVIDITY_MEAN				! mean of avidity distribution (only if fix_avidity = false)
read(nfinp,*) TC_AVIDITY_SHAPE			    ! shape -> 1 gives normal dist with small variance
read(nfinp,*) TC_STIM_RATE_CONSTANT			! rate const for TCR stimulation (-> molecules/min)
read(nfinp,*) TC_STIM_HALFLIFE				! halflife of T cell stimulation (hours)
read(nfinp,*) divide_mean1
read(nfinp,*) divide_shape1
read(nfinp,*) BETA							! speed: 0 < beta < 1		(0.65)
read(nfinp,*) RHO							! persistence: 0 < rho < 1	(0.95)

read(nfinp,*) DC_LIFETIME_MEAN				! days
read(nfinp,*) DC_LIFETIME_SHAPE 			! days

read(nfinp,*) IL2_THRESHOLD					! stimulation needed to initiate IL-2/CD25 production
read(nfinp,*) ACTIVATION_THRESHOLD			! stimulation needed for activation
read(nfinp,*) FIRST_DIVISION_THRESHOLD		! activation level needed for first division
read(nfinp,*) DIVISION_THRESHOLD			! activation level needed for subsequent division
read(nfinp,*) EXIT_THRESHOLD				! activation level below which exit is permitted
read(nfinp,*) STIMULATION_LIMIT				! maximum activation level

read(nfinp,*) NX							! size of cubical region

read(nfinp,*) exit_rule						! 1 = no chemotaxis, 2 = chemotaxis
read(nfinp,*) exit_region					! region for cell exits 1 = capillary, 2 = sinusoid
read(nfinp,*) cross_prob					! probability (/timestep) of monocyte egress to capillary
read(nfinp,*) chemo_radius					! radius of chemotactic influence (sites)
read(nfinp,*) chemo_K_exit					! level of chemotactic influence towards exits
read(nfinp,*) chemo_K_DC					! level of chemotactic influence towards DCs

read(nfinp,*) days							! number of days to simulate
read(nfinp,*) seed(1)						! seed vector(1) for the RNGs
read(nfinp,*) seed(2)						! seed vector(2) for the RNGs
read(nfinp,*) ncpu							! number of threads, not used currently
read(nfinp,*) NT_GUI_OUT					! interval between GUI outputs (timesteps)
read(nfinp,*) SPECIES						! animal species
close(nfinp)
!call logger('Finished reading input file')

if (mod(NX,2) /= 0) NX = NX+1				! ensure that NX is even
NY = NX
NZ = NX
DELTA_X = MONOCYTE_DIAMETER
chemo_radius = chemo_radius/DELTA_X
!write(*,*) 'DC_RADIUS, chemo_radius: ',DC_RADIUS,chemo_radius
!chemo_N = max(3,int(chemo_radius + 0.5))	! convert from um to lattice grids
!write(*,*) 'chemo_N: ',chemo_N
! Note that currently exit chemotaxis and DC chemotaxis are treated in the same way - same decay
!chemo_exp = log(1/CHEMO_MIN)/log(chemo_radius)

!sigma = log(TC_AVIDITY_SHAPE)
!TC_AVIDITY_MEDIAN = TC_AVIDITY_MEAN/exp(sigma*sigma/2)
!sigma = log(DC_ANTIGEN_SHAPE)
!DC_ANTIGEN_MEDIAN = DC_ANTIGEN_MEAN/exp(sigma*sigma/2)
!sigma = log(DC_LIFETIME_SHAPE)
!DC_LIFETIME_MEDIAN = DC_LIFETIME_MEAN/exp(sigma*sigma/2)

!sigma = log(divide_shape1)
!divide_dist1%p1 = log(60*divide_mean1/exp(sigma*sigma/2))
!divide_dist1%p2 = sigma
!sigma = log(divide_shape2)
!divide_dist2%p1 = log(60*divide_mean2/exp(sigma*sigma/2))
!divide_dist2%p2 = sigma
!write(*,*) 'divide_dist2: ',divide_dist2

!if (chemo_K_exit == 0.0) then
!    ep_factor = 2.4
!elseif (chemo_K_exit <= 0.2) then
!    ep_factor = 2.3
!elseif (chemo_K_exit <= 0.4) then
!    ep_factor = 2.2
!elseif (chemo_K_exit <= 0.6) then
!    ep_factor = 2.15
!elseif (chemo_K_exit <= 0.8) then
!    ep_factor = 2.1
!else
!    ep_factor = 2.1
!endif

!call setup_dists

Nsteps = days*60*24/DELTA_T

!open(nfout,file=outputfile,status='replace')
!if (save_input) then
!    call save_inputfile(inputfile)
!    call save_parameters
!	call save_inputfile(fixedfile)
!endif
!write(*,*) 'open resultfile: ',resultfile

ok = .true.
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
subroutine snapshot
character*(128) :: msg
integer :: error

if (use_TCP) then
    if (.not.awp_1%is_open) then
        call logger("snapshot: awp_1 is not open")
        return
    endif
	write(msg,'(4i6)') istep,mono_cnt,nborn,nleft
!	call logger(msg)
    call winsock_send(awp_1,msg,len_trim(msg),error)
endif
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine simulate(ok)
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
	if (mod(istep,240) == 0) then
		call snapshot
!		write(logmsg,'(a,4i6)') 'istep: ',istep,mono_cnt,nmono,nleft
!		call logger(logmsg)
	endif
	if (mod(istep,NT_GUI_OUT) == 0) then
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
subroutine wrapup
call par_zigfree
if (allocated(occupancy)) deallocate(occupancy)
if (allocated(mono)) deallocate(mono)
if (allocated(clast)) deallocate(clast)
if (allocated(stem)) deallocate(stem)
if (allocated(capillary)) deallocate(capillary)
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine terminate(success)
logical :: success
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

call wrapup

if (success) then
	call logger(' Execution successful')
else
	call logger('  === Execution failed ===')
	call sleeper(1)
endif
close(nflog)

if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
!	    call logger("closed PORT_0")
	    call winsock_send(awp_1,quit,8,error)
	    call winsock_close(awp_1)
!	    call logger("closed PORT_1")
	endif
endif

end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine execute(infile)
!DEC$ ATTRIBUTES DLLEXPORT :: EXECUTE
!DEC$ ATTRIBUTES C, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"EXECUTE" :: execute
character*(*) :: infile
logical :: ok

open(nflog,file='bone.log',status='replace')
awp_0%is_open = .false.
awp_1%is_open = .false.

#ifdef GFORTRAN
    write(logmsg,'(a)') 'Built with GFORTRAN'
	call logger(logmsg)
#endif

logmsg = 'OS??'
#ifdef LINUX
    write(logmsg,'(a)') 'OS is Linux'
#endif
#ifdef OSX
    write(logmsg,'(a)') 'OS is OS-X'
#endif
#ifdef _WIN32
    write(logmsg,'(a)') 'OS is Windows'
#endif
#ifdef WINDOWS
    write(logmsg,'(a)') 'OS is Windows'
#endif
call logger(logmsg)

!#ifdef OPENMP
#if defined(OPENMP) || defined(_OPENMP)
    write(logmsg,'(a)') 'Executing with OpenMP'
	call logger(logmsg)
#endif

write(logmsg,*) 'inputfile:  ', infile
call logger(logmsg)
if (use_tcp) then
	call logger('call connecter')
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
	call logger('did connecter')
endif
call setup(infile,ok)
if (ok) then
	call simulate(ok)
	call logger('Ended simulation')
	call logger('Execution successful')
else
	call logger('setup failed')
endif
call terminate(ok)
end subroutine

end module

