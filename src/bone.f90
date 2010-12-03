
module bone_mod
!!DEC$ ATTRIBUTES DLLEXPORT :: BONE_MOD
use, intrinsic :: ISO_C_binding
use global
use fields
use motion

implicit none 
save

integer(c_int),BIND(C) :: success = 12345
!DEC$ ATTRIBUTES DLLEXPORT :: success

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
integer :: x, y, z, del, dx, dy, dz, ix, i, site(3), ndiv, icap, k, na, nbr(3)
integer, allocatable :: templist(:,:)
real :: xc, yc, zc, fac, x2, y2, z2, v(3)
integer :: kpar = 0
real(8) :: R
real :: alfa

ok = .true.
inputfile = infile
call read_inputfile(ok)
if (.not.ok) return
call rng_initialisation
PI = 4*atan(1.0)
call make_reldir

write(logmsg,*) 'NX,NY,NZ: ',NX,NY,NZ
call logger(logmsg)
allocate(occupancy(NX,NY,NZ))
occupancy%region = MARROW
occupancy%indx = 0
occupancy%signal = 0
occupancy%intensity = 0
occupancy%bone_fraction = 0
do y = 1,NBY
	occupancy(:,y,:)%region = BONE
	occupancy(:,y,:)%bone_fraction = 1.0
enddo

! Create a test capillary
ncap = 1
allocate(capillary(ncap))
capillary(1)%radius = CAPILLARY_DIAMETER/2
capillary(1)%pos1 = (/ 0.5, NY/2., NZ/2. /)
capillary(1)%pos2 = (/ NX+0.5, NY/2., NZ/2. /)
v = capillary(1)%pos2 - capillary(1)%pos1
capillary(1)%length = rnorm(v)
capillary(1)%surface_area = PI*capillary(1)%radius**2*capillary(1)%length

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
! Now set up a list of all marrow sites that are adjacent to a blood site
allocate(templist(3,NX*NZ*10))
na = 0
do x = 2,NX-1
	do y = NBY+1,NY-1
		do z = 2,NZ-1
			if (occupancy(x,y,z)%region /= MARROW) cycle
			site = (/x,y,z/)
			do k = 1,6
				nbr = site + neumann(:,k)
				if (occupancy(nbr(1),nbr(2),nbr(3))%region == BLOOD) then
					na = na + 1
					templist(:,na) = site
					exit
				endif
			enddo
		enddo
	enddo
enddo
allocate(entrysite(3,na))
entrysite(:,1:na) = templist(:,1:na)
deallocate(templist)
nentrysites = na

nsignal = 1
signal(1)%site = (/ NX/2,NBY,NZ/2 /)
call init_fields

NMONO_INITIAL = (NX*(NY-NBY)*NZ*DELTA_X**3/1.0e9)*MONO_PER_MM3	! domain as fraction of 1 mm3 x rate of monocytes
NSTEM = (PI*NX*CAPILLARY_DIAMETER*DELTA_X**2/1.0e6)*STEM_PER_MM2	! capillary surface area as fraction of 1 mm2 x rate of stem cells
write(logmsg,*) 'NSTEM, NMONO_INITIAL: ',NSTEM,NMONO_INITIAL
call logger(logmsg)
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
	call addMono('initial',site)
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
!	call random_number(R)
	R = par_uni(kpar)
	stem(i)%dividetime = R*STEM_CYCLETIME	! stem cells are due to divide at random times
	occupancy(x,y,z)%species = STEMCELL
enddo

nclump = 0

call setSignal(1,ON,ok)
if (.not.ok) return

end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine pitrates(pclast)
type(osteoclast_type), pointer :: pclast
integer :: ipit
real :: d, v(3)

do ipit = 1,pclast%npit
	v = pclast%pit(ipit)%site - pclast%site
	d = sqrt(dot_product(v,v))
	pclast%pit(ipit)%rate = resorptionRate(pclast%count,d)
enddo
end subroutine

!------------------------------------------------------------------------------------------------
! For now the rate of entry of osteoclast-precursor monocytes is an input variable, i.e. constant.
! In future it should be dependent on the level of osteoblast signalling, i.e. on work required.
! The inflow rate is specified as cells/hour, in_per_hour
!------------------------------------------------------------------------------------------------
subroutine influx
real :: rate, dr
real(8) :: R
integer :: i, n, kpar = 0

rate = in_per_hour*DELTA_T/60 
n = rate
dr = rate - n
do i = 1,n
	call mono_entry
enddo
if (dr > 0) then
!	call random_number(R)
	R = par_uni(kpar)
	if (R < dr) then
		call mono_entry
	endif
endif
end subroutine

!------------------------------------------------------------------------------------------------
! Choose a random location on the side of a capillary for monocyte initial location
!------------------------------------------------------------------------------------------------
subroutine old_mono_entry
integer :: icap, kpar=0
real(8) :: p(MAX_CAP), R
real :: d, v(3), alfa

do icap = 1,ncap
	p(icap) = capillary(icap)%surface_area
enddo
p(1:ncap) = p(1:ncap)/sum(p(1:ncap))

icap = random_selection(p,ncap)
R = par_uni(kpar)
d = 4 + R*(capillary(icap)%length - 8)		! leave out the ends
alfa = d/capillary(icap)%length
v = (1-alfa)*capillary(icap)%pos1 + alfa*capillary(icap)%pos2
end subroutine

!------------------------------------------------------------------------------------------------
! Choose a random location on the side of a capillary for monocyte initial location
!------------------------------------------------------------------------------------------------
subroutine mono_entry
integer :: i, site(3), kpar=0

do
	i = random_int(1,nentrysites,kpar)
	site = entrysite(:,i)
	if (occupancy(site(1),site(2),site(3))%indx == 0) exit
enddo
call addMono('blood',site)
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine update_S1P1(pmono)
type(monocyte_type), pointer :: pmono
real :: S

S = pmono%S1P1
S = S + rate_S1P1(S)*DELTA_T
pmono%S1P1 = min(S,1.0)
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine update_RANK(pmono)
type(monocyte_type), pointer :: pmono
real :: S, C
integer :: site(3)

S = pmono%RANKSIGNAL
site = pmono%site
C = RANKL_conc(site(1),site(2),site(3))
S = (1-RANKSIGNAL_decayfactor)*S + rate_RANKSIGNAL(S,C)*DELTA_T
pmono%RANKSIGNAL = min(S,1.0)
if (pmono%status == MOTILE .and. pmono%RANKSIGNAL > ST1) then
	pmono%status = CHEMOTACTIC
!	call logger('CHEMOTACTIC')
elseif (pmono%status == CHEMOTACTIC .and. pmono%RANKSIGNAL > ST2) then
	pmono%status = STICKY
!	call logger('STICKY')
endif
if (pmono%status == STICKY) then
	pmono%stickiness = pmono%RANKSIGNAL
endif
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine updater
real :: S
integer :: i, j, k, iclump, iclast, irel, dir, region, kcell, site(3), res, kpar=0
real :: tnow, stickysum
real(8) :: R
type(monocyte_type), pointer :: pmono
type(osteoclast_type), pointer :: pclast
type(occupancy_type), pointer :: pbone
type(clump_type), pointer :: pclump
logical :: on_surface

tnow = istep*DELTA_T
! Monocytes enter from the blood
call influx
! Monocyte S1P1 level grows, RANKL signal accumulates
do i = 1,nmono
	pmono => mono(i)
	if (pmono%status == LEFT) cycle
	if (pmono%region /= MARROW) cycle
	if (pmono%status < 1 .or. pmono%status > FUSING) cycle
	call update_S1P1(pmono)
	call update_RANK(pmono)
enddo
do iclump = 1,nclump
	pclump => clump(iclump)
	if (pclump%status < 0) cycle
! TESTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (pclump%status < FUSED) then
		call consolidate_clump(pclump)
	endif
	if (nclump > 1) then
		call separate_clump(pclump)
	endif
	if (pclump%status < FUSING) then
		stickysum = 0
		do i = 1,pclump%ncells
			stickysum = stickysum + mono(pclump%list(i))%stickiness
		enddo
		if (stickysum/pclump%ncells < ST2) then	! clump breaks up
			write(logmsg,*) 'Clump breaks up: ',iclump
			do i = 1,pclump%ncells
				mono(pclump%list(i))%status = CHEMOTACTIC
			enddo
			pclump%status = DEAD
			cycle
		endif
		if (pclump%ncells > CLUMP_THRESHOLD) then
			pclump%status = FUSING
			pclump%starttime = tnow
			pclump%fusetime = tnow + FUSING_TIME
			do i = 1,pclump%ncells
				mono(pclump%list(i))%status = FUSING
			enddo
		endif
	elseif (pclump%status == FUSING) then
		if (tnow >= pclump%fusetime) then
!			write(logmsg,*) '***fuse clump: ',iclump,pclump%fusetime,tnow
!			call logger(logmsg)
			call fuse_clump(pclump)
		endif
	elseif (pclump%status == FUSED) then
		call lower_clump(pclump,on_surface)
		if (on_surface) then
			call make_osteoclast(pclump)
		endif
	endif
enddo

! Stem cells divide
do i = 1,NSTEM
	if (tnow > stem(i)%dividetime) then
!		call random_number(R)
		R = par_uni(kpar)
		irel = nreldir*R
		do k = 1,nreldir
			irel = irel + 1
			if (irel > nreldir) irel = 1
			dir = reldir(1,irel)
			site = stem(i)%site + jumpvec(:,dir)
			if (free_site(site,region,kcell)) then
				call addMono('stem',site)
				stem(i)%dividetime = tnow + STEM_CYCLETIME
				nborn = nborn + 1
				exit
			endif
		enddo
	endif
enddo
! Osteoclasts complete fusing, dissolve bone, move, or die.
do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status == DEAD) cycle
	if (pclast%status == ALIVE .and. tnow > pclast%dietime) then
		call clastDeath(iclast)
		cycle
	endif
!	if (pclast%status == FUSING) then
!		if (tnow >= pclast%entrytime) then
!			call completeFusing(iclast)
!			cycle
!		endif
!	endif
	if (tnow > pclast%movetime) then
!		write(*,*) 'Move osteoclast: ',tnow,iclast
		call move_clast(pclast,res)
		if (res == 0) then
			pclast%movetime = tnow + CLAST_DWELL_TIME
		elseif (res == 1) then			! osteoclast needs to move fast to get to fresh bone
			pclast%movetime = tnow + CLAST_DWELL_TIME/10
			write(logmsg,*) 'Fast moving osteoclast: ',pclast%movetime,iclast
			call logger(logmsg)
		elseif (res == 2) then		! osteoclast can't move, must die
			call logger('Osteoclast cannot move - dies')
			call clastDeath(iclast)
			cycle
		endif
		call pitrates(pclast)
	endif
	do k = 1,pclast%npit
		site = pclast%pit(k)%site
		pbone => occupancy(site(1),site(2),site(3))
		if (pbone%bone_fraction > 0) then
			pbone%bone_fraction = pbone%bone_fraction - pclast%pit(k)%rate*DELTA_T
!			if (k == pclast%npit) then
!				write(*,*) 'bone_fraction: ',k,site,pbone%bone_fraction
!			endif
			if (pbone%bone_fraction <= 0) then
				pbone%bone_fraction = 0
				pbone%region = PIT
				pbone%indx = 0
				if (site(2) > 1) then
!					write(logmsg,*) 'Created pit site: ',site
!					call logger(logmsg)
!					write(*,*) 'Created pit site: ',site,'  bf: ',pbone%bone_fraction
					pclast%pit(k)%site(2) = site(2) - 1
!					pclast%pit(k)%fraction = 1
!				else
!					pclast%pit(k)%fraction = 0
				endif
			endif
		endif
	enddo
enddo
end subroutine

!------------------------------------------------------------------------------------------------
! After fusing is complete the clump settles down onto the bone surface and becomes an osteoclast.
! Currently no more monocytes can join a clump after it has fused.
!------------------------------------------------------------------------------------------------
subroutine fuse_clump(pclump)
type(clump_type) :: pclump
integer :: i, icm(3)

!call logger('fusing') 
pclump%status = FUSED
do i = 1,pclump%ncells
	mono(pclump%list(i))%status = FUSED
enddo
!icm = pclump%cm + 0.5
!pclump%cm = icm
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine addMono(source,site)
character*(*) :: source
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
mono(nmono)%iclump = 0
mono(nmono)%S1P1 = 0
mono(nmono)%RANKSIGNAL = 0
mono(nmono)%stickiness = 0
!nullify(mono(nmono)%clump)
occupancy(site(1),site(2),site(3))%species = MONOCYTE
!if (.not.use_TCP) then
!	write(*,*) 'added monocyte: ',nmono,site,'  ',source
!else
!	call logger('added monocyte  '//source)
!endif
end subroutine

!------------------------------------------------------------------------------------------------
! S1P1 (S1P receptor) expression is currently assumed to grow at a steady rate.
!------------------------------------------------------------------------------------------------
real function rate_S1P1(S)
real :: S
rate_S1P1 = S1P1_BASERATE
end function

!------------------------------------------------------------------------------------------------
! The rate of RANK signalling depends on the RANKL concentration (C) and the current integrated
! signal level (S).
!------------------------------------------------------------------------------------------------
real function rate_RANKSIGNAL(S,C)
real :: S, C
rate_RANKSIGNAL = RANKSIGNAL_rateconstant*(1-S)*C
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
x1 = site(1) - (SIGNAL_RADIUS+1)
x2 = site(1) + (SIGNAL_RADIUS+1)
y1 = site(2) - (SIGNAL_RADIUS+1)
y2 = site(2) + (SIGNAL_RADIUS+1)
z1 = site(3) - (SIGNAL_RADIUS+1)
z2 = site(3) + (SIGNAL_RADIUS+1)
R2 = SIGNAL_RADIUS**2
if (onoff == ON) then
	if (x1 < 1 .or. x2 > NX .or. z1 < 1 .or. z2 > NZ) then
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
! a specified volume is used to determine the initiation.  It would be better to employ some
! idea of stickiness.
! When the number of monocytes within the region defined by (signal intensity > SIGNAL_THRESHOLD)
! exceeds MTHRESHOLD, monocyte fusing is initiated.
! NOTE: This was the first crude model for fusing initiation
! NOT USED
!------------------------------------------------------------------------------------------------
subroutine checkSignals(ok)
logical :: ok
integer :: isig
integer :: site(3)
integer :: x1, x2, y1, y2, z1, z2
integer :: x, y, z, n, nt
!real :: R2, d2
type(occupancy_type), pointer :: p

ok = .true.
!R2 = 0.33*(SIGNAL_RADIUS)**2
do isig = 1,nsignal
	if (.not.signal(isig)%active) cycle
	site = signal(isig)%site
	x1 = site(1) - (SIGNAL_RADIUS+1)
	x2 = site(1) + (SIGNAL_RADIUS+1)
	y1 = site(2) - (SIGNAL_RADIUS+1)
	y2 = site(2) + (SIGNAL_RADIUS+1)
	z1 = site(3) - (SIGNAL_RADIUS+1)
	z2 = site(3) + (SIGNAL_RADIUS+1)
	n = 0
	nt = 0
	do x = x1,x2
		do y = y1, y2
			do z = z1,z2
				p => occupancy(x,y,z)
!				d2 = (x-site(1))**2 + (y-site(2))**2 + (z-site(3))**2
!				if (d2 <= R2) then
				if (p%signal == isig .and. p%intensity >= SIGNAL_THRESHOLD) then
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
		if (.not.ok) return
		call setSignal(isig,OFF,ok)
		if (.not.ok) return
		signal(isig)%active = .false.
	endif
enddo
end subroutine

!------------------------------------------------------------------------------------------------
! For now it is assumed, for simplicity, that the bone surface lies in an X-Z plane,
! i.e. that it is normal to the Y axis.  This makes it easy to determine which
! bone sites are subject to resorption.
! All monocytes within the high-signal zone (signal intensity >= SIGNAL_THRESHOLD) are joined
! to make an osteoclast.  The list of grid sites below the monocytes is created, and used
! to make the list of pits associated with the osteoclast.
!------------------------------------------------------------------------------------------------
subroutine startFusing(isig,n,ok)
integer :: isig,n
logical :: ok
integer :: site(3)
integer :: x1, x2, y1, y2, z1, z2
integer :: x, y, z, cnt, imin, npit, i, yb
real :: tnow, d
integer :: bonesite(3,100)
logical :: inlist
type(occupancy_type), pointer :: p
type(osteoclast_type), pointer :: pclast

ok = .true.
tnow = istep*DELTA_T
site = signal(isig)%site
nclast = nclast + 1
pclast => clast(nclast)
pclast%ID = nclast
pclast%site = site
pclast%normal = signal(isig)%normal
pclast%status = FUSING
pclast%fusetime = tnow
pclast%entrytime = tnow + FUSING_TIME
pclast%movetime = BIGTIME
x1 = site(1) - (SIGNAL_RADIUS+1)
x2 = site(1) + (SIGNAL_RADIUS+1)
y1 = site(2) - (SIGNAL_RADIUS+1)
y2 = site(2) + (SIGNAL_RADIUS+1)
z1 = site(3) - (SIGNAL_RADIUS+1)
z2 = site(3) + (SIGNAL_RADIUS+1)
cnt = 0
npit = 0
do x = x1,x2
	do y = y1, y2
		do z = z1,z2
			p => occupancy(x,y,z)
			if (p%signal == isig .and. p%intensity >= SIGNAL_THRESHOLD &
				.and. p%region == MARROW .and. p%indx /= 0) then
				cnt = cnt+1
				pclast%mono(cnt) = p%indx
				mono(p%indx)%iclast = nclast
				mono(p%indx)%status = FUSING
				mono(p%indx)%exittime = pclast%entrytime
				
				yb = y
				do
					yb = yb - 1
					if (yb < 1) then
						write(logmsg,*) 'Error: startFusing: yb < 1'
						call logger(logmsg)
						ok = .false.
						return
					endif
					if (occupancy(x,yb,z)%region == BONE .and. occupancy(x,yb,z)%bone_fraction > 0) then
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
				
!				write(logmsg,'(5i6)') cnt,p%indx,x,y,z
!				call logger(logmsg)
			endif
		enddo
	enddo
enddo

pclast%count = cnt
pclast%npit = npit
allocate(pclast%pit(npit))
do i = 1,npit
	pclast%pit(i)%site = bonesite(:,i)
enddo
call pitrates(pclast)
!	clast(nclast)%pit(i)%fraction = 1.0
!	d = sqrt(real((bonesite(1,i)-site(1))**2 + (bonesite(2,i)-site(2))**2 + (bonesite(3,i)-site(3))**2))
!	clast(nclast)%pit(i)%rate = resorptionRate(cnt,d)
!enddo
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine completeFusing(iclast)
integer :: iclast
integer :: i
real :: tnow

tnow = istep*DELTA_T
write(logmsg,*) 'completed fusing: ',iclast,tnow
call logger(logmsg)
clast(iclast)%status = ALIVE
clast(iclast)%movetime = tnow + CLAST_DWELL_TIME 
clast(iclast)%dietime = tnow + clastLifetime()
do i = 1,clast(iclast)%count
	mono(clast(iclast)%mono(i))%status = FUSED
enddo
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine make_osteoclast(pclump)
type(clump_type), pointer :: pclump
integer :: bonesite(3,100), i, j, npit
logical :: inlist
type(occupancy_type), pointer :: p
type(osteoclast_type), pointer :: pclast
real :: tnow
integer :: kpar = 0

tnow = istep*DELTA_T
nclast = nclast + 1
write(logmsg,*) 'make_osteoclast: ',nclast,tnow
call logger(logmsg)
pclast => clast(nclast)
pclast%ID = nclast
pclast%site = pclump%cm
!pclast%cm = pclump%cm
pclast%site(2) = NBY + 1
!pclast%cm(2) = NBY + 0.5
pclast%lastdir = random_int(1,8,kpar)
!pclast%normal = signal(isig)%normal
!pclast%fusetime = tnow
pclast%entrytime = tnow
pclast%status = ALIVE
pclast%movetime = tnow + CLAST_DWELL_TIME 
pclast%dietime = tnow + clastLifetime()
! Now need to create the pit list
npit = 0
do i = 1,pclump%ncells
	j = pclump%list(i)
	if (mono(j)%site(2) == NBY+1) then	! site on the bone surface
		npit = npit+1
		bonesite(:,npit) = mono(j)%site
	endif
	pclast%mono(i) = j
	mono(j)%iclast = nclast
	mono(j)%status = OSTEO
enddo
pclast%npit = npit
pclast%count = pclump%ncells
pclast%cm = 0
allocate(pclast%pit(npit))
do i = 1,npit
	pclast%pit(i)%site = bonesite(:,i)
	pclast%cm = pclast%cm + pclast%pit(i)%site
enddo
call pitrates(pclast)
pclast%cm = pclast%cm/npit
pclast%cm(2) = NBY + 0.5
pclump%status = DEAD
!write(logmsg,'(a,3i4,3f6.1,2i4)') 'clast site, cm, count, npit: ',pclast%site,pclast%cm,pclast%count,pclast%npit
!call logger(logmsg)
end subroutine

!------------------------------------------------------------------------------------------------
! The bone resorption rate at a given (x,z) depends on:
!	n = the number of monocytes that fused to make the osteoclast
!	d = the distance of the target bone site from the osteoclast centre
! The depth factor df decreases linearly to zero as d goes from 0 to MAX_RESORPTION_D
! The nominal maximum rate, MAX_RESORPTION_RATE, is scaled by the depth factor and
! by the monocyte count (relative to a nominal max MAX_RESORPTION_N).
! MAX_RESORPTION_RATE is in um/min, it is converted into grids/min by /DELTA_X
!------------------------------------------------------------------------------------------------
real function resorptionRate(n,d)
integer :: n
real :: d
real :: df

if (d >= MAX_RESORPTION_D) then
	df = 0
else
	df = 1 - d/MAX_RESORPTION_D
endif
resorptionRate = MAX_RESORPTION_RATE*(real(n)/MAX_RESORPTION_N)*df
end function

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine clastDeath(i)
integer :: i
integer :: k, site(3), kcell
integer :: x,y,z

!do k = 1,clast(i)%npit
!	if (clast(i)%pit(k)%fraction > 0 .and. clast(i)%pit(k)%fraction < 0.5) then
!		site = clast(i)%pit(k)%site
!		occupancy(site(1),site(2),site(3))%region = PIT
!		occupancy(site(1),site(2),site(3))%indx = 0
!	endif
!enddo
write(logmsg,*) 'clast death: ',i
call logger(logmsg)
clast(i)%status = DEAD
deallocate(clast(i)%pit)
do k = 1,clast(i)%count
	kcell = clast(i)%mono(k)
	site = mono(kcell)%site
!	write(logmsg,*) 'Mono dies: ',kcell,site
!	call logger(logmsg)
	mono(kcell)%status = DEAD
	occupancy(site(1),site(2),site(3))%indx = 0
enddo
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
real function clastLifetime()

clastLifetime = CLAST_LIFETIME
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
real :: tnow, t1, t2, fraction, ypit
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

open(nfpos,file=fname,status='new')
! Bone section
write(nfpos,'(a2,4i6)') 'B ',NX,NY,NZ,NBY
! Monocyte section
if (nmono > 0) then
	nfused = 0
    do kcell = 1,nmono
		status = mono(kcell)%status
		if (status == DEAD .or. status == LEFT) cycle
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
        if (.not.FAST_DISPLAY .or. mono_state /= 0) then
!		if (status == FUSED) then
	        write(nfpos,'(a2,i6,4i4)') 'M ',kcell-1, site, mono_state
	    endif
!        if (status == FUSING .or. status == FUSED) then
!			write(nflog,*) istep, kcell, site, mono_state
!		endif
    enddo
!    if (nfused > 0) then
!	    write(logmsg,*) 'nfused: ',nfused
!		call logger(logmsg)
!	endif
endif

! Capillary section
do icap = 1,ncap
	write(nfpos,'(a,6f6.1,f5.2)') 'C ',capillary(icap)%pos1,capillary(icap)%pos2,capillary(icap)%radius + 0.25
enddo

! Pit section
do x = 1,NX
	do z = 1,NZ
		do y = 1,NBY
			if (occupancy(x,y,z)%region == PIT .or. &
			   (occupancy(x,y,z)%region == BONE .and. occupancy(x,y,z)%bone_fraction < 1.0)) then
			   ypit = y - 0.5 + occupancy(x,y,z)%bone_fraction
				write(nfpos,'(a,3i4,f7.3)') 'P ',x,y,z,ypit
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

if (.not.clear_to_send) then
	! wait until the file called removefile exists, then remove it
	inquire(file=removefile,exist=ex)
	if (.not.ex) then
		call logger('wait')
		do
!			call logger('wait')
			inquire(file=removefile,exist=ex)
			if (.not.ex) then
!				call sleeper(1)
!				call millisleep(10) ! no good at all
			else
				exit
			endif
		enddo
	endif
	call unlink(removefile)
	clear_to_send = .true.
endif
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(NX_dim,NY_dim,NZ_dim,NBY_dim) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: NX_dim,NY_dim,NZ_dim,NBY_dim

NX_dim = NX
NY_dim = NY
NZ_dim = NZ
NBY_dim = NBY
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_summary(summaryData) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*)

summaryData(1:4) = (/istep,mono_cnt,nborn,nleft/)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_scene(ncap_list,cap_list,nmono_list,mono_list,npit_list,pit_list,nclast_list,clast_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene
use, intrinsic :: iso_c_binding
integer(c_int) :: nmono_list, mono_list(*), ncap_list, npit_list, nclast_list
real(c_float) :: cap_list(*), pit_list(*), clast_list(*)
integer :: k, j, kcell, iclump, site(3), nfused, x, y, z, status, mono_state, icap, iclast
real :: tnow, t1, t2, fraction, ypit, size, lastjump(3)
type(osteoclast_type), pointer :: pclast
real :: clast_diam = 0.9
real :: mono_diam = 0.5

tnow = istep*DELTA_T
! Monocyte section
if (nmono > 0) then
	nfused = 0
	k = 0
    do kcell = 1,nmono
		status = mono(kcell)%status
		if (status == DEAD .or. status == LEFT) cycle
		if (mono(kcell)%region /= MARROW) cycle
		if (status == CROSSING) then
			mono_state = 1
		elseif (status == CLUMPED) then
			mono_state = 2
		elseif (status == FUSING) then
			iclump = mono(kcell)%iclump
			t2 = clump(iclump)%fusetime
			t1 = clump(iclump)%starttime
			fraction = min(1.0,(tnow-t1)/(t2-t1))
			mono_state = 2 + fraction*98
		elseif (status == FUSED) then
			mono_state = 100
			nfused = nfused + 1
		elseif (status == OSTEO) then
			cycle
		else
            mono_state = 0
		endif
        site = mono(kcell)%site
        if (FAST_DISPLAY .and. mono_state < 2) cycle
		k = k+1
		j = 5*(k-1)
		mono_list(j+1) = kcell-1
		mono_list(j+2:j+4) = site
		mono_list(j+5) = mono_state
!		write(logmsg,*) k,mono_state
!		call logger(logmsg)
    enddo
endif
nmono_list = k

! Capillary section
!do icap = 1,ncap
!	write(nfpos,'(a,6f6.1,f5.2)') 'C ',capillary(icap)%pos1,capillary(icap)%pos2,capillary(icap)%radius + 0.25
!enddo
do icap = 1,ncap
	j = (icap-1)*7
	cap_list(j+1:j+3) = capillary(icap)%pos1
	cap_list(j+4:j+6) = capillary(icap)%pos2
	cap_list(j+7) = capillary(icap)%radius+0.25
enddo
ncap_list = ncap

! Pit section
k = 0
do x = 1,NX
	do z = 1,NZ
		do y = 1,NBY
			if (occupancy(x,y,z)%region == PIT .or. &
			   (occupancy(x,y,z)%region == BONE .and. occupancy(x,y,z)%bone_fraction < 1.0)) then
			   ypit = y - 0.5 + occupancy(x,y,z)%bone_fraction
			   k = k+1
			   j = 4*(k-1)
			   pit_list(j+1:j+3) = (/x,y,z/)
			   pit_list(j+4) = ypit
!				write(nfpos,'(a,3i4,f7.3)') 'P ',x,y,z,ypit
				exit
			endif
		enddo
	enddo
enddo
npit_list = k

! Osteoclast section
k = 0
do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status == DEAD) cycle
	k = k+1
	size = clast_size(pclast)
	j = 7*(k-1)
	clast_list(j+1:j+3) = pclast%cm
!	clast_list(j+2) = clast_list(j+2) - 0.5
!	clast_list(j+4:j+6) = (/1,0,0/)		! direction unit vector
	lastjump = dir2D(:,pclast%lastdir)
	clast_list(j+4:j+6) = lastjump/sqrt(dot_product(lastjump,lastjump))
	clast_list(j+7) = size
enddo
nclast_list = k
!write(logmsg,'(a,4i6)') '# of mono, cap, pit, clast: ',nmono_list, ncap_list, npit_list, nclast_list
!call logger(logmsg)
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

!read(nfinp,*) TC_AVIDITY_MEAN				! mean of avidity distribution (only if fix_avidity = false)
!read(nfinp,*) TC_AVIDITY_SHAPE			    ! shape -> 1 gives normal dist with small variance
!read(nfinp,*) TC_STIM_RATE_CONSTANT			! rate const for TCR stimulation (-> molecules/min)
!read(nfinp,*) TC_STIM_HALFLIFE				! halflife of T cell stimulation (hours)
!read(nfinp,*) divide_mean1
!read(nfinp,*) divide_shape1
read(nfinp,*) BETA							! speed: 0 < beta < 1		(0.65)
read(nfinp,*) RHO							! persistence: 0 < rho < 1	(0.95)

!read(nfinp,*) DC_LIFETIME_MEAN				! days
!read(nfinp,*) DC_LIFETIME_SHAPE 			! days

!read(nfinp,*) IL2_THRESHOLD					! stimulation needed to initiate IL-2/CD25 production
!read(nfinp,*) ACTIVATION_THRESHOLD			! stimulation needed for activation
!read(nfinp,*) FIRST_DIVISION_THRESHOLD		! activation level needed for first division
!read(nfinp,*) DIVISION_THRESHOLD			! activation level needed for subsequent division
!read(nfinp,*) EXIT_THRESHOLD				! activation level below which exit is permitted
!read(nfinp,*) STIMULATION_LIMIT				! maximum activation level

read(nfinp,*) X_SIZE						! size of bone region (square)
read(nfinp,*) Y_SIZE						! thickness of slice

read(nfinp,*) CAPILLARY_DIAMETER			! capillary diameter (um) (= 3)
read(nfinp,*) MONOCYTE_DIAMETER				! monocyte diameter (um) = 10	! um
read(nfinp,*) MONO_PER_MM3					! initial (equil) number of monocytes/mm3 (= 2000)
read(nfinp,*) IN_PER_HOUR					! rate of influx of monocytes from the blood
read(nfinp,*) STEM_PER_MM2					! number of stem cells/mm3 (=20) (NEEDS to be normalized)
read(nfinp,*) STEM_CYCLETIME				! stem cell division cycle time (hours) (= 6*60	! 6 hours)
read(nfinp,*) CROSSING_TIME					! time taken for a monocyte to cross into the capillary (mins)

read(nfinp,*) FUSING_TIME					! time taken by monocytes fusing into an osteoclast	(120) (mins)
read(nfinp,*) CLAST_LIFETIME				! lifetime of an osteoclast (4) days -> mins
read(nfinp,*) CLAST_DWELL_TIME				! time an osteoclast spends in one spot (180) (mins)
read(nfinp,*) MAX_RESORPTION_RATE			! maximum bone removal rate (/grid cell) (0.02) (um/min)
read(nfinp,*) MAX_RESORPTION_D				! maximum pit depth (for scaling rate) (5) (um)
read(nfinp,*) MAX_RESORPTION_N				! number of monos in osteoclast corresponding to MAX_RESORPTION_RATE (30)

read(nfinp,*) SIGNAL_RADIUS					! radius of influence of bone signal (um -> grids) (10)
read(nfinp,*) SIGNAL_THRESHOLD				! defines the high-signal region, near the source (0.14)
read(nfinp,*) SIGNAL_AFACTOR				! field amplification factor (0.4)
read(nfinp,*) MTHRESHOLD					! number of monocytes in the high-signal region that triggers fusing (25)

!read(nfinp,*) exit_rule						! 1 = no chemotaxis, 2 = chemotaxis
!read(nfinp,*) exit_region					! region for cell exits 1 = capillary, 2 = sinusoid
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

DELTA_X = MONOCYTE_DIAMETER
NX = X_SIZE/DELTA_X									! convert um -> grids
NY = Y_SIZE/DELTA_X	+NBY							! convert um -> grids
if (mod(NX,2) /= 0) NX = NX+1						! ensure that NX is even (why?)
NZ = NX
CAPILLARY_DIAMETER = CAPILLARY_DIAMETER/DELTA_X		! convert um -> grids
chemo_radius = chemo_radius/DELTA_X					! convert um -> grids
STEM_CYCLETIME = 60*STEM_CYCLETIME					! convert hours -> minutes
CLAST_LIFETIME = CLAST_LIFETIME*24*60				! convert days -> minutes
MAX_RESORPTION_RATE = MAX_RESORPTION_RATE/DELTA_X	! convert um/min -> grids
MAX_RESORPTION_D = MAX_RESORPTION_D/DELTA_X			! convert um -> grids/min
SIGNAL_RADIUS = SIGNAL_RADIUS/DELTA_X				! convert um -> grids

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

write(logmsg,*) 'FUSING_TIME: ',FUSING_TIME					! time taken by monocytes fusing into an osteoclast	(120) (mins)
call logger(logmsg)
write(logmsg,*) 'CLAST_LIFETIME: ',CLAST_LIFETIME				! lifetime of an osteoclast (4) days -> mins
call logger(logmsg)
write(logmsg,*) 'CLAST_DWELL_TIME: ',CLAST_DWELL_TIME				! time an osteoclast spends in one spot (180) (mins)
call logger(logmsg)
write(logmsg,*) 'MAX_RESORPTION_RATE: ',MAX_RESORPTION_RATE			! maximum bone removal rate (/grid cell) (0.02) (um/min)
call logger(logmsg)
write(logmsg,*) 'MAX_RESORPTION_D: ',MAX_RESORPTION_D				! maximum pit depth (for scaling rate) (5) (um)
call logger(logmsg)
write(logmsg,*) 'MAX_RESORPTION_N: ',MAX_RESORPTION_N				! number of monos in osteoclast corresponding to MAX_RESORPTION_RATE (30)
call logger(logmsg)

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

if (use_CPORT1) then
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
endif
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
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step
use, intrinsic :: iso_c_binding
integer(c_int) :: res
logical :: ok
integer :: error

res = 0
!call logger("simulate_step")
ok = .true.
istep = istep + 1
if (mod(istep,1000) == 0) then
	write(logmsg,'(a,i8,6i6)') 'istep: ',istep,mono_cnt,nleft	!,mono(22)%status,mono(22)%site
	call logger(logmsg)
endif
!write(nflog,*) 'call updater'
call updater
!write(nflog,*) 'call mover: ',nmono
call mover
!if (nsignal > 0) then
!	call checkSignals(ok)
!	if (.not.ok) res = 1
!endif
end subroutine


!------------------------------------------------------------------------------------------------
! This is the original version, using files to communicate signals and data with the Qt main.
!------------------------------------------------------------------------------------------------
subroutine simulate(ok)
logical :: ok
integer :: it, error, res

istep = 0
ok = .true.
do it = 1,nsteps
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
	call simulate_step(res)
	if (res /= 0) then
		ok = .false.
		return
	endif
!	call updater
!	call mover
!	if (nsignal > 0) then
!		call checkSignals(ok)
!		if (.not.ok) return
!	endif
enddo
end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine report
integer :: i, k, indx, status

write(logmsg,*) 'nclump: ',nclump
call logger(logmsg)
do i = 1,nclump
	if (clump(i)%status < 0) cycle
	write(logmsg,*) 'clump: ',i,'  status: ',clump(i)%status,'  ncells: ',clump(i)%ncells
	call logger(logmsg)
	do k = 1,clump(i)%ncells
		indx = clump(i)%list(k)
		write(logmsg,*) indx,mono(indx)%site
		call logger(logmsg)
	enddo
enddo
!write(logmsg,*) 'Clumped monocytes  '
!call logger(logmsg)
!do i = 1,nmono
!	status = mono(i)%status
!	if (status == LEFT .or. status == DEAD) cycle
!	if (mono(i)%iclump > 0) then
!		write(logmsg,*) i,mono(i)%iclump,mono(i)%site
!		call logger(logmsg)
!	endif
!enddo
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
if (allocated(entrysite)) deallocate(entrysite)
if (allocated(RANKL_conc)) deallocate(RANKL_conc)
if (allocated(S1P_conc)) deallocate(S1P_conc)
if (allocated(S1P_grad)) deallocate(S1P_grad)

end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run
use, intrinsic :: iso_c_binding
integer(c_int) :: res
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

!call report
call wrapup

if (res == 0) then
	call logger(' Execution successful')
else
	call logger('  === Execution failed ===')
	call sleeper(1)
endif
close(nflog)

if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    if (use_CPORT1) call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
!	    call logger("closed PORT_0")
		if (use_CPORT1) then
		    call winsock_send(awp_1,quit,8,error)
		    call winsock_close(awp_1)
!			call logger("closed PORT_1")
		endif
	endif
endif

end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine execute(infile_array,buflen) BIND(C)
!!DEC$ ATTRIBUTES DLLEXPORT :: EXECUTE
!!DEC$ ATTRIBUTES C, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"execute" :: execute
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128)
integer(c_int) :: buflen
character*(128) :: infile
logical :: ok
integer :: i, res

use_CPORT1 = .false.	! TESTING DIRECT CALLING FROM C++
infile = ''
do i = 1,buflen
	infile(i:i) = infile_array(i)
enddo

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
call logger('did setup')
if (ok) then
	call logger('returned OK')
	clear_to_send = .true.
	simulation_start = .true.
	istep = 0
	
	return
	
	call simulate(ok)
	call logger('Ended simulation')
	call logger('Execution successful')
else
	call logger('setup failed')
endif
if (ok) then
	res = 0
else
	res = 1
endif
call terminate_run(res)
end subroutine

end module

