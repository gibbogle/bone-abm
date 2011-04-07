module motion
use global
use fields
implicit none
save

integer :: reldir(6,26)

integer :: MODEL = MOORE26_MODEL
real, parameter :: Kchemo = 1.0

integer :: njumpdirs, nreldir
integer :: jumpvec(3,MAXRELDIR+1)
real :: unitjump(3,MAXRELDIR+1)
real :: dirprob(0:MAXRELDIR)
logical :: vn_adjacent(MAXRELDIR+1)
integer :: dir2D(3,8) = reshape((/ -1,0,-1, -1,0,0, -1,0,1, 0,0,1, 1,0,1, 1,0,0, 1,0,-1, 0,0,-1/),(/3,8/))


contains

subroutine move_clast1(pclast,res)
type(osteoclast_type), pointer :: pclast
integer :: res
integer :: ipit, iy, x, y, z, site(3), i, imono, kdir, dx, dz, dirmax, iclast, ddir, nsig
real :: prob(0:4) = (/ 0.5, 0.15, 0.07, 0.025, 0.01 /)
integer :: jump(3), lastjump(3), kpar=0
real :: d, bf, v(3), proj, size0, size1, djump, dsum, depth, totsig(0:8), siglim
real(8) :: psum, pmax, R, dp, p(8), ptemp(8)
!real :: p(3) = (/0.25,0.5,0.25/)	! arbitrary, interim
logical :: covered, bdryhit, nearclast, freshbone, possible(8)
type(osteoclast_type), pointer :: pclast1
integer, save :: count = 0
logical :: dbug

if (pclast%ID == -1) then
	dbug = .true.
else
	dbug = .false.
endif

possible = .true.
freshbone = .false.
res = 0
! First choose a direction.  Quantify the new bone available in each direction
! at a distance approx equal to the long axis dimension of the osteoclast, i.e.
! sqrt(count).

site = pclast%cm + 0.5
totsig(0) = TotalSignal(pclast,site)
nsig = 0
size0 = pclast%radius
if (dbug) write(*,'(a,i4,4f8.1)') 'move_clast: ',pclast%ID,pclast%cm,size0
pmax = 0
do kdir = 1,8
	p(kdir) = 0
	jump = dir2D(:,kdir)
!	djump = sqrt(real(dot_product(jump,jump)))
!	dsum = 0
!	do i = 1,20
!		dsum = dsum + djump
!		if (dsum > size0) exit
!	enddo
!	site = pclast%cm + i*jump + 0.5
	site = pclast%cm + jump + 0.5
	totsig(kdir) = TotalSignal(pclast,site)
	if (totsig(kdir) > 0) nsig = nsig + 1
	if (dbug) write(*,*) 'kdir,i,site: ',kdir,i,site
	bdryhit = .false.
	if ((site(1) <= 1 .or. site(1) >= NX) .or. (site(3) <= 1 .or. site(3) >= NZ)) then
		bdryhit = .true.
		possible(kdir) = .false.
		write(*,*) 'bdryhit: ',kdir
		cycle
	endif
!	if (occupancy(site(1),site(2),site(3))%indx /= 0) then
!		possible(kdir) = .false.
!		write(*,*) 'occupancy: ',kdir
!		cycle
!	endif
	! Check for nearby osteoclasts.
	! Treat osteoclast footprint as a circle with radius = clast_size
	nearclast = .false.
	do iclast = 1,nclast
		pclast1 => clast(iclast)
		if (associated(pclast,pclast1)) cycle
		if (pclast1%status == DEAD) cycle
		size1 = pclast1%radius
		v = pclast1%cm - pclast%cm
		d = sqrt(dot_product(v,v))
		if (d < (size0 + size1)) then	! we don't want to move in the direction of v
			proj = dot_product(v,real(jump))
			if (proj > 0.5) then
!				write(*,*) 'Too close: ', iclast
				nearclast = .true.
				exit
!				possible(kdir) = .false.
			endif
		endif
	enddo	
	if (nearclast) then
		possible(kdir) = .false.
		write(*,*) 'nearclast: ',kdir
		cycle
	endif
!	x = site(1)
!	z = site(3)
!	do y = NBY,1,-1
!		bf = occupancy(x,y,z)%bone_fraction
!		if (bf > 0) then
!			depth = (NBY + 0.5) - (y - 0.5 + bf)
!			exit
!		endif
!	enddo

! DO NOT USE DEPTH OR PERSISTENCE OF DIRECTION NOW, USE SIGNAL
!	depth = surface(x,z)%depth
!	p(kdir) = exp(-Kattraction*depth)
!	if (p(kdir) > 0.9) freshbone = .true.	
	! At this stage p(:) contains the relative attractiveness of this direction,
	! without accounting for the direction of the previous jump.
	! Now multiply by the weight that represents persistence of direction
!	ddir = abs(pclast%lastdir - kdir)
!	if (ddir > 4) then
!		ddir = abs(ddir-8)
!	endif
!	if (dbug) write(*,'(5i3,5f8.3)') kdir,ddir,site,depth,p(kdir),prob(ddir)
!	p(kdir) = p(kdir)*prob(ddir)
!	pmax = max(pmax,p(kdir))
enddo

!if (dbug) write(*,'(a,8f7.4)') 'p: ',p

! An OC can move to a new site if the move is possible and if the total signal at the new site
! exceeds (sufficiently) that at the current site, or if the current total signal = 0.  
! The probability of the move is proportional to the signal excess.
! If the current total signal and all neighbour site total signals are zero, the osteoclast dies. 
siglim = min(totsig(0)*1.3,1.0)
pmax = 0
do kdir = 1,8
	if (possible(kdir) .and. totsig(kdir) > siglim) then
		p(kdir) = totsig(kdir)
		pmax = max(pmax,p(kdir))
	endif
enddo

! Now filter out direction probs that are too small in comparison with the prob of the most preferred direction
do kdir = 1,8
	if (p(kdir) < 0.2*pmax) then
		p(kdir) = 0
	endif
enddo
psum = sum(p)
if (psum > 0) then
	p = p/psum
endif
if (pclast%ID == 1) then
	write(*,'(a,9f8.1)') 'totsig: ',totsig
	write(*,'(a,8f8.4)') 'p: ',p
endif
if (psum == 0) then
	if (totsig(0) == 0 .and. nsig == 0) then
		res = 2
	else
		res = 0
	endif
	return
endif
p = p/psum

!write(*,'(9f7.3)') pmax,p
!if (.not.freshbone) then	! no good moves, need to find some fresh bone
!!	write(logmsg,*) 'No good moves for this osteoclast: ', pclast%ID
!!	call logger(logmsg)
!	if (possible(pclast%lastdir)) then	! keep going in the same direction
!		kdir = pclast%lastdir
!		res = 1
!	else								! choose a new direction
!		p = 0
!		do i = 1,8
!			if (possible(i)) p(i) = 1
!		enddo
!		if (sum(p) == 0) then
!			res = 2			! die! (maybe)
!			return
!		endif
!		res = -1
!	endif
!endif

!if (res /= 1) then
	R = par_uni(kpar)
	psum = 0
	do kdir = 1,8
		psum = psum + p(kdir)
		if (psum > R) exit
	enddo
	kdir = min(8,kdir)
!endif
jump = dir2D(:,kdir)
write(*,*) 'clast moves: kdir: ',pclast%ID,kdir,jump
!do ipit = 1,pclast%npit
!	x = pclast%pit(ipit)%site(1) + jump(1)
!	y = pclast%pit(ipit)%site(2)
!	z = pclast%pit(ipit)%site(3) + jump(3)
!	y = 0
!	do iy = NBY,1,-1
!		if (occupancy(x,iy,z)%bone_fraction > 0) then
!			y = iy
!			exit
!		endif
!	enddo
!	if (y == 0) then
!		write(logmsg,*) 'Error: moveclast: y = 0: ',ipit,x,y,z
!		call logger(logmsg)
!		stop
!	endif
!	pclast%pit(ipit)%site = (/x,y,z/)
!enddo
!pclast%site = pclast%site + jump
pclast%cm = pclast%cm + jump
pclast%lastdir = kdir
!do i = 1,pclast%count
!	imono = pclast%mono(i)
!	site = mono(imono)%site
!	occupancy(site(1),site(2),site(3))%indx = 0
!	site = site + jump
!	mono(imono)%site = site
!	occupancy(site(1),site(2),site(3))%indx = imono
!enddo
res = abs(res)
end subroutine


!---------------------------------------------------------------------
!---------------------------------------------------------------------
real function TotalSignal(pclast,site)
type(osteoclast_type), pointer :: pclast
integer :: site(3)
integer :: k, x, z

TotalSignal = 0
do k = 1,pclast%npit
	x = site(1) + pclast%pit(k)%delta(1)
	z = site(3) + pclast%pit(k)%delta(3)
	TotalSignal = TotalSignal + surface(x,z)%signal
enddo
end function

!---------------------------------------------------------------------
! The signal in the OC foorprint is summed, but set to -1 if it
! overlaps another OC footprint.
!---------------------------------------------------------------------
real function ClearSignal(pclast,site)
type(osteoclast_type), pointer :: pclast
integer :: site(3)
integer :: k, x, z, iclast, k1, x1, z1, site1(3)

ClearSignal = 0
do k = 1,pclast%npit
	x = site(1) + pclast%pit(k)%delta(1)
	z = site(3) + pclast%pit(k)%delta(3)
	do iclast = 1,nclast
		if (clast(iclast)%ID == pclast%ID) cycle
		if (clast(iclast)%status /= RESORBING) cycle
		site1 = clast(iclast)%site		! site of another resorbing OC
		do k1 = 1,clast(iclast)%npit
			x1 = site1(1) + clast(iclast)%pit(k1)%delta(1)
			z1 = site1(3) + clast(iclast)%pit(k1)%delta(3)
			if (x == x1 .and. z == z1) then	! footprints overlap
				ClearSignal = -1
				return
			endif
		enddo
	enddo
	ClearSignal = ClearSignal + surface(x,z)%signal
enddo
end function

!---------------------------------------------------------------------
! The osteoclast is moved one lattice jump, and the locations of active
! pit sites are recomputed, together with their resorption rates.
! Chemotaxis
! When osteoclast motion is influenced by the bone signal, the jump 
! probabilities depend on the relative signal strengths of the 
! neighbour sites.
! There are two cases:
! (1) OC is on or near sites needing excavation (signal > 0)
! (2) OC is not near sites needing excavation.  In this case the OC 
!     must first move towards sites with signal.  It's probably safe to
!     assume that the OC has a range for detection of signal.
!
! Result values (res):
!     0 = prefer to stay, continue resorbing
!     1 = prefer or need to move, move possible
!     2 = move necessary (no signal), move blocked
!     3 = no signal in patch
!---------------------------------------------------------------------
subroutine move_clast(pclast,res)
type(osteoclast_type), pointer :: pclast
integer :: res
integer :: ipit, iy, x, y, z, site(3), i, imono, kdir, dx, dz, dirmax, iclast
integer :: ddir, nsig, ocsite(3), kmax, targetsite(3), kmin
real :: prob(0:4) = (/ 0.5, 0.15, 0.07, 0.025, 0.01 /)
integer :: jump(3), lastjump(3), kpar=0
real :: d, bf, v(3), proj, size0, size1, djump, dsum, depth
real :: totsig(0:8), siglim, sig, patchtot, amp(8), cosa, amax, tnow, dmin
real(8) :: psum, pmax, R, dp, p(8), ptemp(8)
!real :: p(3) = (/0.25,0.5,0.25/)	! arbitrary, interim
logical :: covered, bdryhit, nearclast, freshbone, possible(8)
type(osteoclast_type), pointer :: pclast1
integer, save :: count = 0
logical :: dbug
real, parameter :: SIGNAL_EXCESS = 1.05

if (pclast%ID == -1) then
	dbug = .true.
else
	dbug = .false.
endif

tnow = istep*DELTA_T
possible = .true.
freshbone = .false.
res = 0
! First choose a direction.  Quantify the new bone available in each direction
! at a distance approx equal to the long axis dimension of the osteoclast, i.e.
! sqrt(count).

ocsite = pclast%site
if (pclast%status == MOVING) then	! special case, OC moving fast, possibly over other cells
	targetsite = pclast%targetsite
	if (ocsite(1) == targetsite(1) .and. ocsite(3) == targetsite(3)) then	! OC has reached the target site
		pclast%status = RESORBING
		res = 0
		return
	endif
	dmin = 1.0e10
	do kdir = 1,8
		jump = dir2D(:,kdir)
		site = ocsite + jump
		v = targetsite - site
		v(2) = 0
		d = sqrt(dot_product(v,v))
		if (d < dmin) then
			dmin = d
			kmin = kdir
		endif
	enddo
	call ocmove(pclast,kmin)
	if (dmin == 0) then
		pclast%status = RESORBING
	endif
	res = 1
	return
endif

totsig(0) = TotalSignal(pclast,ocsite)
nsig = 0
size0 = pclast%radius
!if (dbug) write(*,'(a,i4,4f8.1)') 'move_clast: ',pclast%ID,pclast%cm,size0
do kdir = 1,8
	jump = dir2D(:,kdir)
	site = ocsite + jump
	totsig(kdir) = TotalSignal(pclast,site)
	if (totsig(kdir) > 0) nsig = nsig + 1
!	if (dbug) write(*,*) 'kdir,i,site: ',kdir,i,site
	bdryhit = .false.
	if ((site(1) <= 1 .or. site(1) >= NX) .or. (site(3) <= 1 .or. site(3) >= NZ)) then
		bdryhit = .true.
		possible(kdir) = .false.
		cycle
	endif
	! Check for nearby osteoclasts.
	! Treat osteoclast footprint as a circle with radius = clast_size
	nearclast = .false.
	do iclast = 1,nclast
		pclast1 => clast(iclast)
		if (pclast%ID == pclast1%ID) cycle
		if (pclast1%status == DEAD) cycle
		size1 = pclast1%radius
		v = pclast1%site - pclast%site
		d = sqrt(dot_product(v,v))
		if (d < (size0 + size1)) then	! we don't want to move in the direction of v
			proj = dot_product(v,real(jump))
			if (proj > 0) then
!				write(*,*) 'Too close: ', iclast
				nearclast = .true.
				exit
!				possible(kdir) = .false.
			endif
		endif
	enddo	
	if (nearclast) then
		if (dbug) write(*,*) 'nearclast: ',kdir,iclast,size0,size1
		possible(kdir) = .false.
		cycle
	endif
enddo

if (dbug) write(*,'(9f8.2)') totsig
if (sum(totsig) == 0) then	
	! OC is not near sites needing excavation.
	! Try to move towards a signalling site within range.
	! If no such movement is possible, turn off resorption
!	write(*,*) 'sum(totsig): ',sum(totsig)
	patchtot = 0
	amp = 0
	dmin = 999
	do i = 1,nsignal
		site = signal(i)%site
		sig = surface(site(1),site(3))%signal
		patchtot = patchtot + sig
		if (sig > 0) then
			v = site - ocsite
			d = sqrt(dot_product(v,v))
			dmin = min(d,dmin)
			if (d <= OC_SIGNAL_SENSING_RANGE) then
				v = v/d		! unit offset vector
				! Look at jump directions that reduce distance to the signalling site
				do kdir = 1,8
!					if (.not.possible(kdir)) cycle
					cosa = dot_product(v,dir2D(:,kdir))
					if (stuck) then
						write(logmsg,'(3i5,L2,2f8.4)') pclast%ID,i,kdir,possible(kdir),sig,cosa
						call logger(logmsg)
					endif
					if (possible(kdir) .and. cosa > 0) then
						amp(kdir) = amp(kdir) + cosa*sig/d
					endif
				enddo
			endif
		endif
	enddo
	if (dmin < OC_SIGNAL_SENSING_RANGE .and. sum(amp) == 0) then
		stuck = .true.
		write(logmsg,'(a,2f8.3)') 'stuck: patchtot,dmin: ',patchtot,dmin
		call logger(logmsg)
		write(logmsg,'(a,8f8.4)') 'amp: ',amp
		call logger(logmsg)
	endif
	if (patchtot == 0) then		! Excavation is complete
		res = 3
		return
	endif
	amax = 0
	kmax = 0
	do kdir = 1,8
		if (amp(kdir) > amax) then
			amax = amp(kdir)
			kmax = kdir
		endif
	enddo
	if (kmax == 0) then		! No useful movement is currently possible.  Need to climb over.
		call ChooseTarget(pclast,targetsite)
		if (targetsite(1) == 0) then
			res = 2
			return
		else
			pclast%targetsite = targetsite
			pclast%status = MOVING
			pclast%movetime = tnow + DT_FAST_MOVE
			res = 0
			return
		endif
	endif
	call ocmove(pclast,kmax)
	res = 1
	return
endif

! An OC can move to a new site if the move is possible and if the total signal at the new site
! exceeds (sufficiently) that at the current site, or if the current total signal = 0.  
! The probability of the move is proportional to the signal excess.
! If the current total signal and all neighbour site total signals are zero, the osteoclast dies. 

if (totsig(0) == 0) then	! any level of signal is attractive
	siglim = 0
else						! need sufficient signal to force a move
	siglim = max(totsig(0)*SIGNAL_EXCESS,1.0)
endif
if (dbug) write(*,*) 'siglim: ',siglim,' possible: ',possible
p = 0
pmax = 0
do kdir = 1,8
	if (possible(kdir) .and. (totsig(kdir) > siglim)) then
		p(kdir) = totsig(kdir)
		if (p(kdir) > pmax) then
			pmax = p(kdir)
			kmax = kdir
		endif
	endif
enddo

if (totsig(0) == 0) then	! need to move
	if (pmax == 0) then		! no useful move is possible, try a random move
		do kdir = 1,8
			if (possible(kdir)) then
				R = par_uni(kpar)
				if (R < 0.1) then
					call ocmove(pclast,kdir)
					res = 1
					return
				endif
			endif
		enddo
		call ChooseTarget(pclast,targetsite)
		if (targetsite(1) == 0) then
			res = 2
			return
		else
			pclast%targetsite = targetsite
			pclast%status = MOVING
			pclast%movetime = tnow + DT_FAST_MOVE
			res = 0
			return		
		endif
	else					! make a move
		call ocmove(pclast,kmax)
		res = 1
		return
	endif
elseif (pmax == 0) then		! do not move
	res = 0
	return
endif	

! Now filter out direction probs that are too small in comparison with the prob of the most preferred direction
do kdir = 1,8
	if (p(kdir) < 0.2*pmax) then
		p(kdir) = 0
	endif
enddo
psum = sum(p)
if (psum > 0) then
	p = p/psum
endif

R = par_uni(kpar)
psum = 0
do kdir = 1,8
	psum = psum + p(kdir)
	if (psum >= R) exit
enddo
call ocmove(pclast,kdir)
res = 1
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine ChooseTarget(pclast,targetsite)
type(osteoclast_type), pointer :: pclast
integer :: targetsite(3)
integer :: ocsite0(3), ocsite1(3), x, z
real :: v(3), d, sig, amp, amax

ocsite0 = pclast%site
targetsite = 0
amax = 0
do x = 1,NX
	do z = 1,NZ
		if (surface(x,z)%target_depth > 0) then	! consider this as a possible OC target site
			ocsite1 = (/ x, NBY+1, z /)
			v = ocsite1 - ocsite0
			v(2) = 0
			sig = ClearSignal(pclast,ocsite1)
			if (sig <= 0) cycle
			d = sqrt(dot_product(v,v))
			if (d <= OC_SIGNAL_SENSING_RANGE) then
				amp = sig/d
				if (amp > amax) then
					amax = amp
					targetsite = ocsite1
				endif
			endif
		endif
	enddo
enddo

end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine ocmove(pclast,kdir)
type(osteoclast_type), pointer :: pclast
integer :: kdir
integer :: jump(3), site(3), i, imono

jump = dir2D(:,kdir)
if (pclast%ID == 1) then
	write(logmsg,'(a,5i4)') 'clast moves: kdir: ',pclast%ID,kdir,jump
	call logger(logmsg)
endif
pclast%site = pclast%site + jump
pclast%lastdir = kdir
!do i = 1,pclast%count
!	imono = pclast%mono(i)
!	site = mono(imono)%site
!	occupancy(site(1),site(2),site(3))%indx = 0
!	site = site + jump
!	mono(imono)%site = site
!	occupancy(site(1),site(2),site(3))%indx = imono
!enddo
call UpdateSurface
end subroutine

!------------------------------------------------------------------------------------------------
! Refresh the locations of OCs on the surface.
!------------------------------------------------------------------------------------------------
subroutine UpdateSurface
type(osteoclast_type), pointer :: pclast
integer :: iclast, i, site(3), x, z

surface%iclast = 0
do iclast = 1,nclast
	pclast => clast(iclast)
	if (pclast%status == DEAD) cycle 
	do i = 1,pclast%npit
		site = pclast%site + pclast%pit(i)%delta
		x = site(1)
		z = site(3)
		if (surface(x,z)%iclast == 0) then
			surface(x,z)%iclast = iclast
		elseif (surface(x,z)%iclast < 100) then
			surface(x,z)%iclast = 100*surface(x,z)%iclast + iclast
		endif
	enddo
enddo
end subroutine

!------------------------------------------------------------------------------------------------
! We need a way to determine whether a monocyte has moved close enough to the bone surface to
! make contact with an OB - needed to initiate clustering.
! Base the criterion on nearness to a signalling site.
!------------------------------------------------------------------------------------------------
logical function NearOB(pmono)
type(monocyte_type), pointer :: pmono
real :: v(3), d2, d2lim
integer :: msite(3), ssite(3), isig
real, parameter :: OB_REACH = 40	! microns

NearOB = .false.
d2lim = (OB_REACH/DELTA_X)**2
msite = pmono%site
do isig = 1,nsignal
	ssite = signal(isig)%site
	if (surface(ssite(1),ssite(3))%signal == 0) cycle
	v = msite - ssite
	d2 = dot_product(v,v)
	if (d2 < d2lim) then
		NearOB = .true.
		exit
	endif
enddo
end function

!---------------------------------------------------------------------
! Currently monocytes move around while in the marrow, and may pass
! into the blood, effectively removed from further consideration.
!---------------------------------------------------------------------
subroutine mover
type (monocyte_type), pointer :: pmono
integer :: kcell, status
logical :: go
integer :: kpar = 0
integer :: kdbug = -1
integer :: site(3)
real :: tnow

tnow = istep*DELTA_T
do kcell = 1,nmono
	pmono => mono(kcell)
	status = pmono%status
	if (status == LEFT .or. status == DEAD) cycle
!	if (pmono%site(2) == NBY+1 .and. status < FUSED) then
!		write(*,*) 'y = NBY+1: ',istep,pmono%ID,pmono%site,status
!		stop
!	endif
	if (status == CROSSING) then
		if (tnow >= pmono%exittime) then
			pmono%status = LEFT
			pmono%region = BLOOD
			site = pmono%site
			occupancy(site(1),site(2),site(3))%indx = 0
			mono_cnt = mono_cnt - 1
			nleft = nleft + 1
!			write(*,'(a,2i6,2f6.3,i6)') 'monocyte leaves: ',istep,kcell,mono(kcell)%S1P1,mono_cnt
		endif
	elseif (status >= MOTILE .and. pmono%iclump == 0) then	! interim criterion
		call mono_jumper(kcell,go,kpar)
		if (kcell == kdbug) then
			write(*,'(3i4)') pmono%site
		endif
	endif
	if (pmono%stickiness > 0 .and. pmono%status < FUSED .and. NearOB(pmono)) then
		call sticker(kcell)
	endif
enddo
	
end subroutine

!--------------------------------------------------------------------------------
! When a monocyte is within the field of influence of a signal site, the
! jump probabilities are modified to increase the prob of jumps towards
! the site.  When close enough to the signal site (intensity above a threshold)
! the jump probabilities are determined solely by the relative intensities.
! This is to ensure that monocytes cluster.
!--------------------------------------------------------------------------------
subroutine mono_jumper(kcell,go,kpar)
integer :: kpar,kcell
logical :: go
type (monocyte_type), pointer :: cell
integer :: site1(3),site2(3)
integer :: region, kcell2
integer :: irel,dir1,lastdir1,status
integer :: savesite2(3,26), jmpdir(26)
real(8) :: psum, p(26), R, wS1P(27), g(3), gamp, S1Pfactor, wRANKL(27), RANKLfactor
real :: f0, f, motility_factor
logical :: free, cross, field

cell => mono(kcell)
status = cell%status
if (status == MOTILE) then
	motility_factor = 1
elseif (status == CHEMOTACTIC) then
	motility_factor = 0.6
elseif (status == STICKY) then
	motility_factor = 0.4
endif
!if (istep - cell%lastmovestep > 1050) then
!	write(logmsg,*) 'No move: ',cell%ID,istep,kcell,istep - cell%lastmovestep,cell%site
!	call logger(logmsg)
!	stop
!elseif (istep - cell%lastmovestep > 1000) then
!	write(logmsg,*) 'No move: ',cell%ID,istep,kcell,istep - cell%lastmovestep,cell%site
!	call logger(logmsg)
!endif
site1 = cell%site
go = .false.
field = .false.
f0 = 0
!if (cell%status >= CHEMOTACTIC .and. occupancy(site1(1),site1(2),site1(3))%signal /= 0) then
!	field = .true.
!	f0 = occupancy(site1(1),site1(2),site1(3))%intensity
!endif

! TESTING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
field = .false.

R = par_uni(kpar)
!call random_number(R)
if (R <= dirprob(0)) then    ! case of no jump
	return
endif

! Now we must jump (if possible)

! Set up weights in the 6 principal axis directions corresponding to S1P_grad(:,:,:,:)
if (S1P_chemotaxis) then
	g = S1P_grad(:,site1(1),site1(2),site1(3))
	gamp = sqrt(dot_product(g,g))	! Magnitude of S1P gradient
	g = g/gamp
	call chemo_weights(g,wS1P)		! w(:) now holds the normalized gradient vector components
	S1Pfactor = S1P_CHEMOLEVEL*min(1.0,gamp/S1P_GRADLIM)*cell%S1P1
!	write(*,*) 'g,gamp: ',g,gamp
!	write(*,*) 'wS1P: ',wS1P
!	write(*,*) 'S1P1: ',cell%S1P1
endif
! Set up weights in the 6 principal axis directions corresponding to RANKL_grad(:,:,:,:)
if (RANKL_chemotaxis) then
	g = RANKL_grad(:,site1(1),site1(2),site1(3))
	gamp = sqrt(dot_product(g,g))	! Magnitude of RANKL gradient
	g = g/gamp
	call chemo_weights(g,wRANKL)		! w(:) now holds the normalized gradient vector components
	RANKLfactor = RANKL_CHEMOLEVEL*min(1.0,gamp/RANKL_GRADLIM)*cell%RANKSIGNAL
!	write(*,*) 'g,gamp: ',g,gamp
!	write(*,*) 'wRANKL: ',wRANKL
!	write(*,*) 'RANK: ',cell%RANKSIGNAL
endif
!write(*,*) 'S1Pfactor, RANKLfactor: ',S1Pfactor, RANKLfactor
!stop
lastdir1 = cell%lastdir
p = 0
psum = 0
do irel = 1,nreldir
    p(irel) = 0
	dir1 = reldir(lastdir1,irel)
	site2 = site1 + jumpvec(:,dir1)
	! With the call to free_site() returning region and kcell, decisions can be made
	! about transition to BLOOD, for example.
	free = free_site(site2,region,kcell2)
	if (free) then
		p(irel) = motility_factor*dirprob(irel)
		if (field) then
			! The probability of a jump is modified by the relative signal intensities of the two sites
			! This is a very crude interim treatment
!			f = occupancy(site2(1),site2(2),site2(3))%intensity
			if (f >= SIGNAL_THRESHOLD) then
				p(irel) = max(0.0,f-f0)
			else
				p(irel) = max(0.0,p(irel) + SIGNAL_AFACTOR*(f-f0))
			endif
		endif
		if (S1P_chemotaxis) then
			! the increment to the probability depends on S1Pfactor and wS1P(dir1)
			p(irel) = p(irel) + S1Pfactor*wS1P(dir1)
		endif
		if (RANKL_chemotaxis) then
			! the increment to the probability depends on RANKLfactor and wRANKL(dir1)
			p(irel) = p(irel) + RANKLfactor*wRANKL(dir1)
		endif
		jmpdir(irel) = dir1
		psum = psum + p(irel)
		savesite2(:,irel) = site2
	elseif (region == BLOOD .and. vn_adjacent(dir1)) then
		if (CrossToBlood(kcell,site1)) then
			return
		endif
!		write(*,*) 'Not free_site: ',site2
	endif
enddo

if (psum == 0) then
!	cell%lastdir = random_int(1,6,kpar)
!	if (field) then
!		write(*,*) 'clustering: ',kcell
!	endif
	write(logmsg,*) 'psum = 0: ',cell%ID
	call logger(logmsg)
	return
else
    go = .true.
endif

! Now choose a direction on the basis of these probs p()
R = par_uni(kpar)
!call random_number(R)
R = psum*R
psum = 0
do irel = 1,nreldir
   	psum = psum + p(irel)
   	if (R <= psum) then
   		exit
   	endif
enddo
if (irel > nreldir) then
	irel = nreldir
endif
site2 = savesite2(:,irel)
if (irel <= nreldir) then
	dir1 = reldir(lastdir1,irel)
	if (dir1 == 0) then
		write(logmsg,*) 'dir1 = 0: ',kcell,lastdir1,irel,dir1
		call logger(logmsg)
	endif
endif
if (diagonal_jumps) then
	dir1 = fix_lastdir(dir1,kpar)
elseif (dir1 == 0) then
	dir1 = random_int(1,6,kpar)
endif

cell%site = site2
cell%lastdir = dir1
cell%lastmovestep = istep
occupancy(site2(1),site2(2),site2(3))%indx = kcell
occupancy(site1(1),site1(2),site1(3))%indx = 0
end subroutine

!---------------------------------------------------------------------
! Use for both S1P and RANKL chemotaxis
!---------------------------------------------------------------------
subroutine chemo_weights(g,w)
real(8) :: g(3), w(27)

w = 0
if (g(1) < 0) then		! -x
	w(5) = -g(1)		! reldir26(1,1)
else					! +x
	w(23) = g(1)		! reldir26(2,1)
endif
if (g(2) < 0) then		! -y
	w(11) = -g(2)		! reldir26(3,1)
else					! +y
	w(17) = g(2)		! reldir26(4,1)
endif
if (g(3) < 0) then		! -z
	w(13) = -g(3)		! reldir26(5,1)
else					! +z
	w(15) = g(3)		! reldir26(6,1)
endif
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
logical function canstick(stick1,stick2)
real :: stick1, stick2
real :: R
integer :: kpar = 0

R = par_uni(kpar)
if (R < stick1*stick2) then
	canstick = .true.
else
	canstick = .false.
endif
end function

!---------------------------------------------------------------------
! mono(icell1) is sticky, may stick to an adjacent sticky cell to
! either join a clump or initiate a new clump.
!---------------------------------------------------------------------
subroutine sticker(icell1)
integer :: icell1
type(monocyte_type), pointer :: cell1, cell2
!type(clump_type), pointer :: clump1, clump2
integer :: site1(3), site2(3), i, icell2, iclump1, iclump2
!logical :: isclump1, isclump2

cell1 => mono(icell1)
site1 = cell1%site
if (site1(2) == NBY+1) then
	write(*,*) 'sticker (a): ',cell1%ID,site1
	stop
endif
iclump1 = cell1%iclump
do i = 1,27
	if (i == 14) cycle
	site2 = site1 + jumpvec(:,i)
	if (site2(1) < 1 .or. site2(1) > NX) cycle
!	if (site2(2) < NBY+1 .or. site2(2) > NY) cycle
	if (site2(2) <= NBY+1 .or. site2(2) > NY) cycle
	if (site2(3) < 1 .or. site2(3) > NZ) cycle
	icell2 = occupancy(site2(1),site2(2),site2(3))%indx
	if (icell2 > 0) then
		cell2 => mono(icell2)
		if (cell2%stickiness > 0 .and. cell2%status < FUSED) then
			iclump2 = cell2%iclump
			if ((iclump1 > 0) .and. (iclump1 == iclump2)) cycle		! already stuck together
			if (canstick(cell1%stickiness,cell2%stickiness)) then
				call stick(icell1,icell2)
				iclump1 = cell1%iclump
			endif
		endif
	endif
enddo
if (cell1%site(2) == NBY+1) then
	write(*,*) 'sticker (b): ',cell1%ID,site1
	stop
endif
end subroutine

!------------------------------------------------------------------------------------------------
! A monocyte integrates RANK receptor signal (which also decays) and when it reaches a threshold level
! the cell is activated and becomes sticky.  
! Two threshold approach: 
! Cell integrates RANK receptor signal (S).
! When S > ST1 the cell starts to respond to the attracting chemotactic signal from OBs.
! As S increases the response to chemotactic signal increases, and the S1P1 level declines.
! When S > ST2 the cell becomes sticky.
! When two sticky cells meet they stick (with a probability determined by their stickiness levels).
! While stuck together they lose random motility.
! The clump of stuck cells grows, and when the number in the clump exceeds another threshold 
! the cells fuse.
! The fused cells move towards the bone under the influence of the attracting signal.  
! Additional sticky monocytes coming into contact with the ball of fused cells will also fuse with it.
! On reaching the bone surface, the fused cells enter the final stage of differentiation into a mature osteoclast.
!------------------------------------------------------------------------------------------------
subroutine stick(icell1,icell2)
integer :: icell1, icell2
type(monocyte_type), pointer :: cell1, cell2
type(clump_type), pointer :: pclump
real :: tnow

cell1 => mono(icell1)
cell2 => mono(icell2)
tnow = istep*DELTA_T
if (cell1%iclump > 0) then
	if (cell2%iclump > 0) then
		! both cell1 and cell2!
!		write(logmsg,*) 'stick: did this really happen? ',icell1,icell2,cell1%iclump,cell2%iclump
!		call logger(logmsg)
		if (clump(cell1%iclump)%ncells + clump(cell2%iclump)%ncells <= MAX_CLUMP_CELLS) then
			call join_clump(icell1,icell2)
		endif
	else
		! just cell1
		pclump => clump(cell1%iclump)
		if (pclump%ncells == MAX_CLUMP_CELLS) then
!			write(logmsg,*) 'Too many cells in clump (1)'
!			call logger(logmsg)
			return
!			stop
		endif
		cell2%status = CLUMPED
		cell2%iclump = cell1%iclump
		pclump%ncells = pclump%ncells + 1
		pclump%list(pclump%ncells) = icell2
		call centre_of_mass(pclump)
		if (pclump%status >= FUSING) then
			cell2%status = pclump%status
		endif
!		write(logmsg,*) 'Cell added to clump: status: ',cell1%iclump,pclump%status,pclump%ncells
!		call logger(logmsg)
	endif
elseif (cell2%iclump > 0) then
		! just cell2
	pclump => clump(cell2%iclump)
	if (pclump%ncells == MAX_CLUMP_CELLS) then
!		write(logmsg,*) 'Too many cells in clump (2)'
!		call logger(logmsg)
		return
!		stop
	endif
	cell1%status = CLUMPED
	cell1%iclump = cell2%iclump
	pclump%ncells = pclump%ncells + 1
	pclump%list(pclump%ncells) = icell1
	call centre_of_mass(pclump)
	if (pclump%status >= FUSING) then
		cell1%status = pclump%status
	endif
!	write(logmsg,*) 'Cell added to clump: status: ',cell2%iclump,pclump%status,pclump%ncells
!	call logger(logmsg)
else
		! neither cell is in a clump
	nclump = nclump + 1
	if (nclump > MAX_NCLUMP) then
		write(logmsg,*) 'Too many clumps'
		call logger(logmsg)
		stop
	endif
	pclump => clump(nclump)
	pclump%ID = nclump
	cell1%status = CLUMPED
	cell2%status = CLUMPED
	cell1%iclump = nclump
	cell2%iclump = nclump
	pclump%ncells = 2
	pclump%list(1) = icell1
	pclump%list(2) = icell2
	call centre_of_mass(pclump)
	pclump%starttime = tnow
	pclump%status = ALIVE
!	write(logmsg,*) 'New clump'
!	call logger(logmsg)
endif
end subroutine

!---------------------------------------------------------------------
! Cells icell1 and icell2 are both in separate clumps, which must be 
! joined.  One clump is extended, the other has status -> -1
!---------------------------------------------------------------------
subroutine join_clump(icell1,icell2)
integer :: icell1, icell2
integer :: i, j, kcell
type(clump_type), pointer :: pclump1, pclump2

pclump1 => clump(mono(icell1)%iclump)
pclump2 => clump(mono(icell2)%iclump)
write(logmsg,'(a,6i6)')'Joining two clumps: ',icell1,icell2,mono(icell1)%iclump,mono(icell2)%iclump,pclump1%ncells,pclump2%ncells
call logger(logmsg)
do i = 1,pclump2%ncells
	kcell = pclump2%list(i)
	! First check that this isn't in pclump1 - that would be an error
	do j = 1,pclump1%ncells
		if (pclump1%list(j) == kcell) then
			write(logmsg,*) 'Error: joinclumps: cell in both clumps: ',kcell
			call logger(logmsg)
			stop
		endif
	enddo
	! Now add the cell to pclump1
	pclump1%ncells = pclump1%ncells + 1
	pclump1%list(pclump1%ncells) = kcell
	mono(kcell)%iclump = mono(icell1)%iclump
enddo
call centre_of_mass(pclump1)
pclump2%status = -1
end subroutine

!---------------------------------------------------------------------
! The cells in a clump are encouraged to aggregate more closely together.
!---------------------------------------------------------------------
subroutine consolidate_clump(pclump)
type(clump_type), pointer :: pclump
real :: cm(3), dist2(MAX_CLUMP_CELLS), r(3), d2, d2min
integer :: i, k, kmin, kcell, site1(3), site2(3)
type(monocyte_type), pointer :: pmono

!cm = 0
!do i = 1,pclump%ncells
!	cm = cm + mono(pclump%list(i))%site
!enddo
!cm = cm/pclump%ncells
cm = pclump%cm
do i = 1,pclump%ncells
	r = cm - mono(pclump%list(i))%site
	dist2(i) = dot_product(r,r)
enddo
do i = 1,pclump%ncells
	kcell = pclump%list(i)
	site1 = mono(kcell)%site
	d2min = dist2(i)
	kmin = 0
	do k = 1,27
		if (k == 14) cycle
		site2 = site1 + jumpvec(:,k)
		if (site2(1) < 1 .or. site2(1) > NX) cycle
		if (site2(2) <= NBY+1 .or. site2(2) > NY) cycle
		if (site2(3) < 1 .or. site2(3) > NZ) cycle
		if (occupancy(site2(1),site2(2),site2(3))%indx == 0) then	! this site is available
			r = cm - site2
			d2 = dot_product(r,r)
			if (d2 < d2min) then
				d2min = d2
				kmin = k
			endif
		endif
	enddo
	if (kmin /= 0) then		! move the cell to this site
		site2 = site1 + jumpvec(:,kmin)
!		if (pclump%list(i) == 138) then
!			write(*,'(a,7i5)') 'consolidate mono(138): ',i,site1,site2
!		endif
		mono(kcell)%site = site2
		occupancy(site1(1),site1(2),site1(3))%indx = 0
		occupancy(site2(1),site2(2),site2(3))%indx = kcell
	endif
	if (mono(kcell)%site(2) == NBY+1) then
		write(*,*) 'consolidate_clump: ',mono(kcell)%ID,mono(kcell)%site
		stop
	endif
enddo
call centre_of_mass(pclump)
			
end subroutine

!---------------------------------------------------------------------
! Move a clump (ncells = n0) if it is too close to another (big) clump.
! Move away from another clump (ncells = n) if:
!	n0 + n > MAX_CLUMP_CELLS
!---------------------------------------------------------------------
subroutine separate_clump(pclump0)
type(clump_type), pointer :: pclump0
integer :: iclump0, iclump, i, imax, jump(3), site2(3), site0(3), indx, iclast
type(clump_type), pointer :: pclump
real :: r(3), d, rsum(3), proj, pmax
logical :: ok

iclump0 = pclump0%ID
rsum = 0
do iclump = 1,nclump
	if (iclump == iclump0) cycle
	pclump => clump(iclump)
	if (pclump%status == DEAD) cycle
	if (pclump0%ncells + pclump%ncells <= MAX_CLUMP_CELLS) cycle	! the clumps can join
	r = pclump0%cm - pclump%cm
	r(2) = 0
	d = sqrt(dot_product(r,r))
	if (d > CLUMP_SEPARATION) cycle
!	write(logmsg,'(a,2i4,4f6.1)') 'near clump: ',iclump0,iclump,r,d
!	call logger(logmsg)
	rsum = rsum + r
enddo
! Turn off movement away from clasts
!do iclast = 1,nclast
!	if (clast(iclast)%status == DEAD) cycle
!	r = pclump0%cm - clast(iclast)%cm
!	r(2) = 0
!	d = sqrt(dot_product(r,r))
!	if (d > CLUMP_SEPARATION) cycle
!	rsum = rsum + r
!enddo
	
if (rsum(1) == 0 .and. rsum(2) == 0 .and. rsum(3) == 0) return
! Need to try to move pclump0.  Get unit vector in the direction to move
r = rsum/sqrt(dot_product(rsum,rsum))
! Choose jump direction (parallel to bone surface, i.e. dy = 0) closest to r
pmax = 0
imax = 0
do i = 1,27
	if (i == 14) cycle
	if (unitjump(2,i) /= 0) cycle	! consider only jumps with dy = 0
	proj = dot_product(r,unitjump(:,i))
	if (proj > pmax) then
		pmax = proj
		imax = i
	endif
enddo
if (imax == 0) return
jump = jumpvec(:,imax)
! Now see if it will be possible to move all cells in the clump by this jump
do i = 1,pclump0%ncells
	site2 = mono(pclump0%list(i))%site + jump
	ok = .false.
	if (site2(1) < 1 .or. site2(1) > NX) exit
	if (site2(2) < NBY+1 .or. site2(2) > NY) exit
	if (site2(3) < 1 .or. site2(3) > NZ) exit
	if (occupancy(site2(1),site2(2),site2(3))%region /= MARROW) exit
	ok = .true.
	indx = occupancy(site2(1),site2(2),site2(3))%indx
	if (indx == 0) cycle	! this site is available
	if (indx > 0) then
		if (mono(indx)%iclump == iclump0) then	! cell in this clump
			cycle
		else									! another cell in the way
			ok = .false.
			exit
		endif
	endif
enddo
if (.not.ok) return
! We can move all the cells by jump
do i = 1,pclump0%ncells
	site0 = mono(pclump0%list(i))%site
	occupancy(site0(1),site0(2),site0(3))%indx = 0
	mono(pclump0%list(i))%site = site0 + jump
!	if (pclump0%list(i) == 138) then
!		write(*,'(a,7i5)') 'separate mono(138): ',i,site0,site0 + jump
!	endif
enddo
do i = 1,pclump0%ncells
	site0 = mono(pclump0%list(i))%site
	occupancy(site0(1),site0(2),site0(3))%indx = pclump0%list(i)
enddo
pclump0%cm = pclump0%cm + jump
!write(logmsg,'(a,i3,3f5.1,3i3)') 'Moved clump: ',iclump0,r,jump
!call logger(logmsg)
end subroutine

!---------------------------------------------------------------------
! A clump moves towards the bone surface until it is just touching
! (least value of a mono(:)%site(2) = NBY+1), then it transforms into
! an osteoclast
!---------------------------------------------------------------------
subroutine lower_clump(pclump,on_surface)
type(clump_type), pointer :: pclump
logical :: on_surface
real :: R
integer :: i, j, site(3), site2(3), y, ylo, yhi, nhit,hitcell(MAX_CLUMP_CELLS)
integer :: kpar = 0

on_surface = .false.
R = par_uni(kpar)
if (R > CLUMP_FALL_PROB) return
ylo = NY
yhi = 0
nhit = 0
do i = 1,pclump%ncells
	j = pclump%list(i)
	y = mono(j)%site(2)
	if (y == NBY+1) then
		nhit = nhit+1
		hitcell(nhit) = j
	endif
	ylo = min(y,ylo)
	yhi = max(y,yhi)
enddo
if (nhit > 0) then			! there are cells sitting on the surface
	if (yhi > ylo+2) then	! move contacting cells up first before moving clump down
		do i = 1,nhit
			j = hitcell(i)
			site = mono(j)%site
			occupancy(site(1),site(2),site(3))%indx = 0		! clear the old site
			! Need to find a suitable site to move the cell to
			call get_vacant_site(pclump,site(2)+1,site2)
			mono(j)%site = site2
			occupancy(site2(1),site2(2),site2(3))%indx = j
!			write(logmsg,*) 'moved up: ',j,site,' -> ',site2
!			call logger(logmsg)
		enddo
	else					! ready to change to complete transformation into an osteoclast
!		write(logmsg,*) 'reshape_clump: ylo,yhi: ',ylo,yhi
!		call logger(logmsg)
		call reshape_clump(pclump)	! ensure that cells are well placed to make the osteoclast footprint
		on_surface = .true.
		return
	endif
endif
! Now move the clump down
do i = 1,pclump%ncells
	j = pclump%list(i)
	site = mono(j)%site
	occupancy(site(1),site(2),site(3))%indx = 0
	mono(j)%site(2) = site(2) - 1
enddo
do i = 1,pclump%ncells
	j = pclump%list(i)
	site = mono(j)%site
	occupancy(site(1),site(2),site(3))%indx = j
enddo
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine break_clumps
type(clump_type), pointer :: pclump
integer :: iclump, i, kcell

do iclump = 1,nclump
	pclump => clump(iclump)
	if (pclump%status == DEAD) cycle
	do i = 1,pclump%ncells
		kcell = pclump%list(i)
		mono(kcell)%status = DEAD
		mono(kcell)%iclump = 0
	enddo
	pclump%status = DEAD
enddo
nclump = 0
end subroutine

!---------------------------------------------------------------------
! Need a vacant site in the plane at y near the clump
!---------------------------------------------------------------------
subroutine get_vacant_site(pclump,y,site)
type(clump_type), pointer :: pclump
integer :: y, site(3)
integer :: x0, z0, dx, dz, x, z, d2, d2min

x0 = pclump%cm(1)
z0 = pclump%cm(3)
d2min = 99999
do dx = -4,4
	do dz = -4,4
		x = x0 + dx
		z = z0 + dz
		if (x < 1 .or. x > NX) cycle
		if (z < 1 .or. z > NZ) cycle
		if (occupancy(x,y,z)%indx /= 0) cycle
		d2 = dx*dx + dz*dz
		if (d2 < d2min) then
			d2min = d2
			site = (/ x,y,z /)
		endif
	enddo
enddo
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine reshape_clump(pclump)
type(clump_type), pointer :: pclump
integer :: n(2), nt(2), i, j, k, site(3), x0, y0, z0, x, y, z, layer, rng1, rng2

x0 = pclump%cm(1)
y0 = NBY+1
z0 = pclump%cm(3)
n(1) = pclump%ncells/2 + 5
n(2) = pclump%ncells - n(1)
nt(1) = n(1)
nt(2) = n(1) + n(2)
!write(*,*) 'reshape_clump: ',pclump%ncells,n
do i = 1,pclump%ncells
	j = pclump%list(i)
	site = mono(j)%site
	occupancy(site(1),site(2),site(3))%indx = 0
enddo
do layer = 1,2
	y = y0 + layer - 1
	if (n(layer) < 25) then
		rng1 = 1
		rng2 = 2
	else
		rng1 = 2
		rng2 = 3
	endif
	if (layer == 1) then
		i = 0
	else
		i = n(1)
	endif
!	write(*,*) 'layer: ',layer,rng1,rng2
	! Place the first block of 9 or 25
	do x = x0-rng1,x0+rng1
		do z = z0-rng1,z0+rng1
			if (x < 1 .or. x > NX) cycle
			if (z < 1 .or. z > NZ) cycle
			if (occupancy(x,y,z)%indx /= 0) cycle
			site = (/ x,y,z /)
			i = i+1
			j = pclump%list(i)
!			write(*,*) i,j,pclump%ncells
			mono(j)%site = site
			occupancy(site(1),site(2),site(3))%indx = j
			if (i == nt(layer)) exit
		enddo
		if (i == nt(layer)) exit
	enddo
	if (i == nt(layer)) cycle
	! Place the rest of this layer
	do x = x0-rng2,x0+rng2
		do z = z0-rng2,z0+rng2
			if (abs(x-x0) < rng2 .and. abs(z-z0) < rng2) cycle
			if (x < 1 .or. x > NX) cycle
			if (z < 1 .or. z > NZ) cycle
			if (occupancy(x,y,z)%indx /= 0) cycle
			site = (/ x,y,z /)
			i = i+1
			j = pclump%list(i)
!			write(*,*) i,j,pclump%ncells,nt(layer)
			mono(j)%site = site
			occupancy(site(1),site(2),site(3))%indx = j
			if (i == nt(layer)) exit
		enddo
		if (i == nt(layer)) exit
	enddo
enddo
!do i = 1,pclump%ncells
!	j = pclump%list(i)
!	write(*,*) i,j,mono(j)%site
!enddo
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine centre_of_mass(pclump)
type(clump_type), pointer :: pclump
real :: cm(3)
integer :: i

cm = 0
do i = 1,pclump%ncells
	cm = cm + mono(pclump%list(i))%site
enddo
pclump%cm = cm/pclump%ncells
end subroutine

!---------------------------------------------------------------------
! For a candidate destination site, returns .true. if the site is free,
! else returns .false.  Wrapping at the boundaries is handled.
! A site may be unavailable because region /= MARROW, 
! or because kcell /= 0
!---------------------------------------------------------------------
logical function free_site(site,region,kcell)
integer :: site(3), region, kcell

if (site(1) < 1) then
	free_site = wrapx1(site,region,kcell)
	return
endif
if (site(1) > NX) then
	free_site = wrapx2(site,region,kcell)
	return
endif
if (site(2) > NY) then
	free_site = wrapy2(site,region,kcell)
	return
endif
if (site(3) < 1) then
	free_site = wrapz1(site,region,kcell)
	return
endif
if (site(3) > NZ) then
	free_site = wrapz2(site,region,kcell)
	return
endif
!if (site(2) <= NBY+1) then	!???????????  what's wrong with this?
!	kcell = 0
!	region = LAYER
!	free_site = .false.
!	return
!endif
region = occupancy(site(1),site(2),site(3))%region
if (region /= MARROW) then
	kcell = 0
	free_site = .false.
!	write(*,*) 'bad region'
	return
endif
kcell = occupancy(site(1),site(2),site(3))%indx
if (kcell == 0) then
	free_site = .true.
else
	free_site = .false.
!	write(*,*) 'occupied'
endif
end function

!---------------------------------------------------------------------
! site(1) = 0, and the location is wrapped to site(1) = NX
! The wrapped site may be an inaccessible region (non-MARROW), 
! or it may be occupied by a cell.
! Returns .true. for a vacant site, .false. otherwise
! site is replaced by the wrapped location, and the region and 
! cell index (kcell) are also returned as arguments.
!---------------------------------------------------------------------
logical function wrapx1(site,region,kcell)
integer :: site(3), region, kcell

site(1) = NX
site(2) = min(site(2),NY)
site(3) = min(site(3),NZ)
site(3) = max(site(3),1)
region = occupancy(site(1),site(2),site(3))%region
if (region /= MARROW) then
	kcell = 0
	wrapx1 = .false.
	return
endif
kcell = occupancy(site(1),site(2),site(3))%indx
if (kcell == 0) then
	wrapx1 = .true.
else
	wrapx1 = .false.
endif
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
logical function wrapx2(site,region,kcell)
integer :: site(3), region, kcell

site(1) = 1
site(2) = min(site(2),NY)
site(3) = min(site(3),NZ)
site(3) = max(site(3),1)
region = occupancy(site(1),site(2),site(3))%region
if (region /= MARROW) then
	kcell = 0
	wrapx2 = .false.
	return
endif
kcell = occupancy(site(1),site(2),site(3))%indx
if (kcell == 0) then
	wrapx2 = .true.
else
	wrapx2 = .false.
endif
end function

!---------------------------------------------------------------------
! This case is different because instead of wrapping to the opposite 
! side of the cube, the cell enters the top face (y = NY) at another
! location.  The selection is a bit arbitrary.  We choose to map 
! half of the surface, divided at x = NX/2, into the other half.
! z is preserved.
!---------------------------------------------------------------------
logical function wrapy2(site,region,kcell)
integer :: site(3), region, kcell

site(2) = NY
site(1) = min(site(1),NX)
site(1) = max(site(1),1)
site(3) = min(site(3),NZ)
site(3) = max(site(3),1)
if (site(1) > NX/2) then
	site(1) = site(1) - NX/2
else
	site(1) = site(1) + NX/2
endif	
region = occupancy(site(1),site(2),site(3))%region
if (region /= MARROW) then
	kcell = 0
	wrapy2 = .false.
	return
endif
kcell = occupancy(site(1),site(2),site(3))%indx
if (kcell == 0) then
	wrapy2 = .true.
else
	wrapy2 = .false.
endif
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
logical function wrapz1(site,region,kcell)
integer :: site(3), region, kcell

site(3) = NZ
site(2) = min(site(2),NY)
site(1) = min(site(1),NX)
site(1) = max(site(1),1)
region = occupancy(site(1),site(2),site(3))%region
if (region /= MARROW) then
	kcell = 0
	wrapz1 = .false.
	return
endif
kcell = occupancy(site(1),site(2),site(3))%indx
if (kcell == 0) then
	wrapz1 = .true.
else
	wrapz1 = .false.
endif
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
logical function wrapz2(site,region,kcell)
integer :: site(3), region, kcell

site(3) = 1
site(2) = min(site(2),NY)
site(1) = min(site(1),NX)
site(1) = max(site(1),1)
region = occupancy(site(1),site(2),site(3))%region
if (region /= MARROW) then
	kcell = 0
	wrapz2 = .false.
	return
endif
kcell = occupancy(site(1),site(2),site(3))%indx
if (kcell == 0) then
	wrapz2 = .true.
else
	wrapz2 = .false.
endif
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
logical function crossToBlood(kcell,site)
real :: R, prob
integer :: kcell, site(3), kpar=0
type(monocyte_type), pointer :: cell
real :: tnow

cell => mono(kcell)
crossToBlood = .false.
if (cell%S1P1 > S1P1_THRESHOLD) then
	prob = CROSS_PROB*(cell%S1P1 - S1P1_THRESHOLD)/(1 - S1P1_THRESHOLD)
	R = par_uni(kpar)
	if (R < prob) then
		tnow = istep*DELTA_T
		crossToBlood = .true.
		cell%status = CROSSING
		cell%exittime = tnow + CROSSING_TIME
!		write(logmsg,*) 'S1P1, prob, R: ',cell%S1P1,prob,R,cell%exittime
!		call logger(logmsg)
	endif
endif
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine chemo(S,site,p,jmpdir,pstar)
integer :: site(3), jmpdir(MAXRELDIR)
real :: S, p(0:MAXRELDIR), pstar(0:MAXRELDIR)
integer :: u(3), i
real :: g(3), gdotu, psum

g = S1P_grad(:,site(1),site(2),site(3))
do i = 1,MAXRELDIR
	if (p(i) == 0) cycle
	u = jumpvec(:,jmpdir(i))
	gdotu = g(1)*u(1) + g(2)*u(2) + g(3)*u(3)
	pstar(i) = p(i) + Kchemo*S*gdotu
	pstar(i) = max(pstar(i),0.0)
	pstar(i) = min(pstar(i),1.0)
	psum = psum + pstar(i)
enddo
if (psum > 1) then
	pstar(1:MAXRELDIR) = pstar(1:MAXRELDIR)/psum
	pstar(0) = 0
else
	pstar(0) = 1 - psum
endif
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
real function clast_size(count)
integer :: count
clast_size = CLAST_RADIUS_FACTOR*sqrt(real(count))
end function

!---------------------------------------------------------------------
! Choose one of the 6 axis directions to allocate to the direction of
! previous step given by jump, which takes values from 0 - 27.  
! This could be already on an axis, or could be a D2 or D3 diagonal, 
! in which case the choice is random.
!---------------------------------------------------------------------
integer function fix_lastdir(jump,kpar)
integer :: jump, kpar
integer :: k,nax
!                            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
integer :: naxes(0:27) = (/  0, 3, 2, 3, 2, 1, 2, 3, 2, 3, 2, 1, 2, 1, 0, 1, 2, 1, 2, 3, 2, 3, 2, 1, 2, 3, 2, 3 /)
integer :: axes(3,27) = reshape( (/ &
!   1      2      3      4      5      6      7      8      9     10     11     12     13     14
1,3,5, 1,3,0, 1,3,6, 1,5,0, 1,0,0, 1,6,0, 1,4,5, 1,4,0, 1,4,6, 3,5,0, 3,0,0, 3,6,0, 5,0,0, 0,0,0, &
!  15     16     17     18     19     20     21     22     23     24     25     26     27 
6,0,0, 4,5,0, 4,0,0, 4,6,0, 2,3,5, 2,3,0, 2,3,6, 2,5,0, 2,0,0, 2,6,0, 2,4,5, 2,4,0, 2,4,6 /), (/3,27/) )

if (diagonal_jumps) then
	nax = naxes(jump)
	if (nax == 0) then
	    write(*,*) 'Should not get here: fix_lastdir: nax=0: jump: ',jump
	    stop
		fix_lastdir = random_int(1,6,kpar)
	elseif (nax == 1) then
		fix_lastdir = axes(1,jump)
	else
		k = random_int(1,nax,kpar)
		fix_lastdir = axes(k,jump)
	endif
else
	write(*,*) 'fix_lastdir: Not MOORE model: ',jump
	stop
endif
end function

!---------------------------------------------------------------------
! Need to precompute array reldir(:,:)
! FOR NRELDIR = 6 (diagonal_jumps = .false.)
! The directions 1 - 6 are defined according to the neumann array,
! i.e. 1 is -x, 2 is +x, 3 is -y, 4 is +y, 5 is -z, 6 is +z
! The relative directions are numbered as follows:
! irel = 1 gives the same direction as lastdir
! irel = 6 gives the opposite direction to lastdir
! irel = 2,3,4,5 cover the remaining directions - order not important
! at the moment since all directions normal to lastdir are equally likely
! FOR NRELDIR = 17 (diagonal_jumps = .true.)
! The reldir(:,:) array uses the same set of 6 previous jump directions,
! restricted to the axes.  Diagonal previous jumps are accommodated by
! making a random selection of an axis direction from either 2 or 3
! possibilities (depending on whether it was a D2 or D3 jump).
! The possible jumps, now including diagonal moves, are stored in jumpvec(:)
! The jumpvec(:) entries are ordered with z varying fastest, then y, then x,
! each ranging -1, 0, +1.
! For dirprob(:) calculation:
! There are 3 groups of jumps, with sets of probs
! Group 1: Q(1) = P(1), D1 jump in the same as last direction
! Group 2: Q(2) = P(2)+P(3)+P(4)+P(5) (D2 jumps) + P(6)+P(7)+P(8)+P(9) (D3 jumps)
! Group 3: Q(3) = P(10)+P(11)+P(12)+P(13) (D1 jumps) + P(14)+P(15)+P(16)+P(17) (D2 jumps)
! i.e. Q(1) = prob of no direction change, Q(2) = prob of direction change < 90 deg,
! Q(3) = prob of direction change = 90 deg.
! Setting D1 jumps (von Neumann) to 1, the D2 jumps have length L2 = sqrt(2)
! and the D3 jumps have length L3 = sqrt(3)
! Therefore jump probs must be scaled appropriately to give symmetry within a set
! In other words:
!		P(i)/P(j) = L2/L3 for j=6,7,8,9 j=2,3,4,5
! and	P(i)/P(j) = 1/L2  for i=14,15,16,17 j=10,11,12,13
! Then choosing P(2) to represent D2 jump prob in Group 2, and P(10) to
! represent D1 jump prob in Group 3:
!		Q(2) = 4(1 + L2/L3).P(2)
!		Q(3) = 4(1 + 1/L2).P(10)
! It is always the case that Q(1) + Q(2) + Q(3) = BETA
!
! When RHO = 1, Q(1) = BETA = p1, Q(2) + Q(3) = 0
!
! When RHO = 0, P(10) = P(1), and P(2) = P(1)/L2
! therefore:	Q(3) = 4(1 + 1/L2).P(1)
! and			Q(2) = 4(1 + L2/L3).P(1)/L2
! giving:		P(1).[1 + 4(1 + L2/L3)/L2 + 4(1 + 1/L2)] = BETA
!				P(1) = BETA/[1 + 4(1 + L2/L3)/L2 + 4(1 + 1/L2)] = p0
! and:			Q(3)/Q(2) = (1 + L2)/(1 + L2/L3) = q320
!
! Now for 0 < RHO < 1
! first set Q(1) = P(1) to vary linearly with RHO from 
! the value at 0 (p0) to the value at 1 (p1)
!	P(1) = p0 + RHO*(p1-p0)
! Now we know that Q(2) + Q(3) = BETA - Q(1)
! and the third parameter ALPHA is used to determine how Q(3)/Q(2)
! varies with RHO, by making it change linearly from q23 at RHO=0
! to ALPHA*q23 at RHO = 1
! Now we have all the info to solve for Q(1), Q(2), Q(3) and all P(:)
!
! Revised Moore model.
! Now allow reverse directions, but disallow D3 diagonal jumps, giving
! 18 possible directions.
! In the preferred direction the surface is a prolate spheroid, with
! a = 1, b = 1-RHO.
! In the reverse direction the surface is either a sphere (a = b) or
! an oblate spheroid (a = b^2).
!---------------------------------------------------------------------
subroutine make_reldir
integer :: lastdir,irel,k,ix,iy,iz,i,site(3)
integer :: reldir18(6,18),reldir26(6,26)
real :: Q(3),L2,L3,p0,p1,q320,q321,q32,qsum,a,b,c,d(MAXRELDIR),R(MAXRELDIR),E, theta, p(8), psum, u2, u(3)

if (MODEL == NEUMANN_MODEL) then
	diagonal_jumps = .false.
	njumpdirs = 6
	nreldir = 6
	reldir = 0
	do lastdir = 1,6
		reldir(lastdir,1) = lastdir
		if (mod(lastdir,2) == 0) then
			reldir(lastdir,6) = lastdir - 1
		else
			reldir(lastdir,6) = lastdir + 1
		endif
		irel = 2
		do k = 1,6
			if (k /= reldir(lastdir,1) .and. k /= reldir(lastdir,6)) then
				reldir(lastdir,irel) = k
				irel = irel + 1
			endif
		enddo		
	enddo
!   	write(*,*) 'reldir'
!    do lastdir = 1,6
!	    write(*,'(6i6)') (reldir(lastdir,k),k=1,6)
!    enddo
	jumpvec(:,1:6) = neumann(:,1:6)

else
	if (MODEL == MOORE18_MODEL) then
    	nreldir = 18
	elseif (MODEL == MOORE26_MODEL) then
    	nreldir = 26
    endif
    njumpdirs = 27
	diagonal_jumps = .true.
	k = 0
	do ix = -1,1
		do iy = -1,1
			do iz = -1,1
				k = k+1
				jumpvec(:,k) = (/ix,iy,iz/)
				u = jumpvec(:,k)
				u2 = dot_product(u,u)
				if (u2 > 0) then
					unitjump(:,k) = u/sqrt(u2)
				else
					unitjump(:,k) = 0
				endif
			enddo
		enddo
	enddo

! Data for revised Moore18 model.  D3 jumps are excluded, but reverse directions are allowed.
	reldir18(1,:) = (/  5,  2, 4, 6, 8, 11,13,15,17, 10,12,16,18, 22,24,26,20, 23 /)	! -x
	reldir18(2,:) = (/ 23, 20,22,24,26, 11,13,15,17, 10,12,16,18,  2, 4, 6, 8,  5 /)	! +x
	reldir18(3,:) = (/ 11,  2,10,12,20,  5,13,15,23,  4, 6,22,24,  8,16,18,26, 17 /)	! -y
	reldir18(4,:) = (/ 17,  8,16,18,26,  5,13,15,23,  4, 6,22,24,  2,10,12,20, 11 /)	! +y
	reldir18(5,:) = (/ 13,  4,10,16,22,  5,11,17,23,  2, 8,20,26,  6,12,18,24, 15 /)	! -z
	reldir18(6,:) = (/ 15,  6,12,18,24,  5,11,17,23,  2, 8,20,26,  4,10,16,22, 13 /)	! +z

! Data for revised Moore26 model.  D3 jumps are included, but reverse directions are allowed.
	reldir26(1,:) = (/  5,  2, 4, 6, 8,  1, 3, 7, 9, 11,13,15,17, 10,12,16,18, 19,21,27,25, 20,22,24,26, 23 /)	! -x
	reldir26(2,:) = (/ 23, 20,22,24,26, 19,21,27,25, 11,13,15,17, 10,12,16,18,  1, 3, 7, 9,  2, 4, 6, 8,  5 /)	! +x
	reldir26(3,:) = (/ 11,  2,10,12,20,  1, 3,19,21,  5,13,15,23,  4, 6,22,24,  7, 9,25,27,  8,16,18,26, 17 /)	! -y
	reldir26(4,:) = (/ 17,  8,16,18,26,  7, 9,25,27,  5,13,15,23,  4, 6,22,24,  1, 3,19,21,  2,10,12,20, 11 /)	! +y
	reldir26(5,:) = (/ 13,  4,10,16,22,  1, 7,19,25,  5,11,17,23,  2, 8,20,26,  3, 9,21,27,  6,12,18,24, 15 /)	! -z
	reldir26(6,:) = (/ 15,  6,12,18,24,  3, 9,21,27,  5,11,17,23,  2, 8,20,26,  1, 7,19,25,  4,10,16,22, 13 /)	! +z

    if (MODEL == MOORE18_MODEL) then
        reldir(1:6,1:18) = reldir18
    elseif (MODEL == MOORE26_MODEL) then
        reldir = reldir26
    endif

endif
call compute_dirprobs
!write(*,*)  'dirprob: '
!write(*,'(7f6.3)') dirprob(0:nreldir)
qsum = 0
do i = 0,nreldir
    qsum = qsum + dirprob(i)
enddo
call make_vn_adjacent
!write(*,*) 'sum: ',qsum
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine make_vn_adjacent
integer :: dir, v(3)

do dir = 1,njumpdirs
	v = jumpvec(:,dir)
	if ((v(1) == 0 .and. v(2) == 0 .and. v(3) /= 0) .or. &
	    (v(1) == 0 .and. v(2) /= 0 .and. v(3) == 0) .or. &
	    (v(1) /= 0 .and. v(2) == 0 .and. v(3) == 0)) then
	   vn_adjacent(dir) = .true.
	else
		vn_adjacent(dir) = .false.
	endif
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine compute_dirprobs
integer :: i
real :: L2,L3,b,c,e,d(MAXRELDIR),R(MAXRELDIR)

if (MODEL == NEUMANN_MODEL) then
!	diagonal_jumps = .false.
	nreldir = 6
	! New Neumann model (with reverses)
	b = 1-RHO
	b = b*b
	e = BETA/(1 + 4*b + b**2)
	dirprob(0) = 1 - BETA	! this is the prob of no jump
	dirprob(1) = e
	dirprob(2:5) = e*b
	dirprob(6) = e*b**2
else
    ! New prolate spheroid approach
!	diagonal_jumps = .true.
    if (MODEL == MOORE18_MODEL) then
	    nreldir = 18
	    L2 = sqrt(2.0)
	    b = 1 - RHO
	    b = b*b		! revised version
	    c = b*b		! oblate spheroid on reverse side
	    d(1) = 1
	    R(1) = 1
	    d(2:5) = L2
	    R(2:5) = L2*b/sqrt(b*b+1)
	    d(6:9) = 1
	    R(6:9) = b
	    d(10:13) = L2
	    R(10:13) = b
	    d(14:17) = L2
	    R(14:17) = L2*b*c/sqrt(b*b+c*c)
	    d(18) = 1
	    R(18) = c
	    e = BETA/(1 + 4*b/sqrt(b*b+1) + 4*b + 4*b/L2 + 4*b*c/sqrt(b*b+c*c) + c)
	    dirprob(0) = 1 - BETA
	    dirprob(1:nreldir) = e*R(1:nreldir)/d(1:nreldir)
    elseif (MODEL == MOORE26_MODEL) then
    	nreldir = 26
	    L2 = sqrt(2.0)
	    L3 = sqrt(3.0)
	    b = 1 - RHO
	    b = b*b		! revised version
	    c = b*b		! oblate spheroid on reverse side
	    d(1) = 1
	    R(1) = 1
	    d(2:5) = L2
	    R(2:5) = L2*b/sqrt(b*b+1)
	    d(6:9) = L3
	    R(6:9) = L2*b/sqrt(b*b+1)
	    d(10:13) = 1
	    R(10:13) = b
	    d(14:17) = L2
	    R(14:17) = b
	    d(18:21) = L3
	    R(18:21) = L2*b*c/sqrt(b*b+c*c)
	    d(22:25) = L2
	    R(22:25) = L2*b*c/sqrt(b*b+c*c)
	    d(26) = 1
	    R(26) = c
	    e = BETA/(1 + 4*(1+L2/L3)*b/sqrt(b*b+1) + 4*b + 4*b/L2 + 4*(1+L2/L3)*b*c/sqrt(b*b+c*c) + c)
	    dirprob(0) = 1 - BETA
	    dirprob(1:nreldir) = e*R(1:nreldir)/d(1:nreldir)
    endif
endif
end subroutine

end module

