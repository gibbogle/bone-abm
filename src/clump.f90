! Monocyte clumping procedures
! A clump grows near an OB that is secreting a chemotactic molecule (currently this is
! referred to as RANKL, although possibly a different molecule is involved).
! A monocyte experiencing a RANKL concentration C receives signal from its RANK receptors,
! which it integrates and stores as S.  The rate of signalling is a function of C and S,
! computed by rate_RANKSIGNAL(S,C).  It is the level of S, and the proximity of the monocyte
! to an OB, that determines monocyte differentiation and clumping.
! The stages of this process are:
! MOTILE      --> CHEMOTACTIC  when S > ST1
! CHEMOTACTIC --> STICKY       when S > ST2 and monocyte is near OB (stickiness = S)
! STICKY      --> CLUMPED      when monocyte encounters another sticky mono,
!                              with probability = product of stickiness levels
! Once a clump has formed (at least two monocytes) the monocytes that compose it are
! no longer mobile, although there may be small adjustments to the shape of the clump.
! When the clump has grown big enough the fusing process is initiated, and once it is complete
! the clump moves down onto the bone surface and becomes an osteoclast:
! ALIVE  --> FUSING            when ncells > CLUMP_THRESHOLD
! FUSING --> FUSED             after a time FUSING_TIME
! FUSED  --> DEAD              when the osteoclast is formed, clump is removed
!
! The mechanisms that control the fate of monocytes that fail to join a clump,
! either through insufficient RANK stimulation or because there is no lump to join, 
! and clumps that fail to reach critical size and fuse into osteoclasts,
! both need to be improved.
!---------------------------------------------------------------------

module clumping

use global

implicit none 
save

contains

!---------------------------------------------------------------------
! Not working yet
!---------------------------------------------------------------------
subroutine detacher(pclump)
type(monocyte_type), pointer :: pmono
type(clump_type), pointer :: pclump
integer :: imono, k, site(3)
real(8) :: R
integer :: kpar = 0
real, parameter :: DETACH_PROB = 0.01

R = par_uni(kpar)
if (R > DELTA_T/60) then
	return
endif
imono = pclump%list(pclump%ncells)	! look at last monocyte added to the clump
pmono => mono(imono)
site = pmono%site
R = par_uni(kpar)
if (R < DETACH_PROB) then	! try to release this mono
	do k = 3,5
		if (occupancy(site(1),site(2)+k,site(3))%indx == 0) then	! to this site
			pclump%ncells = pclump%ncells - 1
			if (pclump%ncells == 0) then	! remove the clump
				call RemoveClump(pclump)
			endif
			pmono%site(2) = site(2)+k
			pmono%status = MOTILE
			pmono%iclump = 0
		endif
	enddo
endif
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
logical function canstick(stick1,stick2)
real :: stick1, stick2
real(8) :: R
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
! NOTE:  Need to ensure that an OB is associated with only one clump!!!
! After a clump has transformed into an OC the OB is available to
! start another clump.
!---------------------------------------------------------------------
subroutine sticker(icell1,iblast)
integer :: icell1, iblast
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
				call stick(icell1,icell2,iblast)
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
subroutine stick(icell1,icell2,iblast)
integer :: icell1, icell2, iblast
type(monocyte_type), pointer :: cell1, cell2
type(clump_type), pointer :: pclump
real :: tnow

cell1 => mono(icell1)
cell2 => mono(icell2)
tnow = istep*DELTA_T
if (cell1%iclump > 0) then
	if (cell2%iclump > 0) then
		! both cell1 and cell2 are in a clump - these clumps are too close together!
!		write(logmsg,*) 'stick: did this really happen? ',icell1,icell2,cell1%iclump,cell2%iclump, &
!			clump(cell1%iclump)%iblast,clump(cell2%iclump)%iblast
!		call logger(logmsg)
		return
		if (clump(cell1%iclump)%ncells + clump(cell2%iclump)%ncells <= MAX_CLUMP_CELLS) then
			if (clump(cell1%iclump)%iblast == clump(cell2%iclump)%iblast) then
				call join_clump(icell1,icell2)
			endif
		endif
	else
		! just cell1 is in a clump
		pclump => clump(cell1%iclump)
		if (pclump%ncells == MAX_CLUMP_CELLS) then
			write(logmsg,*) 'Too many cells in clump (1)'
			call logger(logmsg)
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
		! just cell2 is in a clump
	pclump => clump(cell2%iclump)
	if (pclump%ncells == MAX_CLUMP_CELLS) then
		write(logmsg,*) 'Too many cells in clump (2)'
		call logger(logmsg)
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
	if (blast(iblast)%iclump /= 0) then
		! This OB already has a clump attached
		return
	endif
	nclump = nclump + 1
	if (nclump > MAX_NCLUMP) then
		write(logmsg,*) 'Too many clumps'
		call logger(logmsg)
		stop
	endif
	pclump => clump(nclump)
	pclump%ID = nclump
	pclump%iblast = iblast			! create correspondence between clump and OB
	blast(iblast)%iclump = nclump
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
pclump2%ncells = 0
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
	call RemoveClump(pclump)
enddo
nclump = 0
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


!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine RemoveClump(pclump)
type(clump_type), pointer :: pclump
integer :: iblast, k, kcell

! Is this OK?
!do k = 1,pclump%ncells
!	kcell = pclump%list(k)
!	mono(kcell)%status = DEAD
!enddo

pclump%status = DEAD
pclump%ncells = 0
iblast = pclump%iblast
blast(iblast)%iclump = 0
end subroutine

end module