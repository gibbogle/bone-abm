module motion
use global
use fields
implicit none

integer :: reldir(6,26)

integer :: MODEL = MOORE26_MODEL
real, parameter :: Kchemo = 1.0

integer :: njumpdirs, nreldir
integer :: jumpvec(3,MAXRELDIR+1)
real :: dirprob(0:MAXRELDIR)


contains

!---------------------------------------------------------------------
! Currently monocytes move around while in the marrow, and may pass
! into the blood, effectively removed from further consideration.
!---------------------------------------------------------------------
subroutine mover
integer :: kcell
logical :: go
integer :: kpar = 0
integer :: kdbug = -1
integer :: site(3)
real :: tnow

tnow = istep*DELTA_T
do kcell = 1,nmono
	if (mono(kcell)%status == MOTILE) then	! interim criterion
		call mono_jumper(kcell,go,kpar)
		if (kcell == kdbug) then
			write(*,'(3i4)') mono(kcell)%site
		endif
	elseif (mono(kcell)%status == CROSSING) then
		if (tnow >= mono(kcell)%exittime) then
			mono(kcell)%status = LEFT
			mono(kcell)%region = BLOOD
			site = mono(kcell)%site
			occupancy(site(1),site(2),site(3))%indx = 0
			mono_cnt = mono_cnt - 1
			nleft = nleft + 1
!			write(*,'(a,2i6,2f6.3,i6)') 'monocyte leaves: ',istep,kcell,cell%S1P1,prob,mono_cnt
		endif
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
integer :: irel,dir1,lastdir1
integer :: savesite2(3,26), jmpdir(26)
real :: psum, p(26), R, f0, f
logical :: free, cross, field

cell => mono(kcell)
site1 = cell%site
go = .false.
field = .false.
if (occupancy(site1(1),site1(2),site1(3))%signal /= 0) then
	field = .true.
	f0 = occupancy(site1(1),site1(2),site1(3))%intensity
endif

!R = par_uni(kpar)
call random_number(R)
if (R <= dirprob(0)) then    ! case of no jump
	return
endif

! Now we must jump (if possible)

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
		p(irel) = dirprob(irel)
		if (field) then
			f = occupancy(site2(1),site2(2),site2(3))%intensity
			if (f >= FTHRESHOLD) then
				p(irel) = max(0.0,f-f0)
			else
				p(irel) = max(0.0,p(irel) + AFACTOR*(f - f0))
			endif
		endif
		jmpdir(irel) = dir1
		psum = psum + p(irel)
		savesite2(:,irel) = site2
	elseif (region == BLOOD) then
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
	return
else
    go = .true.
endif

! Now choose a direction on the basis of these probs p()
!R = par_uni(kpar)
call random_number(R)
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
		write(*,*) 'dir1 = 0: ',kcell,lastdir1,irel,dir1
	endif
endif
if (diagonal_jumps) then
	dir1 = fix_lastdir(dir1,kpar)
elseif (dir1 == 0) then
	dir1 = random_int(1,6,kpar)
endif

cell%site = site2
cell%lastdir = dir1
occupancy(site2(1),site2(2),site2(3))%indx = kcell
occupancy(site1(1),site1(2),site1(3))%indx = 0
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
integer :: kcell, site(3)
type(monocyte_type), pointer :: cell
real :: tnow

cell => mono(kcell)
crossToBlood = .false.
if (cell%S1P1 > S1P1_THRESHOLD) then
	prob = CROSS_PROB*(cell%S1P1 - S1P1_THRESHOLD)/(1 - S1P1_THRESHOLD)
	call random_number(R)
	if (R < prob) then
		tnow = istep*DELTA_T
		crossToBlood = .true.
		cell%status = CROSSING
		cell%exittime = tnow + CROSSING_TIME
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
real :: Q(3),L2,L3,p0,p1,q320,q321,q32,qsum,a,b,c,d(MAXRELDIR),R(MAXRELDIR),E, theta, p(8), psum, dd

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
			enddo
		enddo
	enddo

! Data for revised Moore18 model.  D3 jumps are excluded, but reverse directions are allowed.
	reldir18(1,:) = (/  5,  2, 4, 6, 8, 11,13,15,17, 10,12,16,18, 22,24,26,20, 23 /)	! -x
	reldir18(2,:) = (/ 23, 20,22,24,26, 11,13,15,17, 10,12,16,18,  2, 4, 6, 8,  5 /)	! +x
	reldir18(3,:) = (/ 11,  2,10,12,20,  5,13,15,23,  4, 6,22,24,  8,16,18,26, 17 /)	! -y
	reldir18(4,:) = (/ 17,  8,16,18,26,  5,13,15,23,  4, 6,22,24,  2,10,12,20, 11/)	! +y
	reldir18(5,:) = (/ 13,  4,10,16,22,  5,11,17,23,  2, 8,20,26,  6,12,18,24, 15 /)	! -z
	reldir18(6,:) = (/ 15,  6,12,18,24,  5,11,17,23,  2, 8,20,26,  4,10,16,22, 13 /)	! +z

! Data for revised Moore26 model.  D3 jumps are excluded, but reverse directions are allowed.
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
!write(*,*) 'sum: ',qsum
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

