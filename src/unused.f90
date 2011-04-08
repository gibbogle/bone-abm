! UNUSED procedures

!================================================================================================
! OLD CODE
!================================================================================================
!------------------------------------------------------------------------------------------------
! To avoid problems with wrapping, no signal site is allowed to occur near the boundary.
! The bone normal is used to give higher preference to sites near the bone. 
! NOT USED
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
!						occupancy(x,y,z)%signal = isig
!						occupancy(x,y,z)%intensity = 1/(d2 + dot_product(v,vn))
					else
!						occupancy(x,y,z)%signal = 0
!						occupancy(x,y,z)%intensity = 0
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
!				if (p%signal == isig .and. p%intensity >= SIGNAL_THRESHOLD) then
!					if (p%region == MARROW) then
!						nt = nt + 1
!						if (p%indx /= 0) then
!							n = n+1
!						endif
!					endif
!				endif
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
! to make an osteoclast.  
! The list of grid sites below the monocytes is created, and used  to make the list of pits 
! associated with the osteoclast.
! NOT USED
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
!pclast%site = site
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
!			if (p%signal == isig .and. p%intensity >= SIGNAL_THRESHOLD &
!				.and. p%region == MARROW .and. p%indx /= 0) then
			if (p%region == MARROW .and. p%indx /= 0) then
				cnt = cnt+1
!				pclast%mono(cnt) = p%indx
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
!allocate(pclast%pit(npit))
!do i = 1,npit
!	pclast%pit(i)%site = bonesite(:,i)
!enddo
!call pitrates(pclast)
!	clast(nclast)%pit(i)%fraction = 1.0
!	d = sqrt(real((bonesite(1,i)-site(1))**2 + (bonesite(2,i)-site(2))**2 + (bonesite(3,i)-site(3))**2))
!	clast(nclast)%pit(i)%rate = resorptionRate(cnt,d)
!enddo
end subroutine
!================================================================================================

!----------------------------------------------------------------------------------------
! To formulate algebraic system to solve for steady-state solution of diffusion equation
! by finite differences.
! The index of the vector entry corresponding to C at site (x,y,z) is precomputed as
! vindex(x,y,z).  This is a positive integer for (x,y,z) a marrow site.
! Otherwise (a site is in a capillary) it is 0.
! NOT USED
!----------------------------------------------------------------------------------------
subroutine formulate(Kdiffusion,Kdecay)
real :: Kdiffusion, Kdecay
integer, allocatable :: vindex(:,:,:), influx(:,:,:)
real, allocatable :: A(:,:), b(:), V(:)
integer :: x, y, z, k, NV, xn, yn, zn, row, col, col0, offset(3)
logical :: use_nbr(6)
real :: Area, Vol

Area = 1
Vol = 1
allocate(influx(NX,NY,NZ))
allocate(vindex(NX,NY,NZ))
k = 0
do x = 1,NX
	do y = NBY+1,NY
		do z = 1,NZ
			if (occupancy(x,y,z)%region /= MARROW) cycle
			k = k+1
			vindex(x,y,z) = k
		enddo
	enddo
enddo
NV = k
write(*,*) 'NV: ',NV,4*NV*NV
allocate(V(NV))
allocate(A(NV,NV))
allocate(b(NV))
A = 0
b = 0
influx = 0
y = NBY+1
do x = 1,NX
	do z = 1,NZ
		influx(x,y,z) = RANK_BONE_RATIO*surface(x,z)%signal
	enddo
enddo
row = 0
do x = 1,NX
	if (x > 1) then
		use_nbr(1) = .true.
	else
		use_nbr(1) = .false.
	endif
	if (x < NX) then
		use_nbr(2) = .true.
	else
		use_nbr(2) = .false.
	endif
	do y = NBY+1,NY
		if (y > NBY+1) then
			use_nbr(3) = .true.
		else
			use_nbr(3) = .false.
		endif
		if (y < NY) then
			use_nbr(4) = .true.
		else
			use_nbr(4) = .false.
		endif
		do z = 1,NZ
			if (z > 1) then
				use_nbr(5) = .true.
			else
				use_nbr(5) = .false.
			endif
			if (z < NZ) then
				use_nbr(6) = .true.
			else
				use_nbr(6) = .false.
			endif
			if (occupancy(x,y,z)%region /= MARROW) cycle
			do k = 1,6
				xn = x + neumann(1,k)
				yn = y + neumann(2,k)
				zn = z + neumann(3,k)
				if (vindex(xn,yn,zn) == 0) then
					use_nbr(k) = .false.
				endif
			enddo
			row = row + 1
			col0 = vindex(x,y,z)
			A(row,col0) = A(row,col) + Kdecay*Vol
			b(row) = influx(x,y,z)
			do k = 1,6
				if (use_nbr(k)) then
					offset = neumann(:,k)
					A(row,col0) = A(row,col0) + Kdiffusion*Area
					col = vindex(x+offset(1),y+offset(2),z+offset(3))
					if (col == 0) then
						write(*,*) 'vindex offset error: ',x,y,z,use_nbr,offset
						stop
					endif
					A(row,col) = A(row,col) - Kdiffusion*Area
				endif
			enddo
		enddo
	enddo
enddo

stop
end subroutine
			
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

!----------------------------------------------------------------------------------------
! Solve for steady-state concentration field for a constituent that is defined by:
!	influx(:,:,:)	if > 0, the rate of influx into grid cell (x,y,z)
!					if < 0, (x,y,z) is not MARROW
!	Kdiffusion		diffusion coefficient (um^2.min^-1)
!	Kdecay			decay coefficient (min^-1)
!----------------------------------------------------------------------------------------
subroutine steadystate1(influx,Kdiffusion,Kdecay,C)
real :: influx(:,:,:), C(:,:,:)
real :: Kdiffusion, Kdecay
real :: dx2diff, total, maxchange, dC, sum, dV
real, parameter :: alpha = 0.5
real, parameter :: tol = 1.0e-6		! max change in C at any site as fraction of average C
integer :: x, y, z, xx, yy, zz, nb, nc, k, it
integer :: nt = 10000
real, allocatable :: Ctemp(:,:,:)

dx2diff = DELTA_X**2/Kdiffusion
dV = DELTA_X**3
allocate(Ctemp(NX,NY,NZ))
maxchange = 1.0e10
total = 0
nc = 1
C = 0
do it = 1,nt
	if (mod(it,100) == 0) then
		x = NX/2
		z = NZ/2
!		write(*,'(i6,20e10.3)') it,(C(x,y,z),y=NBY+1,NY/2-NBY,4),maxchange,total/nc
	endif
	if (maxchange < tol*total/nc) then
		write(logmsg,*) 'Convergence reached: it: ',it
		call logger(logmsg)
		exit
	endif
	maxchange = 0
	total = 0
	nc = 0
	do z = 1,NZ
		do y = NBY+1,NY
			do x = 1,NX
				if (influx(x,y,z) < 0) cycle
				if (influx(x,y,z) > 0) then
					dC = influx(x,y,z)*dx2diff
				else
					dC = 0
				endif
				sum = 0
				nb = 0
				do k = 1,6
					xx = x + neumann(1,k)
					yy = y + neumann(2,k)
					zz = z + neumann(3,k)
					if (outside(xx,yy,zz)) cycle
					if (influx(xx,yy,zz) < 0) cycle
					nb = nb + 1
					sum = sum + C(xx,yy,zz)
				enddo
!				Ctemp(x,y,z) = alpha*(sum + dC)/(Kdecay*dx2diff + dC + nb) + (1-alpha)*C(x,y,z)
				Ctemp(x,y,z) = alpha*(DELTA_X*Kdiffusion*sum + influx(x,y,z))/(Kdecay*dV + nb*DELTA_X*Kdiffusion) + (1-alpha)*C(x,y,z)
				dC = abs(Ctemp(x,y,z) - C(x,y,z))
				maxchange = max(dC,maxchange)
				nc = nc + 1
				total = total + Ctemp(x,y,z)
			enddo
		enddo
	enddo
	C = Ctemp
enddo
deallocate(Ctemp)

end subroutine

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
logical function NearOB1(pmono)
type(monocyte_type), pointer :: pmono
real :: v(3), d2, d2lim
integer :: msite(3), ssite(3), isig

NearOB1 = .false.
d2lim = (OB_REACH/DELTA_X)**2
msite = pmono%site
do isig = 1,nsignal
	ssite = signal(isig)%site
	if (surface(ssite(1),ssite(3))%signal == 0) cycle
	v = msite - ssite
	d2 = dot_product(v,v)
	if (d2 < d2lim) then
		NearOB1 = .true.
		exit
	endif
enddo
end function

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine createOsteoclast1(pclump)
type(clump_type), pointer :: pclump
integer :: bonesite(3,100), i, j, npit, ocsite(3)
logical :: inlist
type(occupancy_type), pointer :: p
type(osteoclast_type), pointer :: pclast
real :: tnow
integer :: kpar = 0

tnow = istep*DELTA_T
nclast = nclast + 1
nliveclast = nliveclast + 1
write(logmsg,*) 'createOsteoclast: ',nclast,tnow
call logger(logmsg)
pclast => clast(nclast)
pclast%ID = nclast
!pclast%site = pclump%cm + 0.5
!pclast%cm = pclump%cm
!pclast%site(2) = NBY + 1
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
!		write(*,*) 'createOsteoclast: pit: ',npit,bonesite(:,npit)
	endif
!	pclast%mono(i) = j
	mono(j)%iclast = nclast
	mono(j)%status = OSTEO
enddo
pclast%npit = npit
pclast%count = pclump%ncells
pclast%cm = 0
allocate(pclast%pit(npit))
do i = 1,npit
	pclast%cm = pclast%cm + bonesite(:,i)
enddo
pclast%cm = pclast%cm/npit
pclast%cm(2) = NBY + 0.5
! Make pit location an offset from the cm
do i = 1,npit
	pclast%pit(i)%delta = bonesite(:,i) - (pclast%cm + 0.5)
enddo
call RemoveClump(pclump)
call pitrates(pclast)
write(logmsg,'(a,3f6.1,2i4)') 'clast site, cm, count, npit: ',pclast%cm,pclast%count,pclast%npit
call logger(logmsg)
end subroutine
