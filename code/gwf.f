c===============================================
	subroutine gwf
	1    (c, b, eps, kh, z, rho, u, bf, nz, nc, iz0, gwfrc, tlev)
c===============================================
c...AD parameterizaion with arbitrary tabulated 
c...momentum flux (phase speed) spectrum: rho*u'*v' = b(c)*dc
c...computes force gwfrc in a single azimuthal direction
c
c...gwfrc(z) = 1/rho * dc/cz * B(c(z)) * eps
c...where c(z) is the min phase speed of all waves
c...still propagating at level z
c...(back reflection, or turning points, are accounted for)
c
c...Input:
c...c(1:nc) phase speed (m/s)
c...b(1:nc) mom flux density (kg/m**2/s)
c...kh=horizontal wavenumber (1/m)
c...eps=intermittency factor
c...rho(1:nz)=density profile (kg/m**3)
c...u(1:nz)=wind profile (m/s)
c...bf(1:nz)=bouyancy frequency (1/s)
c...z(1:nz)=altitude (km)
c...iz0=source grid level
c
c...gwfrc(1:nz)=output force
c;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	implicit none

	integer maxz, maxc
	parameter (maxz=150, maxc=300)

c...arguements
	integer nz, nc, iz0
	real c(maxc), b(maxc)
	real z(maxz)
	real rho(maxz)
	real u(maxz)
	real bf(maxz)
	real gwfrc(maxz)
	real eps		!intermittency factor
	real kh
	real tlev		!turbopause level

c...local variables
	real cz1(maxz)		!min c of waves breaking at z
	real cz2(maxz)		!max c of waves breaking at z
	real crefl(maxz)	!phase speed of turning point at z
	real amin, amax		!g(cm), g(cp)
	integer i, j, jj
	real hb(maxz), omc, k2, alp2, cm
	real sat(maxc)		!sat(c) > 0 if unsaturated at z
	real bcz1(maxz)		!b(c1(z))
	real bcz2(maxz)		!b(c2(z))
	real dc1,dc2,dz
	real ds, epsi, const,delc,dc,db
	integer ic1, ic2, izbot

c------------------------------------------------------------------

	do i=1,nz
	  gwfrc(i)=0.
	enddo

	if (nc .eq. 0) return

c...scale height
	do i=2,nz-1
	  hb(i)=-(z(i+1)-z(i-1))/(alog(rho(i+1))-alog(rho(i-1)))
	enddo

	hb(1) =-(z(2)  -z(1))  /(alog(rho(2)) -alog(rho(1)))
	hb(nz)=-(z(nz)-z(nz-1))/(alog(rho(nz))-alog(rho(nz-1)))

	k2=kh*kh
	cm=2000.
	izbot=nz+1			!first breaking level

c...find phase speed at turning point
c...must decrease with altitude
	do i=iz0,nz
	  alp2=1./(4e6*Hb(i)*Hb(i))
	  omc=sqrt(bf(i)*bf(i)*k2/(k2+alp2))
	  crefl(i)=min(u(i)+omc/kh,cm)
	  cm=crefl(i)
	enddo

	do i=iz0,nz
	  const=kh*rho(i)/bf(i)/2.

c...find first and last saturated wave in the spectrum at z(i)
	  ic1=0
	  ic2=0
	  do j=1,nc,1
	    sat(j)=const*(c(j)-u(i))**3-b(j)
	  enddo

	  do j=1,nc
	    if (sat(j) .lt. 0.) ic2=j
	    jj=nc-j+1
	    if (sat(jj) .lt. 0.) ic1=jj
	  enddo

c...interpolate to find c and b(c) where sat(c)=0
 	  if (ic1 .ne. 0) then
	    if (izbot .eq. (nz+1)) izbot=i

	    if (ic1 .eq. 1) then
	      cz1(i)=c(1)
	      bcz1(i)=0.
	    else
	      j=ic1
	      ds=sat(j-1)/(sat(j)-sat(j-1))
	      delc=c(j)-c(j-1)
	      dc=-delc*ds
	      dc=min(dc,delc)
	      dc=max(dc,-delc)
	      db=(b(j)-b(j-1))/delc*dc
	      cz1(i) =c(j-1)+dc
	      bcz1(i)=b(j-1)+db
	    endif

	    if (ic2 .eq. nc) then
	      cz2(i)=c(nc)
	      bcz2(i)=0.
	    else
	      j=ic2+1
	      ds=sat(j-1)/(sat(j)-sat(j-1))
	      delc=c(j)-c(j-1)
	      dc=-delc*ds
	      dc=min(dc,delc)
	      dc=max(dc,-delc)
	      db=(b(j)-b(j-1))/delc*dc
	      cz2(i) =c(j-1)+dc
	      bcz2(i)=b(j-1)+db
	    endif
	    else
c...cz1 > cz2 is a signal no waves break at this level
c...This will allow  c1/c2 to be properly  put in decr/incr order
	      cz1(i)=c(nc)
	      bcz1(i)=0.
	      cz2(i)=c(1)
	      bcz2(i)=0.
	  endif
	enddo

c...remove waves reflected at turning points
c...at or below the current level
	do i=izbot,nz
	  cz1(i)=min(cz1(i),crefl(i))
	  cz2(i)=min(cz2(i),crefl(i))
	enddo

c...c2(z) must increase with altitude, c1(z) decreases
	do i=izbot+1,nz
	  cz1(i)=min(cz1(i),cz1(i-1))
	  cz2(i)=max(cz2(i),cz2(i-1))
	enddo

	epsi=eps*1e-3		!convert from km to m

c...calculate the force
	do i=izbot+1,nz-1
	  dc2=cz2(i+1)-cz2(i-1)
	  dc1=cz1(i+1)-cz1(i-1)
	  dz=z(i+1)-z(i-1)
	  gwfrc(i)=(dc2*bcz2(i)-dc1*bcz1(i))/(dz*rho(i))*epsi
	  if (z(i) .gt. tlev)
	1      gwfrc(i) = gwfrc(i) * exp(-(z(i)-tlev)/3.) 
				!make a turbopause
	enddo

	i=izbot
	dc2=cz2(i+1)-cz2(i)
	dc1=cz1(i+1)-cz1(i)
	dz=z(i+1)-z(i)
	gwfrc(i)=(dc2*bcz2(i)-dc1*bcz1(i))/(dz*rho(i))*epsi

	dc2=cz2(nz-1)-cz2(nz)
	dc1=cz1(nz-1)-cz1(nz)
	dz=z(nz-1)-z(nz)
	gwfrc(nz)=(dc2*bcz2(nz)-dc1*bcz1(nz))/(dz*rho(nz))*epsi

	if (z(nz) .gt. tlev)
	1    gwfrc(nz) = gwfrc(nz) * exp(-(z(nz)-tlev)/3.)

	end
