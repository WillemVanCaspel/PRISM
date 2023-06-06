SUBROUTINE F04AEF(A,IA,B,IB,N,M,C,IC,WKSPCE,AA,IAA,BB,IBB,IFAIL)
!     Simulates Nag routine F04AEF using LINPACK routines, provided B
!     is the identity matrix and N.eq.M. i.e. It inverts a general
!     real square matrix A of order N, and returns the inverse C.
!     WKSPCE should be an array of dimension at least N. AA and BB
!     are not used. Unless the routine detects an error, IFAIL
!     contains 0 on exit.
!
        integer IA, IB, N, M, IC, IAA, IBB, IFAIL
        REAL*8 A(IA,N), B(IB,M), C(IC,M), WKSPCE(N)
        REAL*8 AA(IAA,N), BB(IBB,N)
!
        integer inmax
        parameter (inmax=300)
        integer i, j, job, ipvt(inmax)
        REAL*8 rcond, detmnt(2)
!
        if (N.gt.inmax) stop 'F04AEF'
!     Increase value of inmax if above statement is executed
!
!
        if (N.ne.M) then
        write(6,*) 'N ne M in f04aef'
        IFAIL=1
        return
        endif
!     
        do 10 j=1,N
        do 10 i=1,N
    10     C(i,j)= A(i,j)
!
        call sgeco(C, IC, N, ipvt, rcond, WKSPCE)
        if (1.0+rcond.eq.1.0) then
        IFAIL=1
        write(6,*) 'matrix ill conditioned in f04aef'
        return
        endif
!
        job= 01
        call sgedi(C, IC, N, ipvt, detmnt, WKSPCE, job)
!
        IFAIL=0
        return
        end

        SUBROUTINE F03AAF(A,IA,N,DET,WKSPCE,IFAIL)
!     Simulates Nag routine F03AAF using LINPACK routines. i.e. It
!     computes the determinant DET of a general real square matrix
!     A of order N. WKSPCE should be an array of dimension at least
!     N. Unless the routine detects an error, IFAIL contains 0 on
!     exit.
! **Note** This routine destroys the matrix A
!
        integer IA, N, IFAIL
        REAL*8 A(IA,N), DET, WKSPCE(N)
!
        integer inmax
        parameter (inmax=100)
        integer job, ipvt(inmax)
        REAL*8 rcond, detmnt(2)
!
        if (N.gt.inmax) stop 'F03AAF'
!     Increase value of inmax if above statement is executed
!
        call sgeco(A, IA, N, ipvt, rcond, WKSPCE)
        if (1.0+rcond.eq.1.0) then
        IFAIL=1
        return
        endif
!
        job= 10
        call sgedi(A, IA, N, ipvt, detmnt, WKSPCE, job)
!
        DET= detmnt(1) * 10.0**detmnt(2)
        IFAIL=0
        return
        end
!
! *** The following routines are from the LINPACK library
!
        integer function isamax(n,sx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!
        REAL*8 sx(1),smax
        integer i,incx,ix,n
!
        isamax = 0
        if( n .lt. 1 ) return
        isamax = 1
        if(n.eq.1)return
        if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
        ix = 1
        smax = abs(sx(1))
        ix = ix + incx
        do 10 i = 2,n
            if(abs(sx(ix)).le.smax) go to 5
            isamax = i
            smax = abs(sx(ix))
    5    ix = ix + incx
    10 continue
        return
!
!        code for increment equal to 1
!
    20 smax = abs(sx(1))
        do 30 i = 2,n
            if(abs(sx(i)).le.smax) go to 30
            isamax = i
            smax = abs(sx(i))
    30 continue
        return
        end
        REAL*8 function sasum(n,x,incx)
!
!     takes the sum of the absolute values.
!     jack dongarra, linpack, 3/11/78.
!
        REAL*8 x(1),temp
        integer i,incx,m,mp1,n,nincx
!
        sasum = 0.0e0
        temp = 0.0e0
        if(n.le.0)return
        if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
        nincx = n*incx
        do 10 i = 1,nincx,incx
        temp = temp + abs(x(i))
    10 continue
        sasum = temp
        return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
    20 m = mod(n,6)
        if( m .eq. 0 ) go to 40
        do 30 i = 1,m
        temp = temp + abs(x(i))
    30 continue
        if( n .lt. 6 ) go to 60
    40 mp1 = m + 1
        do 50 i = mp1,n,6
        temp = temp + abs(x(i)) + abs(x(i + 1)) + abs(x(i + 2)) &
          + abs(x(i + 3)) + abs(x(i + 4)) + abs(x(i + 5))
    50 continue
    60 sasum = temp
        return
        end
        subroutine saxpy(n,sa,sx,incx,sy,incy)
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
        REAL*8 sx(1),sy(1),sa
        integer i,incx,incy,ix,iy,m,mp1,n
!
        if(n.le.0)return
        if (sa .eq. 0.0e0) return
        if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
        ix = 1
        iy = 1
        if(incx.lt.0)ix = (-n+1)*incx + 1
        if(incy.lt.0)iy = (-n+1)*incy + 1
        do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
    10 continue
        return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
    20 m = mod(n,4)
        if( m .eq. 0 ) go to 40
        do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
    30 continue
        if( n .lt. 4 ) return
    40 mp1 = m + 1
        do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
    50 continue
        return
        end
        REAL*8 function sdot(n,sx,incx,sy,incy)
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!
        REAL*8 sx(1),sy(1),temp
        integer i,incx,incy,ix,iy,m,mp1,n
!
        sdot = 0.0e0
        temp = 0.0e0
        if(n.le.0)return
        if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
        ix = 1
        iy = 1
        if(incx.lt.0)ix = (-n+1)*incx + 1
        if(incy.lt.0)iy = (-n+1)*incy + 1
        do 10 i = 1,n
        temp = temp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
    10 continue
        sdot = temp
        return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
    20 m = mod(n,5)
        if( m .eq. 0 ) go to 40
        do 30 i = 1,m
        temp = temp + sx(i)*sy(i)
    30 continue
        if( n .lt. 5 ) go to 60
    40 mp1 = m + 1
        do 50 i = mp1,n,5
        temp = temp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) + &
            sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
    50 continue
    60 sdot = temp
        return
        end
        subroutine sgeco(a,lda,n,ipvt,rcond,z)
        integer lda,n,ipvt(1)
        REAL*8 a(lda,1),z(1)
        REAL*8 rcond
!
!     sgeco factors a real matrix by gaussian elimination
!     and estimates the condition of the matrix.
!
!     if  rcond  is not needed, sgefa is slightly faster.
!     to solve  a*x = b , follow sgeco by sgesl.
!     to compute  inverse(a)*c , follow sgeco by sgesl.
!     to compute  determinant(a) , follow sgeco by sgedi.
!     to compute  inverse(a) , follow sgeco by sgedi.
!
!     on entry
!
!        a       REAL*8(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        rcond   REAL*8
!                an estimate of the reciprocal condition of  a .
!                for the system  a*x = b , relative perturbations
!                in  a  and  b  of size  epsilon  may cause
!                relative perturbations in  x  of size  epsilon/rcond .
!                if  rcond  is so small that the logical expression
!                           1.0 + rcond .eq. 1.0
!                is true, then  a  may be singular to working
!                precision.  in particular,  rcond  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        z       REAL*8(n)
!                a work vector whose contents are usually unimportant.
!                if  a  is close to a singular matrix, then  z  is
!                an approximate null vector in the sense that
!                norm(a*z) = rcond*norm(a)*norm(z) .
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     linpack sgefa
!     blas saxpy,sdot,sscal,sasum
!     fortran abs,max,sign
!
!     internal variables
!
        REAL*8 sdot,ek,t,wk,wkm
        REAL*8 anorm,s,sasum,sm,ynorm
        integer info,j,k,kb,kp1,l
!
!
!     compute 1-norm of a
!
        anorm = 0.0e0
        do 10 j = 1, n
            anorm = max(anorm,sasum(n,a(1,j),1))
    10 continue
!
!     factor
!
        call sgefa(a,lda,n,ipvt,info)
!
!     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
!     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
!     trans(a)  is the transpose of a .  the components of  e  are
!     chosen to cause maximum local growth in the elements of w  where
!     trans(u)*w = e .  the vectors are frequently rescaled to avoid
!     overflow.
!
!     solve trans(u)*w = e
!
        ek = 1.0e0
        do 20 j = 1, n
            z(j) = 0.0e0
    20 continue
        do 100 k = 1, n
            if (z(k) .ne. 0.0e0) ek = sign(ek,-z(k))
            if (abs(ek-z(k)) .le. abs(a(k,k))) go to 30
            s = abs(a(k,k))/abs(ek-z(k))
            call sscal(n,s,z,1)
            ek = s*ek
    30    continue
            wk = ek - z(k)
            wkm = -ek - z(k)
            s = abs(wk)
            sm = abs(wkm)
            if (a(k,k) .eq. 0.0e0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
            go to 50
    40    continue
            wk = 1.0e0
            wkm = 1.0e0
    50    continue
            kp1 = k + 1
            if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
                sm = sm + abs(z(j)+wkm*a(k,j))
                z(j) = z(j) + wk*a(k,j)
                s = s + abs(z(j))
    60       continue
            if (s .ge. sm) go to 80
                t = wkm - wk
                wk = wkm
                do 70 j = kp1, n
                    z(j) = z(j) + t*a(k,j)
    70          continue
    80       continue
    90    continue
            z(k) = wk
    100 continue
        s = 1.0e0/sasum(n,z,1)
        call sscal(n,s,z,1)
!
!     solve trans(l)*y = w
!
        do 120 kb = 1, n
            k = n + 1 - kb
            if (k .lt. n) z(k) = z(k) + sdot(n-k,a(k+1,k),1,z(k+1),1)
            if (abs(z(k)) .le. 1.0e0) go to 110
            s = 1.0e0/abs(z(k))
            call sscal(n,s,z,1)
    110    continue
            l = ipvt(k)
            t = z(l)
            z(l) = z(k)
            z(k) = t
    120 continue
        s = 1.0e0/sasum(n,z,1)
        call sscal(n,s,z,1)
!
        ynorm = 1.0e0
!
!     solve l*v = y
!
        do 140 k = 1, n
            l = ipvt(k)
            t = z(l)
            z(l) = z(k)
            z(k) = t
            if (k .lt. n) call saxpy(n-k,t,a(k+1,k),1,z(k+1),1)
            if (abs(z(k)) .le. 1.0e0) go to 130
            s = 1.0e0/abs(z(k))
            call sscal(n,s,z,1)
            ynorm = s*ynorm
    130    continue
    140 continue
        s = 1.0e0/sasum(n,z,1)
        call sscal(n,s,z,1)
        ynorm = s*ynorm
!
!     solve  u*z = v
!
        do 160 kb = 1, n
            k = n + 1 - kb
            if (abs(z(k)) .le. abs(a(k,k))) go to 150
            s = abs(a(k,k))/abs(z(k))
            call sscal(n,s,z,1)
            ynorm = s*ynorm
    150    continue
            if (a(k,k) .ne. 0.0e0) z(k) = z(k)/a(k,k)
            if (a(k,k) .eq. 0.0e0) z(k) = 1.0e0
            t = -z(k)
            call saxpy(k-1,t,a(1,k),1,z(1),1)
    160 continue
!     make znorm = 1.0
        s = 1.0e0/sasum(n,z,1)
        call sscal(n,s,z,1)
        ynorm = s*ynorm
!
        if (anorm .ne. 0.0e0) rcond = ynorm/anorm
        if (anorm .eq. 0.0e0) rcond = 0.0e0
        return
        end
        subroutine sgedi(a,lda,n,ipvt,det,work,job)
        integer lda,n,ipvt(1),job
        REAL*8 a(lda,1),det(2),work(1)
!
!     sgedi computes the determinant and inverse of a matrix
!     using the factors computed by sgeco or sgefa.
!
!     on entry
!
!        a       REAL*8(lda, n)
!                the output from sgeco or sgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from sgeco or sgefa.
!
!        work    REAL*8(n)
!                work vector.  contents destroyed.
!
!        job     integer
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     on return
!
!        a       inverse of original matrix if requested.
!                otherwise unchanged.
!
!        det     REAL*8(2)
!                determinant of original matrix if requested.
!                otherwise not referenced.
!                determinant = det(1) * 10.0**det(2)
!                with  1.0 .le. abs(det(1)) .lt. 10.0
!                or  det(1) .eq. 0.0 .
!
!     error condition
!
!        a division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        it will not occur if the subroutines are called correctly
!        and if sgeco has set rcond .gt. 0.0 or sgefa has set
!        info .eq. 0 .
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas saxpy,sscal,sswap
!     fortran abs,mod
!
!     internal variables
!
        REAL*8 t
        REAL*8 ten
        integer i,j,k,kb,kp1,l,nm1
!
!
!     compute determinant
!
        if (job/10 .eq. 0) go to 70
            det(1) = 1.0e0
            det(2) = 0.0e0
            ten = 10.0e0
            do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
!        ...exit
            if (det(1) .eq. 0.0e0) go to 60
    10       if (abs(det(1)) .ge. 1.0e0) go to 20
                det(1) = ten*det(1)
                det(2) = det(2) - 1.0e0
            go to 10
    20       continue
    30       if (abs(det(1)) .lt. ten) go to 40
                det(1) = det(1)/ten
                det(2) = det(2) + 1.0e0
            go to 30
    40       continue
    50    continue
    60    continue
    70 continue
!
!     compute inverse(u)
!
        if (mod(job,10) .eq. 0) go to 150
            do 100 k = 1, n
            a(k,k) = 1.0e0/a(k,k)
            t = -a(k,k)
            call sscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
                t = a(k,j)
                a(k,j) = 0.0e0
                call saxpy(k,t,a(1,k),1,a(1,j),1)
    80       continue
    90       continue
    100    continue
!
!        form inverse(u)*inverse(l)
!
            nm1 = n - 1
            if (nm1 .lt. 1) go to 140
            do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
                work(i) = a(i,k)
                a(i,k) = 0.0e0
    110       continue
            do 120 j = kp1, n
                t = work(j)
                call saxpy(n,t,a(1,j),1,a(1,k),1)
    120       continue
            l = ipvt(k)
            if (l .ne. k) call sswap(n,a(1,k),1,a(1,l),1)
    130    continue
    140    continue
    150 continue
        return
        end
        subroutine sgefa(a,lda,n,ipvt,info)
        integer lda,n,ipvt(1),info
        REAL*8 a(lda,1)
!
!     sgefa factors a real matrix by gaussian elimination.
!
!     sgefa is usually called by sgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
!
!     on entry
!
!        a       REAL*8(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that sgesl or sgedi will divide by zero
!                     if called.  use  rcond  in sgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas saxpy,sscal,isamax
!
!     internal variables
!
        REAL*8 t
        integer isamax,j,k,kp1,l,nm1
!
!
!     gaussian elimination with partial pivoting
!
        info = 0
        nm1 = n - 1
        if (nm1 .lt. 1) go to 70
        do 60 k = 1, nm1
            kp1 = k + 1
!
!        find l = pivot index
!
            l = isamax(n-k+1,a(k,k),1) + k - 1
            ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
            if (a(l,k) .eq. 0.0e0) go to 40
!
!           interchange if necessary
!
            if (l .eq. k) go to 10
                t = a(l,k)
                a(l,k) = a(k,k)
                a(k,k) = t
    10       continue
!
!           compute multipliers
!
            t = -1.0e0/a(k,k)
            call sscal(n-k,t,a(k+1,k),1)
!
!           row elimination with column indexing
!
            do 30 j = kp1, n
                t = a(l,j)
                if (l .eq. k) go to 20
                    a(l,j) = a(k,j)
                    a(k,j) = t
    20          continue
                call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
    30       continue
            go to 50
    40    continue
            info = k
    50    continue
    60 continue
    70 continue
        ipvt(n) = n
        if (a(n,n) .eq. 0.0e0) info = n
        return
        end
        subroutine  sscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!
        REAL*8 da,dx(1)
        integer i,incx,m,mp1,n,nincx
!
        if(n.le.0)return
        if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
        nincx = n*incx
        do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
    10 continue
        return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
    20 m = mod(n,5)
        if( m .eq. 0 ) go to 40
        do 30 i = 1,m
        dx(i) = da*dx(i)
    30 continue
        if( n .lt. 5 ) return
    40 mp1 = m + 1
        do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
    50 continue
        return
        end
        subroutine sswap (n,sx,incx,sy,incy)
!
!     interchanges two vectors.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.
!
        REAL*8 sx(1),sy(1),stemp
        integer i,incx,incy,ix,iy,m,mp1,n
!
        if(n.le.0)return
        if(incx.eq.1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
        ix = 1
        iy = 1
        if(incx.lt.0)ix = (-n+1)*incx + 1
        if(incy.lt.0)iy = (-n+1)*incy + 1
        do 10 i = 1,n
        stemp = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp
        ix = ix + incx
        iy = iy + incy
    10 continue
        return
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
    20 m = mod(n,3)
        if( m .eq. 0 ) go to 40
        do 30 i = 1,m
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
    30 continue
        if( n .lt. 3 ) return
    40 mp1 = m + 1
        do 50 i = mp1,n,3
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
        stemp = sx(i + 1)
        sx(i + 1) = sy(i + 1)
        sy(i + 1) = stemp
        stemp = sx(i + 2)
        sx(i + 2) = sy(i + 2)
        sy(i + 2) = stemp
    50 continue
        return
        end
    