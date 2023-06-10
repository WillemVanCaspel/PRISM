      subroutine pwl_interp_2d ( nxd, nyd, xd, yd, zd, ni, xi, yi, zi )

c*********************************************************************72
c
cc PWL_INTERP_2D: piecewise linear interpolant to data defined on a 2D grid.
c
c  Discussion:
c
c    Thanks to Adam Hirst for pointing out an error in the formula that
c    chooses the interpolation triangle, 04 February 2018.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 February 2018
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NXD, NYD, the number of X and Y data values.
c
c    Input, double precision XD(NXD), YD(NYD), the sorted X and Y data.
c
c    Input, double precision ZD(NXD,NYD), the Z data.
c
c    Input, integer NI, the number of interpolation points.
c
c    Input, double precision XI(NI), YI(NI), the coordinates of the
c    interpolation points.
c
c    Output, double precision ZI(NI), the value of the interpolant.
c
      implicit none

      integer ni
      integer nxd
      integer nyd

      double precision alpha
      double precision beta
      double precision det
      double precision dxa
      double precision dxb
      double precision dxi
      double precision dya
      double precision dyb
      double precision dyi
      double precision gamma
      integer i
      integer j
      integer k
      double precision r8_huge
      integer r8vec_bracket5
      double precision xd(nxd)
      double precision xi(ni)
      double precision yd(nyd)
      double precision yi(ni)
      double precision zd(nxd,nyd)
      double precision zi(ni)

      do k = 1, ni
c
c  For interpolation point (xi(k),yi(k)), find data intervals I and J so that:
c
c    xd(i) <= xi(k) <= xd(i+1),
c    yd(j) <= yi(k) <= yd(j+1).
c
c  But if the interpolation point is not within a data interval, 
c  assign the dummy interpolant value zi(k) = infinity.
c
        i = r8vec_bracket5 ( nxd, xd, xi(k) )
        if ( i .eq. -1 ) then
          write(6,*) 'interpolation error'
          zi(k) = r8_huge ( )
          cycle
        end if

        j = r8vec_bracket5 ( nyd, yd, yi(k) )
        if ( j .eq. -1 ) then
          write(6,*) 'interpolation error'
          zi(k) = r8_huge ( )
          cycle
        end if
c
c  The rectangular cell is arbitrarily split into two triangles.
c  The linear interpolation formula depends on which triangle 
c  contains the data point.
c
c    (I,J+1)--(I+1,J+1)
c      |\       |
c      | \      |
c      |  \     |
c      |   \    |
c      |    \   |
c      |     \  |
c    (I,J)---(I+1,J)
c
        if ( yi(k) .lt. yd(j+1) + ( yd(j) - yd(j+1) ) 
     &    * ( xi(k) - xd(i) ) / ( xd(i+1) - xd(i) ) ) then

          dxa = xd(i+1) - xd(i)
          dya = yd(j)   - yd(j)

          dxb = xd(i)   - xd(i)
          dyb = yd(j+1) - yd(j)

          dxi = xi(k)   - xd(i)
          dyi = yi(k)   - yd(j)

          det = dxa * dyb - dya * dxb

          alpha = ( dxi * dyb - dyi * dxb ) / det
          beta =  ( dxa * dyi - dya * dxi ) / det
          gamma = 1.0D+00 - alpha - beta

          zi(k) = alpha * zd(i+1,j) + beta * zd(i,j+1) + gamma * zd(i,j)

        else

          dxa = xd(i)   - xd(i+1)
          dya = yd(j+1) - yd(j+1)

          dxb = xd(i+1) - xd(i+1)
          dyb = yd(j)   - yd(j+1)

          dxi = xi(k)   - xd(i+1)
          dyi = yi(k)   - yd(j+1)

          det = dxa * dyb - dya * dxb

          alpha = ( dxi * dyb - dyi * dxb ) / det
          beta =  ( dxa * dyi - dya * dxi ) / det
          gamma = 1.0D+00 - alpha - beta

          zi(k) = alpha * zd(i,j+1) + beta * zd(i+1,j) 
     &      + gamma * zd(i+1,j+1)

        end if

      end do

      return
      end
