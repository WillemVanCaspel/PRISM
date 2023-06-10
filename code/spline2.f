CD ---------------------------------------------------------------------------
      SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,MD,ND,Y2A,IER)
CD ---------------------------------------------------------------------------
CD this is a slightly different version to support init_wind_temp.for

CD    SUBROUTINE SPLIE2
CD
CD    Purpose:
CD      constructs one-dimensional natural cubic splines and returns 2nd derivatives
CD
CD    Author: W. H. Press, B. P. Flannery, S. A. Teukolsky, and
CD            W. T. Vetterling
CD
CD    Modifier:
CD      S. L. Osburn
CD
CD    Usage:
CD      CALL SPLIE2(X1A,X2A,YA,M,N,Y2A,IER)
CD
CD   Arguments:
CD
CD  Input
CD Parameters   Type    Units   Dimensions               Description
CD   M          Integer   -     none                     dimensions
CD   N          Integer   -     none                     
CD   X1A        Real      -     M                        data points
CD   X2A        Real      -     N                         (x1,x2,y)
CD   YA         Real      -     M,N                      
CD
CD  Output
CD Parameters   Type    Units   Dimensions               Description
CD   Y2A        Real      -     M,N                      second derivatives at each point
CD   IER        Integer   -     none                   error code
CD                                                     = 0, no error
CD
CD
CD
CD  Example:
CD
CD  Restrictions:
CD
CD  Subroutine and functions required:
CD    SPLINE
CD
CD  INCLUDE code files required:
CD    none
CD
CD  Method:
CD    See "Numerical Recipies", 1986 edition, section 3.6,
CD  "Interpolation in Two or More Dimensions", page 95.
CD
CD  References:
CD    Code is a slightly modified version of routine SPLIE2 found
CD in "Numerical Recipes, The Art of Scientific Computing", 
CD by William H. Press, Brian P. Flannery, Saul A. Teukolsky, and
CD William T. Vetterling.  Published by Cambridge University Press, 1986.
CD Code found on pages 100-101.
CD    Additional references are:
CD
CD  Change History:
CD   ver        date         by                        description
CD   1.0 September 3, 1987  SLO                      documented
CD
CD---------------------------------------------------------------------------
CD
c
c    NN...max expected value for N or M
c
      PARAMETER (NN=500)
      DIMENSION X1A(M),X2A(N),YA(MD,ND),Y2A(MD,ND),YTMP(NN),Y2TMP(NN)
      IER = 0
      DO 13 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
11      CONTINUE
c
c    values 1E+30 signal a natural spline
        CALL SPLINE(X2A,YTMP,N,1.E30,1.E30,Y2TMP,IER)
        DO 12 K=1,N
          Y2A(J,K)=Y2TMP(K)
12      CONTINUE
13    CONTINUE
      RETURN
      END

CD---------------------------------------------------------------------------
      SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,MD,ND,X1,X2,Y,IER)
CD
CD---------------------------------------------------------------------------
CD
CD    SUBROUTINE SPLIN2
CD
CD    Purpose:
CD      return an interpolated Y value using bicubic interpolation for X1 and X2 
CD
CD    Author: W. H. Press, B. P. Flannery, S. A. Teukolsky, and
CD            W. T. Vetterling
CD
CD    Modifier:
CD      S. L. Osburn
CD
CD    Usage:
CD      CALL SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y,IER)
CD
CD   Arguments:
CD
CD  Input
CD Parameters   Type    Units   Dimensions               Description
CD   M          Integer   -     none                     number of coordinates in X1
CD   N          Integer   -     none                     number of coordinates in X2
CD   X1         Real      -     none                     coordinate to
CD   X2         Real      -     none                       interpolate
CD   X1A        Real      -     M                        data points in form:
CD   X2A        Real      -     N                         ( X1,X2,Y )
CD   YA         Real      -     M,N                      
CD
CD  Output
CD Parameters   Type    Units   Dimensions               Description
CD   Y          Real      -     none                     interpolated result
CD   Y2A        Real      -     M,N                      second derivatives at each data point 
CD   IER        Integer   -     none                   error code
CD                                                     = 0, no error
CD                                                     ERROR IN SPLINT
CD                                                     = 10 BAD XA INPUT
CD
CD
CD
CD  Example:
CD
CD  Restrictions:
CD
CD  Subroutine and functions required:
CD    SPLINE
CD    SPLINT
CD
CD  INCLUDE code files required:
CD    none
CD
CD  Method:
CD    See "Numerical Recipies", 1986 edition, section 3.6,
CD  "Interpolation in Two or More Dimensions", page 95.
CD
CD  References:
CD    Code is a slightly modified version of routine SPLIN2 found
CD in "Numerical Recipes, The Art of Scientific Computing", 
CD by William H. Press, Brian P. Flannery, Saul A. Teukolsky, and 
CD William T. Vetterling.  Published by Cambridge University Press, 1986.
CD Code found on page 101.
CD    Additional references are:
CD
CD  Change History:
CD   ver        date         by                        description
CD   1.0 September 3, 1987  SLO                      documented
CD
CD---------------------------------------------------------------------------
CD
c
c    NN...maximum expected dimension of N or M
c
      PARAMETER (NN=500)
      DIMENSION X1A(M),X2A(N),YA(MD,ND),Y2A(MD,ND),YTMP(NN),
     $     Y2TMP(NN),YYTMP(NN)
      IER = 0
c
c    do M evaluations of row splines constructed with SPLIE2
c       using 1-dimensional spline evaluator SPLINT
c
      DO 12 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
          Y2TMP(K)=Y2A(J,K)
11      CONTINUE
        CALL SPLINT(X2A,YTMP,Y2TMP,N,X2,YYTMP(J),IER)
        IER = 10 * IER
        IF (IER.EQ.10) RETURN
12    CONTINUE
c
c    construct a one dimensional column spline
c
      CALL SPLINE(X1A,YYTMP,M,1.E30,1.E30,Y2TMP,IER)
c                           
c    evaluate it
c
      CALL SPLINT(X1A,YYTMP,Y2TMP,M,X1,Y,IER)
      IER = IER * 10
      RETURN
      END
CD---------------------------------------------------------------------------
      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2,IER)
CD
CD---------------------------------------------------------------------------
CD
CD    SUBROUTINE SPLINE
CD
CD    Purpose:
CD      given a set of points cooresponding to an interpolated function, and 
CD    its first derivative at each end of the function, returns the second
CD    derivative of the interpolated function at each point
CD
CD    Author: W. H. Press, B. P. Flannery, S. A. Teukolsky, and
CD            W. T. Vetterling
CD
CD    Modifier:
CD      S. L. Osburn
CD
CD    Usage:
CD      CALL SPLINE(X,Y,N,YP1,YPN,Y2,IER)
CD
CD   Arguments:
CD
CD  Input
CD Parameters   Type    Units   Dimensions               Description
CD   N          Integer   -     none                     number of data points
CD   X          Real      -     N                        data points in
CD   Y          Real      -     N                           ascending order of X
CD   YP1        Real      -     none                     1st derivatives at X(1)
CD   YPN        Real      -     none                       and X(N)
CD
CD  Output
CD Parameters   Type    Units   Dimensions               Description
CD   Y2         Real      -     N                        2nd derivatives at each point
CD   IER        Integer   -     none                   error code
CD                                                     = 0, no error
CD
CD
CD
CD  Example:
CD
CD  Restrictions:
CD
CD  Subroutine and functions required:
CD    none
CD
CD  INCLUDE code files required:
CD    none
CD
CD  Method:
CD    See "Numerical Recipies", 1986 edition, section 3.3,
CD  "Cubic Spline Interpolation", page 86.
CD
CD  References:
CD    Code is a slightly modified version of routine SPLINE found
CD in "Numerical Recipes, The Art of Scientific Computing", 
CD by William H. Press, Brian P. Flannery, Saul A. Teukolsky, and 
CD William T. Vetterling.  Published by Cambridge University Press, 1986.
CD Code found on page 88.
CD    Additional references are:
CD
CD  Change History:
CD   ver        date         by                        description
CD   1.0 September 3, 1987  SLO                      documented
CD
CD---------------------------------------------------------------------------
CD
c
c    NMAX...maximum expected N
c
      PARAMETER (NMAX=500)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)
      IER = 0
c
c    lower boundary condition set to be natural or...
c
      IF (YP1.GT..99E30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
c
c    ...have a specified 1st derivative
c
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
c
c    decomposition loop of tridiagonal algorithm (Y2 and U are used as temp. storage)
c
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
c
c    upper boundary condition set to be natural or...
c
      IF (YPN.GT..99E30) THEN
        QN=0.
        UN=0.
      ELSE
c
c    ...have a specified first derivative
c
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
c
c    backsubstitution loop of tridiagonal algorithm
c
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END
CD---------------------------------------------------------------------------
      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y,IER)
CD
CD---------------------------------------------------------------------------
CD
CD    SUBROUTINE SPLINT
CD
CD    Purpose:
CD      preform a cubic-spline interpolation
CD
CD    Author: W. H. Press, B. P. Flannery, S. A. Teukolsky, and
CD            W. T. Vetterling
CD
CD    Modifier:
CD      S. L. Osburn
CD
CD    Usage:
CD      CALL SPLINT(XA,YA,Y2A,N,X,Y,IER)
CD
CD   Arguments:
CD
CD  Input
CD Parameters   Type    Units   Dimensions               Description
CD   N          Integer   -     none                     number of points
CD   X          Real      -     none                     point to interpolate
CD   XA         Real      -     N                        data points
CD   YA         Real      -     N                        
CD   Y2A        Real      -     N                        array output from SPLINE
CD
CD  Output
CD Parameters   Type    Units   Dimensions               Description
CD   Y          Real      -     none                     interpolation
CD   IER        Integer   -     none                     error code
CD                                                       = 0, no error
CD                                                       = 1, bad XA input
CD
CD
CD
CD  Example:
CD
CD  Restrictions:
CD    The code has been slightly modified to return an error code
CD  instead of exicuting a PAUSE if an unusual circumstance is
CD  encountered.  Error arguments must be included in subroutine call.
CD
CD  Subroutine and functions required:
CD    none
CD
CD  INCLUDE code files required:
CD    none
CD
CD  Method:
CD    See "Numerical Recipies", 1986 edition, section 3.3,
CD  "Cubic Spline Interpolation", page 86.
CD
CD  References:
CD    Code is a slightly modified version of routine SPLINT found
CD in "Numerical Recipes, The Art of Scientific Computing", 
CD by William H. Press, Brian P. Flannery, Saul A. Teukolsky, and 
CD William T. Vetterling.  Published by Cambridge University Press, 1986.
CD Code found on page 89.
CD    Additional references are:
CD
CD  Change History:
CD   ver        date         by                        description
CD   1.0 August 28, 1987    SLO                      documented and added error
CD                                                   return
CD
CD---------------------------------------------------------------------------
CD
      DIMENSION XA(N),YA(N),Y2A(N)
      IER = 0
c
c    find place by bisection (good random calls of points)
c
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
c
c    KLO and KHI bracket the point
c
      H=XA(KHI)-XA(KLO)
c
c    bad XA
c
      IF (H.EQ.0.) THEN
         IER = 1
         RETURN
      ENDIF
c
c    cubic spline polynomial
c
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
c
c    and the value is
c
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END

