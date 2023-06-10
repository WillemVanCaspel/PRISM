C...heat code has 2 km grid. some RADPARAM values hardwired here

C...STUFF RELATED TO RADIATIVE HEATING CODE SETUP
	INTEGER KM
	PARAMETER (KM=76)

	REAL*8 YKH2O(10)/1.0D-3,1.33D-2,4.22D-2,1.334D-1,4.217D-1  ! k-coefficients
	1    ,1.334D0,5.623D0,31.62D0,177.8D0,1.0D3/               ! for water vapor
	REAL*8 DGKH2O(10)/3.0015D-1,5.019D-2,4.56D-2,3.827D-2      ! heating rate in 
	1    ,2.968D-2,2.284D-2,2.324D-2,1.232D-2,5.16D-3,2.33D-3/ ! the troposphere
	COMMON /KDGH2O/ YKH2O, DGKH2O

c...f107 =(70,145,220): 60 < F107 < 270 ==> 10.7 cm solar flux index
	DOUBLE PRECISION F107 /70./
	DOUBLE PRECISION theta27 /0./! = phase of 27-day periodic variation: sin(theta27)
	DOUBLE PRECISION OME0 /.25/	!EFFECTIVE ALBEDO
	DOUBLE PRECISION DEC	!SOLAR DECLINATION

c...mdlrd=1, no heat; 2 => diurnal avg heat rate; 3 => local heat rate;
	INTEGER MDLRD		
	INTEGER KLBD /6/		! effective cloud-top level
	DOUBLE PRECISION XMU/0.5D0/	! COS(solar zenith angle) = DUMMY VARIABLE
	INTEGER LTMU/1/		! 1 => input=DEC, PHI, TIME(hr); 0 => direct input = XMU

	COMMON /RADPARAM/ F107, THETA27, OME0, DEC, MDLRD, KLBD
