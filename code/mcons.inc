C MCONS.INC
C  This file is included in the routines of PROGNOS.FOR and should also be
C included in a main driver.
C
C Parameters (array dimensions):
C L1MAX    -- maximum no. of pressure-levels (>=2)
C NPGQ	   -- total # tracers (>=3) (VOR, DIV, POT, passive)
C NTRACE   -- number of "passive" tracers (<= NPGQ-3)
C             (Default value is 0). Set in main program.
C
C Array pointers:
C Prognostic quantities
C JVOR   -- denotes relative vorticity
C JDIV   -- denotes horizontal divergence
C JPOT   -- denotes potential temperature
C JTR1   -- denotes tracer 1
C
C Vector components
C JX     -- denotes zonal direction
C JY     -- denotes meridional direction

	INTEGER L1MAX
	PARAMETER (L1MAX=201)

	INTEGER NPGQ
	PARAMETER (NPGQ=3)

	INTEGER NTRACE
	COMMON /NTR/ NTRACE

	INTEGER JVOR, JDIV, JPOT, JTR1, JX, JY
	PARAMETER (JVOR=1, JDIV=2, JPOT=3, JTR1=4, JX=1, JY=2)
