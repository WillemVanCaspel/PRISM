C SPFFTB.INC
C
C Common module for f.f.t. table used by sptransforms
      INTEGER NFFTAB
      REAL FFTAB
      PARAMETER (NFFTAB=512)
      COMMON/SPFFTB/ FFTAB(4,2,NFFTAB)
