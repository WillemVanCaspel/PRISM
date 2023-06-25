module Config_mod

  implicit none

  integer, parameter :: config_IO = 42, TXTLEN_FILE = 64  

  character(len=TXTLEN_FILE), public, save :: &
    ! paths to folders containing input and output files
     INPATH     = 'NOTSET' &
    ,OUTPATH    = 'NOTSET' & ! can be changed to local folder later on

    ! file names for tidal heating, background atmosphere, and ocean tides in INPATH
    ,TIDEHEAT   = 'NOTSET' &
    ,ATMOSPHERE = 'NOTSET' &
    ,OCEANTIDE  = 'NOTSET' &

    ! input files for 1-D pressure and geom. altitude grids,
    ! used to initialize sigma-coords and dissipation profs. 
    ,GRIDP      = 'NOTSET' & 
    ,GRIDZ      = 'NOTSET' &

    ! 2-D temperature and zonal wind files for initialization
    ,WINDINIT   = 'NOTSET' &
    ,TEMPINIT   = 'NOTSET' 
  
  logical, public, save :: &   
    ! flags controlling final output dump, and start from continuation file
     write_continuation     = .false. &
    ,start_continuation     = .false.

  real, public, save ::  &
    ! vertical grid parameters
     zbot_conf = 170. & ! so and so
    ,zlid_conf = 170. & ! 
    ,dz0_conf  = 0.   &
    ,delz_conf = 0.   &
    ,zlow_conf = 0.     
    
  real, public, save :: &
    ! time resolution and run time config
     model_dt = 0.2   & ! time step in hours
    ,runtime  = 161.  & ! total run time in days
    ,dt_out   = 1.0   & ! time resolution output (hrs)
    ,start_t  = 0.      ! start time (doy)

  ! meridional legendre polynom. (ntrunc) and zonal wave nr. (mtrunc) spectral resolution
  integer, public, save :: ntrunc = 3, mtrunc = 16

  ! to-do, get rid of ray-fac

  ! make stratosphere and mesosphere limits pressure-dependent

     

  public :: Config_Constants

  ! -------------------------------------------------------------------------- !

  contains

    subroutine Config_Constants()
      ! subroutine to initialize config_prism.nml namelist

      NAMELIST /Model_config/ &
         INPATH  &
        ,OUTPATH &
        ,GRIDP   &
        ,GRIDZ   

      open( config_IO,file='config_prism.nml',delim='APOSTROPHE')
      read( config_IO,NML=Model_config)
      close(config_IO)

    end subroutine 

end module Config_mod