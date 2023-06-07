module Config_mod

    implicit none
    private

    public :: Config_Constants ! subroutine to initialize config_prism.nml namelist settings

    integer, public, parameter :: TXTLEN_FILE = 64

    logical :: &   
      FOREST_FIRES     = .true.   !  Forest fire options
  
    integer, public, save ::  &
      NMET = 2     

    character(len=TXTLEN_FILE), public, save :: &
      GRID = 'NOTSET' ! default grid

    contains
      subroutine Config_Constants()

      NAMELIST /Model_config/ &
        GRID &
       ,NMET &
       ,FOREST_FIRES

      open(42,file='config_prism.nml',delim='APOSTROPHE')
      read(42,NML=Model_config)
      close(42)

      end subroutine 

end module Config_mod