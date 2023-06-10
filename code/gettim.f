         subroutine gettim (ihr, imin, isec, ihs)
         character (len = 10), parameter :: zfmt = '(2i2,f6.3)'
         character (len = 10) :: ztim
         integer, intent (out) :: ihr, imin, isec, ihs     ! elapsed time
         call date_and_time (time = ztim)
         read (ztim, zfmt) ihr, imin, xs
         isec = xs
         ihs = 0
         return
         end subroutine gettim
