MODULE timing
  USE in_out_manager
  USE dom_oce
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: timing_init, timing_finalize
  PUBLIC :: timing_reset
  PUBLIC :: timing_start, timing_stop
  TYPE :: timer
    CHARACTER(LEN = 20) :: cname
    REAL(KIND = wp) :: t_cpu, t_clock, tsum_cpu, tsum_clock, tmax_cpu, tmax_clock, tmin_cpu, tmin_clock, tsub_cpu, tsub_clock
    INTEGER :: ncount, ncount_max, ncount_rate
    INTEGER :: niter
    LOGICAL :: l_tdone
    TYPE(timer), POINTER :: next => NULL()
    TYPE(timer), POINTER :: prev => NULL()
    TYPE(timer), POINTER :: parent_section => NULL()
  END TYPE timer
  TYPE :: alltimer
    CHARACTER(LEN = 20), DIMENSION(:), POINTER :: cname => NULL()
    REAL(KIND = wp), DIMENSION(:), POINTER :: tsum_cpu => NULL()
    REAL(KIND = wp), DIMENSION(:), POINTER :: tsum_clock => NULL()
    INTEGER, DIMENSION(:), POINTER :: niter => NULL()
    TYPE(alltimer), POINTER :: next => NULL()
    TYPE(alltimer), POINTER :: prev => NULL()
  END TYPE alltimer
  TYPE(timer), POINTER :: s_timer_root => NULL()
  TYPE(timer), POINTER :: s_timer => NULL()
  TYPE(timer), POINTER :: s_wrk => NULL()
  REAL(KIND = wp) :: t_overclock, t_overcpu
  LOGICAL :: l_initdone = .FALSE.
  INTEGER :: nsize
  REAL(KIND = wp) :: tot_etime, tot_ctime
  REAL(KIND = wp), DIMENSION(2) :: t_elaps, t_cpu
  REAL(KIND = wp), ALLOCATABLE, DIMENSION(:) :: all_etime, all_ctime
  INTEGER :: nfinal_count, ncount, ncount_rate, ncount_max
  INTEGER, DIMENSION(8) :: nvalues
  CHARACTER(LEN = 8), DIMENSION(2) :: cdate
  CHARACTER(LEN = 10), DIMENSION(2) :: ctime
  CHARACTER(LEN = 5) :: czone
  LOGICAL :: ln_onefile = .TRUE.
  LOGICAL :: lwriter
  CONTAINS
  SUBROUTINE timing_start(cdinfo)
    CHARACTER(LEN = *), INTENT(IN) :: cdinfo
    IF (.NOT. l_initdone) THEN
      CALL timing_ini_var(cdinfo)
    ELSE
      s_timer => s_timer_root
      DO WHILE (TRIM(s_timer % cname) /= TRIM(cdinfo))
        IF (ASSOCIATED(s_timer % next)) s_timer => s_timer % next
      END DO
    END IF
    s_timer % l_tdone = .FALSE.
    s_timer % niter = s_timer % niter + 1
    s_timer % t_cpu = 0.
    s_timer % t_clock = 0.
    CALL CPU_TIME(s_timer % t_cpu)
    CALL SYSTEM_CLOCK(COUNT_RATE = s_timer % ncount_rate, COUNT_MAX = s_timer % ncount_max)
    CALL SYSTEM_CLOCK(COUNT = s_timer % ncount)
  END SUBROUTINE timing_start
  SUBROUTINE timing_stop(cdinfo, csection)
    CHARACTER(LEN = *), INTENT(IN) :: cdinfo
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: csection
    INTEGER :: ifinal_count, iperiods
    REAL(KIND = wp) :: zcpu_end, zmpitime
    s_wrk => NULL()
    CALL SYSTEM_CLOCK(COUNT = ifinal_count)
    CALL CPU_TIME(zcpu_end)
    s_timer => s_timer_root
    DO WHILE (TRIM(s_timer % cname) /= TRIM(cdinfo))
      IF (ASSOCIATED(s_timer % next)) s_timer => s_timer % next
    END DO
    s_timer % t_cpu = zcpu_end - s_timer % t_cpu - t_overcpu - s_timer % tsub_cpu
    iperiods = ifinal_count - s_timer % ncount
    IF (ifinal_count < s_timer % ncount) iperiods = iperiods + s_timer % ncount_max
    s_timer % t_clock = REAL(iperiods) / s_timer % ncount_rate - t_overclock - s_timer % tsub_clock
    IF (.NOT. PRESENT(csection)) THEN
      s_wrk => s_timer
      DO WHILE (ASSOCIATED(s_wrk % parent_section))
        s_wrk => s_wrk % parent_section
        s_wrk % tsub_cpu = s_wrk % tsub_cpu + s_timer % t_cpu
        s_wrk % tsub_clock = s_wrk % tsub_clock + s_timer % t_clock
      END DO
    END IF
    s_timer % tsum_clock = s_timer % tsum_clock + s_timer % t_clock
    s_timer % tsum_cpu = s_timer % tsum_cpu + s_timer % t_cpu
    s_timer % tsub_clock = 0.
    s_timer % tsub_cpu = 0.
    s_timer % l_tdone = .TRUE.
  END SUBROUTINE timing_stop
  SUBROUTINE timing_init
    INTEGER :: iperiods, istart_count, ifinal_count
    REAL(KIND = wp) :: zdum
    LOGICAL :: ll_f
    IF (ln_onefile) THEN
      IF (lwp) CALL ctl_opn(numtime, 'timing.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, .TRUE., narea)
      lwriter = lwp
    ELSE
      CALL ctl_opn(numtime, 'timing.output', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, .FALSE., narea)
      lwriter = .TRUE.
    END IF
    IF (lwriter) THEN
      WRITE(numtime, FMT = *)
      WRITE(numtime, FMT = *) '      CNRS - NERC - Met OFFICE - MERCATOR-ocean - CMCC - INGV'
      WRITE(numtime, FMT = *) '                             NEMO team'
      WRITE(numtime, FMT = *) '                  Ocean General Circulation Model'
      WRITE(numtime, FMT = *) '                        version 4.0  (2018) '
      WRITE(numtime, FMT = *)
      WRITE(numtime, FMT = *) '                        Timing Informations '
      WRITE(numtime, FMT = *)
      WRITE(numtime, FMT = *)
    END IF
    CALL SYSTEM_CLOCK(COUNT_RATE = ncount_rate, COUNT_MAX = ncount_max)
    CALL SYSTEM_CLOCK(COUNT = istart_count)
    CALL SYSTEM_CLOCK(COUNT = ifinal_count)
    iperiods = ifinal_count - istart_count
    IF (ifinal_count < istart_count) iperiods = iperiods + ncount_max
    t_overclock = REAL(iperiods) / ncount_rate
    CALL CPU_TIME(zdum)
    CALL CPU_TIME(t_overcpu)
    t_overcpu = t_overcpu - zdum
    t_overclock = t_overcpu + t_overclock
    CALL DATE_AND_TIME(cdate(1), ctime(1), czone, nvalues)
    CALL CPU_TIME(t_cpu(1))
    CALL SYSTEM_CLOCK(COUNT_RATE = ncount_rate, COUNT_MAX = ncount_max)
    CALL SYSTEM_CLOCK(COUNT = ncount)
  END SUBROUTINE timing_init
  SUBROUTINE timing_finalize
    TYPE(timer), POINTER :: s_temp
    INTEGER :: idum, iperiods, icode
    LOGICAL :: ll_ord, ll_averep
    CHARACTER(LEN = 120) :: clfmt
    ll_averep = .TRUE.
    CALL CPU_TIME(t_cpu(2))
    t_cpu(2) = t_cpu(2) - t_cpu(1) - t_overcpu
    CALL SYSTEM_CLOCK(COUNT = nfinal_count)
    iperiods = nfinal_count - ncount
    IF (nfinal_count < ncount) iperiods = iperiods + ncount_max
    t_elaps(2) = REAL(iperiods) / ncount_rate - t_overclock
    CALL DATE_AND_TIME(cdate(2), ctime(2), czone, nvalues)
    nsize = 0
    s_timer => s_timer_root
    DO WHILE (ASSOCIATED(s_timer))
      nsize = nsize + 1
      s_timer => s_timer % next
    END DO
    idum = nsize
    IF (lk_mpp) CALL mpp_sum(idum)
    IF (idum / jpnij /= nsize) THEN
      IF (lwriter) WRITE(numtime, FMT = *) '        ===> W A R N I N G: '
      IF (lwriter) WRITE(numtime, FMT = *) ' Some CPU have different number of routines instrumented for timing'
      IF (lwriter) WRITE(numtime, FMT = *) ' No detailed report on averaged timing can be provided'
      IF (lwriter) WRITE(numtime, FMT = *) ' The following detailed report only deals with the current processor'
      IF (lwriter) WRITE(numtime, FMT = *)
      ll_averep = .FALSE.
    END IF
    tot_etime = t_elaps(2)
    tot_ctime = t_cpu(2)
    IF (lwriter) WRITE(numtime, FMT = *) 'Total timing (sum) :'
    IF (lwriter) WRITE(numtime, FMT = *) '--------------------'
    IF (lwriter) WRITE(numtime, FMT = "('Elapsed Time (s)  CPU Time (s)')")
    IF (lwriter) WRITE(numtime, FMT = '(5x,f12.3,1x,f12.3)') tot_etime, tot_ctime
    IF (lwriter) WRITE(numtime, FMT = *)
    IF (lwriter) CALL wcurrent_info
    clfmt = '(1X,"Timing started on ",2(A2,"/"),A4," at ",2(A2,":"),A2," MET ",A3,":",A2," from GMT")'
    IF (lwriter) WRITE(numtime, TRIM(clfmt)) cdate(1)(7 : 8), cdate(1)(5 : 6), cdate(1)(1 : 4), ctime(1)(1 : 2), ctime(1)(3 : 4), ctime(1)(5 : 6), czone(1 : 3), czone(4 : 5)
    clfmt = '(1X,  "Timing   ended on ",2(A2,"/"),A4," at ",2(A2,":"),A2," MET ",A3,":",A2," from GMT")'
    IF (lwriter) WRITE(numtime, TRIM(clfmt)) cdate(2)(7 : 8), cdate(2)(5 : 6), cdate(2)(1 : 4), ctime(2)(1 : 2), ctime(2)(3 : 4), ctime(2)(5 : 6), czone(1 : 3), czone(4 : 5)
    IF (lwriter) CLOSE(UNIT = numtime)
  END SUBROUTINE timing_finalize
  SUBROUTINE wcurrent_info
    LOGICAL :: ll_ord
    CHARACTER(LEN = 2048) :: clfmt
    s_wrk => NULL()
    s_timer => s_timer_root
    DO
      ll_ord = .TRUE.
      s_timer => s_timer_root
      DO WHILE (ASSOCIATED(s_timer % next))
        IF (.NOT. ASSOCIATED(s_timer % next)) EXIT
        IF (s_timer % tsum_clock < s_timer % next % tsum_clock) THEN
          ALLOCATE(s_wrk)
          s_wrk = s_timer % next
          CALL insert(s_timer, s_timer_root, s_wrk)
          CALL suppress(s_timer % next)
          ll_ord = .FALSE.
          CYCLE
        END IF
        IF (ASSOCIATED(s_timer % next)) s_timer => s_timer % next
      END DO
      IF (ll_ord) EXIT
    END DO
    WRITE(numtime, FMT = *) 'Detailed timing for proc :', narea - 1
    WRITE(numtime, FMT = *) '--------------------------'
    WRITE(numtime, FMT = *) 'Section             ', 'Elapsed Time (s)  ', 'Elapsed Time (%)  ', 'CPU Time(s)  ', 'CPU Time (%)  ', 'CPU/Elapsed  ', 'Frequency'
    s_timer => s_timer_root
    clfmt = '(1x,a,4x,f12.3,6x,f12.3,x,f12.3,2x,f12.3,6x,f7.3,2x,i9)'
    DO WHILE (ASSOCIATED(s_timer))
      WRITE(numtime, TRIM(clfmt)) s_timer % cname, s_timer % tsum_clock, s_timer % tsum_clock * 100. / t_elaps(2), s_timer % tsum_cpu, s_timer % tsum_cpu * 100. / t_cpu(2), s_timer % tsum_cpu / s_timer % tsum_clock, s_timer % niter
      s_timer => s_timer % next
    END DO
    WRITE(numtime, FMT = *)
  END SUBROUTINE wcurrent_info
  SUBROUTINE timing_ini_var(cdinfo)
    CHARACTER(LEN = *), INTENT(IN) :: cdinfo
    LOGICAL :: ll_section
    IF (.NOT. ASSOCIATED(s_timer_root)) THEN
      ALLOCATE(s_timer_root)
      s_timer_root % cname = cdinfo
      s_timer_root % t_cpu = 0._wp
      s_timer_root % t_clock = 0._wp
      s_timer_root % tsum_cpu = 0._wp
      s_timer_root % tsum_clock = 0._wp
      s_timer_root % tmax_cpu = 0._wp
      s_timer_root % tmax_clock = 0._wp
      s_timer_root % tmin_cpu = 0._wp
      s_timer_root % tmin_clock = 0._wp
      s_timer_root % tsub_cpu = 0._wp
      s_timer_root % tsub_clock = 0._wp
      s_timer_root % ncount = 0
      s_timer_root % ncount_rate = 0
      s_timer_root % ncount_max = 0
      s_timer_root % niter = 0
      s_timer_root % l_tdone = .FALSE.
      s_timer_root % next => NULL()
      s_timer_root % prev => NULL()
      s_timer => s_timer_root
      s_wrk => NULL()
    ELSE
      s_timer => s_timer_root
      DO WHILE (ASSOCIATED(s_timer))
        IF (TRIM(s_timer % cname) .EQ. TRIM(cdinfo)) RETURN
        s_timer => s_timer % next
      END DO
      s_timer => s_timer_root
      DO WHILE (ASSOCIATED(s_timer % next))
        s_timer => s_timer % next
      END DO
      ALLOCATE(s_timer % next)
      s_timer % next % cname = cdinfo
      s_timer % next % t_cpu = 0._wp
      s_timer % next % t_clock = 0._wp
      s_timer % next % tsum_cpu = 0._wp
      s_timer % next % tsum_clock = 0._wp
      s_timer % next % tmax_cpu = 0._wp
      s_timer % next % tmax_clock = 0._wp
      s_timer % next % tmin_cpu = 0._wp
      s_timer % next % tmin_clock = 0._wp
      s_timer % next % tsub_cpu = 0._wp
      s_timer % next % tsub_clock = 0._wp
      s_timer % next % ncount = 0
      s_timer % next % ncount_rate = 0
      s_timer % next % ncount_max = 0
      s_timer % next % niter = 0
      s_timer % next % l_tdone = .FALSE.
      s_timer % next % parent_section => NULL()
      s_timer % next % prev => s_timer
      s_timer % next % next => NULL()
      s_timer => s_timer % next
      s_wrk => s_timer % prev
      ll_section = .FALSE.
      DO WHILE (ASSOCIATED(s_wrk) .AND. .NOT. ll_section)
        IF (.NOT. s_wrk % l_tdone) THEN
          ll_section = .TRUE.
          s_timer % parent_section => s_wrk
        END IF
        s_wrk => s_wrk % prev
      END DO
    END IF
  END SUBROUTINE timing_ini_var
  SUBROUTINE timing_reset
    l_initdone = .TRUE.
  END SUBROUTINE timing_reset
  RECURSIVE SUBROUTINE timing_list(ptr)
    TYPE(timer), POINTER, INTENT(INOUT) :: ptr
    IF (ASSOCIATED(ptr % next)) CALL timing_list(ptr % next)
    IF (lwp) WRITE(numout, FMT = *) '   ', ptr % cname
  END SUBROUTINE timing_list
  SUBROUTINE insert(sd_current, sd_root, sd_ptr)
    TYPE(timer), POINTER, INTENT(INOUT) :: sd_current, sd_root, sd_ptr
    IF (ASSOCIATED(sd_current, sd_root)) THEN
      sd_root => sd_ptr
    ELSE
      sd_current % prev % next => sd_ptr
    END IF
    sd_ptr % next => sd_current
    sd_ptr % prev => sd_current % prev
    sd_current % prev => sd_ptr
    sd_ptr => NULL()
  END SUBROUTINE insert
  SUBROUTINE suppress(sd_ptr)
    TYPE(timer), POINTER, INTENT(INOUT) :: sd_ptr
    TYPE(timer), POINTER :: sl_temp
    sl_temp => sd_ptr
    sd_ptr => sd_ptr % next
    IF (ASSOCIATED(sl_temp % next)) sl_temp % next % prev => sl_temp % prev
    DEALLOCATE(sl_temp)
    sl_temp => NULL()
  END SUBROUTINE suppress
END MODULE timing