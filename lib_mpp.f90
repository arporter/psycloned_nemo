MODULE lib_mpp
  USE in_out_manager
  INTERFACE mpp_sum
    MODULE PROCEDURE mpp_sum_a2s, mpp_sum_as, mpp_sum_ai, mpp_sum_s, mpp_sum_i, mppsum_realdd, mppsum_a_realdd
  END INTERFACE
  INTERFACE mpp_max
    MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
  END INTERFACE
  INTERFACE mpp_min
    MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
  END INTERFACE
  INTERFACE mpp_minloc
    MODULE PROCEDURE mpp_minloc2d, mpp_minloc3d
  END INTERFACE
  INTERFACE mpp_maxloc
    MODULE PROCEDURE mpp_maxloc2d, mpp_maxloc3d
  END INTERFACE
  INTERFACE mpp_max_multiple
    MODULE PROCEDURE mppmax_real_multiple
  END INTERFACE
  LOGICAL, PUBLIC, PARAMETER :: lk_mpp = .FALSE.
  LOGICAL, PUBLIC :: ln_nnogather
  INTEGER :: ncomm_ice
  INTEGER, PUBLIC :: mpi_comm_oce
  CONTAINS
  INTEGER FUNCTION lib_mpp_alloc(kumout)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kumout
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('lib_mpp_alloc', 'r0', 0, 0)
    lib_mpp_alloc = 0
    CALL profile_psy_data0 % PostEnd
  END FUNCTION lib_mpp_alloc
  FUNCTION mynode(ldtxt, ldname, kumnam_ref, knumnam_cfg, kumond, kstop, localComm) RESULT(function_value)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, OPTIONAL, INTENT(IN) :: localComm
    CHARACTER(LEN = *), DIMENSION(:) :: ldtxt
    CHARACTER(LEN = *) :: ldname
    INTEGER :: kumnam_ref, knumnam_cfg, kumond, kstop
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mynode', 'r0', 0, 0)
    IF (PRESENT(localComm)) mpi_comm_oce = localComm
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    function_value = 0
    IF (.FALSE.) ldtxt(:) = 'never done'
    !$ACC END KERNELS
    CALL ctl_opn(kumond, TRIM(ldname), 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', - 1, 6, .FALSE., 1)
  END FUNCTION mynode
  SUBROUTINE mppsync
  END SUBROUTINE mppsync
  SUBROUTINE mpp_sum_as(parr, kdim, kcom)
    REAL, DIMENSION(:) :: parr
    INTEGER :: kdim
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mpp_sum_as: You should not have seen this print! error?', kdim, parr(1), kcom
  END SUBROUTINE mpp_sum_as
  SUBROUTINE mpp_sum_a2s(parr, kdim, kcom)
    REAL, DIMENSION(:, :) :: parr
    INTEGER :: kdim
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mpp_sum_a2s: You should not have seen this print! error?', kdim, parr(1, 1), kcom
  END SUBROUTINE mpp_sum_a2s
  SUBROUTINE mpp_sum_ai(karr, kdim, kcom)
    INTEGER, DIMENSION(:) :: karr
    INTEGER :: kdim
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mpp_sum_ai: You should not have seen this print! error?', kdim, karr(1), kcom
  END SUBROUTINE mpp_sum_ai
  SUBROUTINE mpp_sum_s(psca, kcom)
    REAL :: psca
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mpp_sum_s: You should not have seen this print! error?', psca, kcom
  END SUBROUTINE mpp_sum_s
  SUBROUTINE mpp_sum_i(kint, kcom)
    INTEGER :: kint
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mpp_sum_i: You should not have seen this print! error?', kint, kcom
  END SUBROUTINE mpp_sum_i
  SUBROUTINE mppsum_realdd(ytab, kcom)
    COMPLEX(KIND = wp), INTENT(INOUT) :: ytab
    INTEGER, INTENT(IN), OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mppsum_realdd: You should not have seen this print! error?', ytab
  END SUBROUTINE mppsum_realdd
  SUBROUTINE mppsum_a_realdd(ytab, kdim, kcom)
    INTEGER, INTENT(IN) :: kdim
    COMPLEX(KIND = wp), DIMENSION(kdim), INTENT(INOUT) :: ytab
    INTEGER, INTENT(IN), OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mppsum_a_realdd: You should not have seen this print! error?', kdim, ytab(1), kcom
  END SUBROUTINE mppsum_a_realdd
  SUBROUTINE mppmax_a_real(parr, kdim, kcom)
    REAL, DIMENSION(:) :: parr
    INTEGER :: kdim
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mppmax_a_real: You should not have seen this print! error?', kdim, parr(1), kcom
  END SUBROUTINE mppmax_a_real
  SUBROUTINE mppmax_real(psca, kcom)
    REAL :: psca
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mppmax_real: You should not have seen this print! error?', psca, kcom
  END SUBROUTINE mppmax_real
  SUBROUTINE mppmin_a_real(parr, kdim, kcom)
    REAL, DIMENSION(:) :: parr
    INTEGER :: kdim
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mppmin_a_real: You should not have seen this print! error?', kdim, parr(1), kcom
  END SUBROUTINE mppmin_a_real
  SUBROUTINE mppmin_real(psca, kcom)
    REAL :: psca
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mppmin_real: You should not have seen this print! error?', psca, kcom
  END SUBROUTINE mppmin_real
  SUBROUTINE mppmax_a_int(karr, kdim, kcom)
    INTEGER, DIMENSION(:) :: karr
    INTEGER :: kdim
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mppmax_a_int: You should not have seen this print! error?', kdim, karr(1), kcom
  END SUBROUTINE mppmax_a_int
  SUBROUTINE mppmax_int(kint, kcom)
    INTEGER :: kint
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mppmax_int: You should not have seen this print! error?', kint, kcom
  END SUBROUTINE mppmax_int
  SUBROUTINE mppmin_a_int(karr, kdim, kcom)
    INTEGER, DIMENSION(:) :: karr
    INTEGER :: kdim
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mppmin_a_int: You should not have seen this print! error?', kdim, karr(1), kcom
  END SUBROUTINE mppmin_a_int
  SUBROUTINE mppmin_int(kint, kcom)
    INTEGER :: kint
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mppmin_int: You should not have seen this print! error?', kint, kcom
  END SUBROUTINE mppmin_int
  SUBROUTINE mpp_minloc2d(ptab, pmask, pmin, ki, kj)
    REAL :: pmin
    REAL, DIMENSION(:, :) :: ptab, pmask
    INTEGER :: ki, kj
    WRITE(*, FMT = *) 'mpp_minloc2d: You should not have seen this print! error?', pmin, ki, kj, ptab(1, 1), pmask(1, 1)
  END SUBROUTINE mpp_minloc2d
  SUBROUTINE mpp_minloc3d(ptab, pmask, pmin, ki, kj, kk)
    REAL :: pmin
    REAL, DIMENSION(:, :, :) :: ptab, pmask
    INTEGER :: ki, kj, kk
    WRITE(*, FMT = *) 'mpp_minloc3d: You should not have seen this print! error?', pmin, ki, kj, kk, ptab(1, 1, 1), pmask(1, 1, 1)
  END SUBROUTINE mpp_minloc3d
  SUBROUTINE mpp_maxloc2d(ptab, pmask, pmax, ki, kj)
    REAL :: pmax
    REAL, DIMENSION(:, :) :: ptab, pmask
    INTEGER :: ki, kj
    WRITE(*, FMT = *) 'mpp_maxloc2d: You should not have seen this print! error?', pmax, ki, kj, ptab(1, 1), pmask(1, 1)
  END SUBROUTINE mpp_maxloc2d
  SUBROUTINE mpp_maxloc3d(ptab, pmask, pmax, ki, kj, kk)
    REAL :: pmax
    REAL, DIMENSION(:, :, :) :: ptab, pmask
    INTEGER :: ki, kj, kk
    WRITE(*, FMT = *) 'mpp_maxloc3d: You should not have seen this print! error?', pmax, ki, kj, kk, ptab(1, 1, 1), pmask(1, 1, 1)
  END SUBROUTINE mpp_maxloc3d
  SUBROUTINE mppstop
    STOP
  END SUBROUTINE mppstop
  SUBROUTINE mpp_ini_ice(kcom, knum)
    INTEGER :: kcom, knum
    WRITE(*, FMT = *) 'mpp_ini_ice: You should not have seen this print! error?', kcom, knum
  END SUBROUTINE mpp_ini_ice
  SUBROUTINE mpp_ini_znl(knum)
    INTEGER :: knum
    WRITE(*, FMT = *) 'mpp_ini_znl: You should not have seen this print! error?', knum
  END SUBROUTINE mpp_ini_znl
  SUBROUTINE mpp_comm_free(kcom)
    INTEGER :: kcom
    WRITE(*, FMT = *) 'mpp_comm_free: You should not have seen this print! error?', kcom
  END SUBROUTINE mpp_comm_free
  SUBROUTINE mppmax_real_multiple(ptab, kdim, kcom)
    REAL, DIMENSION(:) :: ptab
    INTEGER :: kdim
    INTEGER, OPTIONAL :: kcom
    WRITE(*, FMT = *) 'mppmax_real_multiple: You should not have seen this print! error?', ptab(1), kdim
  END SUBROUTINE mppmax_real_multiple
  SUBROUTINE ctl_stop(cd1, cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd1, cd2, cd3, cd4, cd5
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd6, cd7, cd8, cd9, cd10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_stop', 'r0', 0, 0)
    nstop = nstop + 1
    IF (lwp) THEN
      WRITE(numout, cform_err)
      IF (PRESENT(cd1)) WRITE(numout, FMT = *) cd1
      IF (PRESENT(cd2)) WRITE(numout, FMT = *) cd2
      IF (PRESENT(cd3)) WRITE(numout, FMT = *) cd3
      IF (PRESENT(cd4)) WRITE(numout, FMT = *) cd4
      IF (PRESENT(cd5)) WRITE(numout, FMT = *) cd5
      IF (PRESENT(cd6)) WRITE(numout, FMT = *) cd6
      IF (PRESENT(cd7)) WRITE(numout, FMT = *) cd7
      IF (PRESENT(cd8)) WRITE(numout, FMT = *) cd8
      IF (PRESENT(cd9)) WRITE(numout, FMT = *) cd9
      IF (PRESENT(cd10)) WRITE(numout, FMT = *) cd10
    END IF
    CALL FLUSH(numout)
    IF (numstp /= - 1) CALL FLUSH(numstp)
    IF (numrun /= - 1) CALL FLUSH(numrun)
    IF (numevo_ice /= - 1) CALL FLUSH(numevo_ice)
    IF (cd1 == 'STOP') THEN
      IF (lwp) WRITE(numout, FMT = *) 'huge E-R-R-O-R : immediate stop'
      CALL mppstop
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ctl_stop
  SUBROUTINE ctl_warn(cd1, cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd1, cd2, cd3, cd4, cd5
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd6, cd7, cd8, cd9, cd10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_warn', 'r0', 0, 0)
    nwarn = nwarn + 1
    IF (lwp) THEN
      WRITE(numout, cform_war)
      IF (PRESENT(cd1)) WRITE(numout, FMT = *) cd1
      IF (PRESENT(cd2)) WRITE(numout, FMT = *) cd2
      IF (PRESENT(cd3)) WRITE(numout, FMT = *) cd3
      IF (PRESENT(cd4)) WRITE(numout, FMT = *) cd4
      IF (PRESENT(cd5)) WRITE(numout, FMT = *) cd5
      IF (PRESENT(cd6)) WRITE(numout, FMT = *) cd6
      IF (PRESENT(cd7)) WRITE(numout, FMT = *) cd7
      IF (PRESENT(cd8)) WRITE(numout, FMT = *) cd8
      IF (PRESENT(cd9)) WRITE(numout, FMT = *) cd9
      IF (PRESENT(cd10)) WRITE(numout, FMT = *) cd10
    END IF
    CALL FLUSH(numout)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ctl_warn
  SUBROUTINE ctl_opn(knum, cdfile, cdstat, cdform, cdacce, klengh, kout, ldwp, karea)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(OUT) :: knum
    CHARACTER(LEN = *), INTENT(IN) :: cdfile
    CHARACTER(LEN = *), INTENT(IN) :: cdstat
    CHARACTER(LEN = *), INTENT(IN) :: cdform
    CHARACTER(LEN = *), INTENT(IN) :: cdacce
    INTEGER, INTENT(IN) :: klengh
    INTEGER, INTENT(IN) :: kout
    LOGICAL, INTENT(IN) :: ldwp
    INTEGER, OPTIONAL, INTENT(IN) :: karea
    CHARACTER(LEN = 80) :: clfile
    INTEGER :: iost
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_opn', 'r0', 0, 0)
    clfile = TRIM(cdfile)
    IF (PRESENT(karea)) THEN
      IF (karea > 1) WRITE(clfile, FMT = "(a,'_',i4.4)") TRIM(clfile), karea - 1
    END IF
    knum = get_unit()
    iost = 0
    IF (cdacce(1 : 6) == 'DIRECT') THEN
      OPEN(UNIT = knum, FILE = clfile, FORM = cdform, ACCESS = cdacce, STATUS = cdstat, RECL = klengh, ERR = 100, IOSTAT = iost)
    ELSE
      OPEN(UNIT = knum, FILE = clfile, FORM = cdform, ACCESS = cdacce, STATUS = cdstat, ERR = 100, IOSTAT = iost)
    END IF
    IF (iost == 0) THEN
      IF (ldwp) THEN
        WRITE(kout, FMT = *) '     file   : ', clfile, ' open ok'
        WRITE(kout, FMT = *) '     unit   = ', knum
        WRITE(kout, FMT = *) '     status = ', cdstat
        WRITE(kout, FMT = *) '     form   = ', cdform
        WRITE(kout, FMT = *) '     access = ', cdacce
        WRITE(kout, FMT = *)
      END IF
    END IF
100 CONTINUE
    IF (iost /= 0) THEN
      IF (ldwp) THEN
        WRITE(kout, FMT = *)
        WRITE(kout, FMT = *) ' ===>>>> : bad opening file: ', clfile
        WRITE(kout, FMT = *) ' =======   ===  '
        WRITE(kout, FMT = *) '           unit   = ', knum
        WRITE(kout, FMT = *) '           status = ', cdstat
        WRITE(kout, FMT = *) '           form   = ', cdform
        WRITE(kout, FMT = *) '           access = ', cdacce
        WRITE(kout, FMT = *) '           iostat = ', iost
        WRITE(kout, FMT = *) '           we stop. verify the file '
        WRITE(kout, FMT = *)
      ELSE
        WRITE(*, FMT = *)
        WRITE(*, FMT = *) ' ===>>>> : bad opening file: ', clfile
        WRITE(*, FMT = *) ' =======   ===  '
        WRITE(*, FMT = *) '           unit   = ', knum
        WRITE(*, FMT = *) '           status = ', cdstat
        WRITE(*, FMT = *) '           form   = ', cdform
        WRITE(*, FMT = *) '           access = ', cdacce
        WRITE(*, FMT = *) '           iostat = ', iost
        WRITE(*, FMT = *) '           we stop. verify the file '
        WRITE(*, FMT = *)
      END IF
      CALL FLUSH(kout)
      STOP 'ctl_opn bad opening'
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ctl_opn
  SUBROUTINE ctl_nam(kios, cdnam, ldwp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(INOUT) :: kios
    CHARACTER(LEN = *), INTENT(IN) :: cdnam
    CHARACTER(LEN = 5) :: clios
    LOGICAL, INTENT(IN) :: ldwp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_nam', 'r0', 0, 0)
    WRITE(clios, FMT = '(I5.0)') kios
    IF (kios < 0) THEN
      CALL ctl_warn('end of record or file while reading namelist ' // TRIM(cdnam) // ' iostat = ' // TRIM(clios))
    END IF
    IF (kios > 0) THEN
      CALL ctl_stop('misspelled variable in namelist ' // TRIM(cdnam) // ' iostat = ' // TRIM(clios))
    END IF
    kios = 0
    RETURN
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ctl_nam
  INTEGER FUNCTION get_unit()
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL :: llopn
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('get_unit', 'r0', 0, 0)
    get_unit = 15
    llopn = .TRUE.
    DO WHILE ((get_unit < 998) .AND. llopn)
      get_unit = get_unit + 1
      INQUIRE(UNIT = get_unit, OPENED = llopn)
    END DO
    IF ((get_unit == 999) .AND. llopn) THEN
      CALL ctl_stop('get_unit: All logical units until 999 are used...')
      get_unit = - 1
    END IF
    CALL profile_psy_data0 % PostEnd
  END FUNCTION get_unit
END MODULE lib_mpp