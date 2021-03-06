MODULE cpl_oasis3
  USE par_oce
  USE dom_oce
  USE in_out_manager
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: cpl_init
  PUBLIC :: cpl_define
  PUBLIC :: cpl_snd
  PUBLIC :: cpl_rcv
  PUBLIC :: cpl_freq
  PUBLIC :: cpl_finalize
  INTEGER, PUBLIC :: OASIS_Rcv = 1
  INTEGER, PUBLIC :: OASIS_idle = 0
  INTEGER :: ncomp_id
  INTEGER :: nerror
  INTEGER :: OASIS_Out = - 1
  INTEGER :: OASIS_REAL = - 1
  INTEGER :: OASIS_Ok = - 1
  INTEGER :: OASIS_In = - 1
  INTEGER :: OASIS_Sent = - 1
  INTEGER :: OASIS_SentOut = - 1
  INTEGER :: OASIS_ToRest = - 1
  INTEGER :: OASIS_ToRestOut = - 1
  INTEGER :: OASIS_Recvd = - 1
  INTEGER :: OASIS_RecvOut = - 1
  INTEGER :: OASIS_FromRest = - 1
  INTEGER :: OASIS_FromRestOut = - 1
  INTEGER :: nrcv
  INTEGER :: nsnd
  INTEGER :: ncplmodel
  INTEGER, PUBLIC, PARAMETER :: nmaxfld = 60
  INTEGER, PUBLIC, PARAMETER :: nmaxcat = 5
  INTEGER, PUBLIC, PARAMETER :: nmaxcpl = 5
  LOGICAL, PARAMETER :: ltmp_wapatch = .TRUE.
  INTEGER :: nldi_save, nlei_save
  INTEGER :: nldj_save, nlej_save
  TYPE, PUBLIC :: FLD_CPL
    LOGICAL :: laction
    CHARACTER(LEN = 8) :: clname
    CHARACTER(LEN = 1) :: clgrid
    REAL(KIND = wp) :: nsgn
    INTEGER, DIMENSION(nmaxcat, nmaxcpl) :: nid
    INTEGER :: nct
    INTEGER :: ncplmodel
  END TYPE FLD_CPL
  TYPE(FLD_CPL), DIMENSION(nmaxfld), PUBLIC :: srcv, ssnd
  REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: exfld
  CONTAINS
  SUBROUTINE cpl_init(cd_modname, kl_comm)
    CHARACTER(LEN = *), INTENT(IN   ) :: cd_modname
    INTEGER, INTENT(  OUT) :: kl_comm
    CALL oasis_init_comp(ncomp_id, TRIM(cd_modname), nerror)
    IF (nerror /= OASIS_Ok) CALL oasis_abort(ncomp_id, 'cpl_init', 'Failure in oasis_init_comp')
    CALL oasis_get_localcomm(kl_comm, nerror)
    IF (nerror /= OASIS_Ok) CALL oasis_abort(ncomp_id, 'cpl_init', 'Failure in oasis_get_localcomm')
  END SUBROUTINE cpl_init
  SUBROUTINE cpl_define(krcv, ksnd, kcplmodel)
    INTEGER, INTENT(IN) :: krcv, ksnd
    INTEGER, INTENT(IN) :: kcplmodel
    INTEGER :: id_part
    INTEGER :: paral(5)
    INTEGER :: ishape(2, 2)
    INTEGER :: ji, jc, jm
    CHARACTER(LEN = 64) :: zclname
    CHARACTER(LEN = 2) :: cli2
    IF (ltmp_wapatch) THEN
      nldi_save = nldi
      nlei_save = nlei
      nldj_save = nldj
      nlej_save = nlej
      IF (nimpp == 1) nldi = 1
      IF (nimpp + jpi - 1 == jpiglo) nlei = jpi
      IF (njmpp == 1) nldj = 1
      IF (njmpp + jpj - 1 == jpjglo) nlej = jpj
    END IF
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) 'cpl_define : initialization in coupled ocean/atmosphere case'
    IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~~~~~~'
    IF (lwp) WRITE(numout, FMT = *)
    ncplmodel = kcplmodel
    IF (kcplmodel > nmaxcpl) THEN
      CALL oasis_abort(ncomp_id, 'cpl_define', 'ncplmodel is larger than nmaxcpl, increase nmaxcpl')
      RETURN
    END IF
    nrcv = krcv
    IF (nrcv > nmaxfld) THEN
      CALL oasis_abort(ncomp_id, 'cpl_define', 'nrcv is larger than nmaxfld, increase nmaxfld')
      RETURN
    END IF
    nsnd = ksnd
    IF (nsnd > nmaxfld) THEN
      CALL oasis_abort(ncomp_id, 'cpl_define', 'nsnd is larger than nmaxfld, increase nmaxfld')
      RETURN
    END IF
    !$ACC KERNELS
    ishape(:, 1) = (/1, nlei - nldi + 1/)
    ishape(:, 2) = (/1, nlej - nldj + 1/)
    !$ACC END KERNELS
    ALLOCATE(exfld(nlei - nldi + 1, nlej - nldj + 1), STAT = nerror)
    IF (nerror > 0) THEN
      CALL oasis_abort(ncomp_id, 'cpl_define', 'Failure in allocating exfld')
      RETURN
    END IF
    paral(1) = 2
    paral(2) = jpiglo * (nldj - 1 + njmpp - 1) + (nldi - 1 + nimpp - 1)
    paral(3) = nlei - nldi + 1
    paral(4) = nlej - nldj + 1
    paral(5) = jpiglo
    IF (ln_ctl) THEN
      WRITE(numout, FMT = *) ' multiexchg: paral (1:5)', paral
      WRITE(numout, FMT = *) ' multiexchg: jpi, jpj =', jpi, jpj
      WRITE(numout, FMT = *) ' multiexchg: nldi, nlei, nimpp =', nldi, nlei, nimpp
      WRITE(numout, FMT = *) ' multiexchg: nldj, nlej, njmpp =', nldj, nlej, njmpp
    END IF
    CALL oasis_def_partition(id_part, paral, nerror, jpiglo * jpjglo)
    ssnd(:) % ncplmodel = kcplmodel
    DO ji = 1, ksnd
      IF (ssnd(ji) % laction) THEN
        IF (ssnd(ji) % nct > nmaxcat) THEN
          CALL oasis_abort(ncomp_id, 'cpl_define', 'Number of categories of ' // TRIM(ssnd(ji) % clname) // ' is larger than nmaxcat, increase nmaxcat')
          RETURN
        END IF
        DO jc = 1, ssnd(ji) % nct
          DO jm = 1, kcplmodel
            IF (ssnd(ji) % nct .GT. 1) THEN
              WRITE(cli2, FMT = '(i2.2)') jc
              zclname = TRIM(ssnd(ji) % clname) // '_cat' // cli2
            ELSE
              zclname = ssnd(ji) % clname
            END IF
            IF (kcplmodel > 1) THEN
              WRITE(cli2, FMT = '(i2.2)') jm
              zclname = 'model' // cli2 // '_' // TRIM(zclname)
            END IF
            IF (ln_ctl) WRITE(numout, FMT = *) "Define", ji, jc, jm, " " // TRIM(zclname), " for ", OASIS_Out
            CALL oasis_def_var(ssnd(ji) % nid(jc, jm), zclname, id_part, (/2, 0/), OASIS_Out, ishape, OASIS_REAL, nerror)
            IF (nerror /= OASIS_Ok) THEN
              WRITE(numout, FMT = *) 'Failed to define transient ', ji, jc, jm, " " // TRIM(zclname)
              CALL oasis_abort(ssnd(ji) % nid(jc, jm), 'cpl_define', 'Failure in oasis_def_var')
            END IF
            IF (ln_ctl .AND. ssnd(ji) % nid(jc, jm) /= - 1) WRITE(numout, FMT = *) "variable defined in the namcouple"
            IF (ln_ctl .AND. ssnd(ji) % nid(jc, jm) == - 1) WRITE(numout, FMT = *) "variable NOT defined in the namcouple"
          END DO
        END DO
      END IF
    END DO
    srcv(:) % ncplmodel = kcplmodel
    DO ji = 1, krcv
      IF (srcv(ji) % laction) THEN
        IF (srcv(ji) % nct > nmaxcat) THEN
          CALL oasis_abort(ncomp_id, 'cpl_define', 'Number of categories of ' // TRIM(srcv(ji) % clname) // ' is larger than nmaxcat, increase nmaxcat')
          RETURN
        END IF
        DO jc = 1, srcv(ji) % nct
          DO jm = 1, kcplmodel
            IF (srcv(ji) % nct .GT. 1) THEN
              WRITE(cli2, FMT = '(i2.2)') jc
              zclname = TRIM(srcv(ji) % clname) // '_cat' // cli2
            ELSE
              zclname = srcv(ji) % clname
            END IF
            IF (kcplmodel > 1) THEN
              WRITE(cli2, FMT = '(i2.2)') jm
              zclname = 'model' // cli2 // '_' // TRIM(zclname)
            END IF
            IF (ln_ctl) WRITE(numout, FMT = *) "Define", ji, jc, jm, " " // TRIM(zclname), " for ", OASIS_In
            CALL oasis_def_var(srcv(ji) % nid(jc, jm), zclname, id_part, (/2, 0/), OASIS_In, ishape, OASIS_REAL, nerror)
            IF (nerror /= OASIS_Ok) THEN
              WRITE(numout, FMT = *) 'Failed to define transient ', ji, jc, jm, " " // TRIM(zclname)
              CALL oasis_abort(srcv(ji) % nid(jc, jm), 'cpl_define', 'Failure in oasis_def_var')
            END IF
            IF (ln_ctl .AND. srcv(ji) % nid(jc, jm) /= - 1) WRITE(numout, FMT = *) "variable defined in the namcouple"
            IF (ln_ctl .AND. srcv(ji) % nid(jc, jm) == - 1) WRITE(numout, FMT = *) "variable NOT defined in the namcouple"
          END DO
        END DO
      END IF
    END DO
    CALL oasis_enddef(nerror)
    IF (nerror /= OASIS_Ok) CALL oasis_abort(ncomp_id, 'cpl_define', 'Failure in oasis_enddef')
    IF (ltmp_wapatch) THEN
      nldi = nldi_save
      nlei = nlei_save
      nldj = nldj_save
      nlej = nlej_save
    END IF
  END SUBROUTINE cpl_define
  SUBROUTINE cpl_snd(kid, kstep, pdata, kinfo)
    INTEGER, INTENT(IN   ) :: kid
    INTEGER, INTENT(  OUT) :: kinfo
    INTEGER, INTENT(IN   ) :: kstep
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN   ) :: pdata
    INTEGER :: jc, jm
    IF (ltmp_wapatch) THEN
      nldi_save = nldi
      nlei_save = nlei
      nldj_save = nldj
      nlej_save = nlej
      IF (nimpp == 1) nldi = 1
      IF (nimpp + jpi - 1 == jpiglo) nlei = jpi
      IF (njmpp == 1) nldj = 1
      IF (njmpp + jpj - 1 == jpjglo) nlej = jpj
    END IF
    DO jc = 1, ssnd(kid) % nct
      DO jm = 1, ssnd(kid) % ncplmodel
        IF (ssnd(kid) % nid(jc, jm) /= - 1) THEN
          CALL oasis_put(ssnd(kid) % nid(jc, jm), kstep, pdata(nldi : nlei, nldj : nlej, jc), kinfo)
          IF (ln_ctl) THEN
            IF (kinfo == OASIS_Sent .OR. kinfo == OASIS_ToRest .OR. kinfo == OASIS_SentOut .OR. kinfo == OASIS_ToRestOut) THEN
              WRITE(numout, FMT = *) '****************'
              WRITE(numout, FMT = *) 'oasis_put: Outgoing ', ssnd(kid) % clname
              WRITE(numout, FMT = *) 'oasis_put: ivarid ', ssnd(kid) % nid(jc, jm)
              WRITE(numout, FMT = *) 'oasis_put:  kstep ', kstep
              WRITE(numout, FMT = *) 'oasis_put:   info ', kinfo
              WRITE(numout, FMT = *) '     - Minimum value is ', MINVAL(pdata(:, :, jc))
              WRITE(numout, FMT = *) '     - Maximum value is ', MAXVAL(pdata(:, :, jc))
              WRITE(numout, FMT = *) '     -     Sum value is ', SUM(pdata(:, :, jc))
              WRITE(numout, FMT = *) '****************'
            END IF
          END IF
        END IF
      END DO
    END DO
    IF (ltmp_wapatch) THEN
      nldi = nldi_save
      nlei = nlei_save
      nldj = nldj_save
      nlej = nlej_save
    END IF
  END SUBROUTINE cpl_snd
  SUBROUTINE cpl_rcv(kid, kstep, pdata, pmask, kinfo)
    INTEGER, INTENT(IN   ) :: kid
    INTEGER, INTENT(IN   ) :: kstep
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pdata
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN   ) :: pmask
    INTEGER, INTENT(  OUT) :: kinfo
    INTEGER :: jc, jm
    LOGICAL :: llaction, llfisrt
    IF (ltmp_wapatch) THEN
      nldi_save = nldi
      nlei_save = nlei
      nldj_save = nldj
      nlej_save = nlej
    END IF
    kinfo = OASIS_idle
    DO jc = 1, srcv(kid) % nct
      IF (ltmp_wapatch) THEN
        IF (nimpp == 1) nldi = 1
        IF (nimpp + jpi - 1 == jpiglo) nlei = jpi
        IF (njmpp == 1) nldj = 1
        IF (njmpp + jpj - 1 == jpjglo) nlej = jpj
      END IF
      llfisrt = .TRUE.
      DO jm = 1, srcv(kid) % ncplmodel
        IF (srcv(kid) % nid(jc, jm) /= - 1) THEN
          CALL oasis_get(srcv(kid) % nid(jc, jm), kstep, exfld, kinfo)
          llaction = kinfo == OASIS_Recvd .OR. kinfo == OASIS_FromRest .OR. kinfo == OASIS_RecvOut .OR. kinfo == OASIS_FromRestOut
          IF (ln_ctl) WRITE(numout, FMT = *) "llaction, kinfo, kstep, ivarid: ", llaction, kinfo, kstep, srcv(kid) % nid(jc, jm)
          IF (llaction) THEN
            kinfo = OASIS_Rcv
            IF (llfisrt) THEN
              !$ACC KERNELS
              pdata(nldi : nlei, nldj : nlej, jc) = exfld(:, :) * pmask(nldi : nlei, nldj : nlej, jm)
              llfisrt = .FALSE.
              !$ACC END KERNELS
            ELSE
              !$ACC KERNELS
              pdata(nldi : nlei, nldj : nlej, jc) = pdata(nldi : nlei, nldj : nlej, jc) + exfld(:, :) * pmask(nldi : nlei, nldj : nlej, jm)
              !$ACC END KERNELS
            END IF
            IF (ln_ctl) THEN
              WRITE(numout, FMT = *) '****************'
              WRITE(numout, FMT = *) 'oasis_get: Incoming ', srcv(kid) % clname
              WRITE(numout, FMT = *) 'oasis_get: ivarid ', srcv(kid) % nid(jc, jm)
              WRITE(numout, FMT = *) 'oasis_get:   kstep', kstep
              WRITE(numout, FMT = *) 'oasis_get:   info ', kinfo
              WRITE(numout, FMT = *) '     - Minimum value is ', MINVAL(pdata(:, :, jc))
              WRITE(numout, FMT = *) '     - Maximum value is ', MAXVAL(pdata(:, :, jc))
              WRITE(numout, FMT = *) '     -     Sum value is ', SUM(pdata(:, :, jc))
              WRITE(numout, FMT = *) '****************'
            END IF
          END IF
        END IF
      END DO
      IF (ltmp_wapatch) THEN
        nldi = nldi_save
        nlei = nlei_save
        nldj = nldj_save
        nlej = nlej_save
      END IF
      IF (.NOT. llfisrt) CALL lbc_lnk(pdata(:, :, jc), srcv(kid) % clgrid, srcv(kid) % nsgn)
    END DO
  END SUBROUTINE cpl_rcv
  INTEGER FUNCTION cpl_freq(cdfieldname)
    CHARACTER(LEN = *), INTENT(IN) :: cdfieldname
    INTEGER :: id
    INTEGER :: info
    INTEGER, DIMENSION(1) :: itmp
    INTEGER :: ji, jm
    INTEGER :: mop
    cpl_freq = 0
    id = - 1
    DO ji = 1, nsnd
      IF (ssnd(ji) % laction) THEN
        DO jm = 1, ncplmodel
          IF (ssnd(ji) % nid(1, jm) /= - 1) THEN
            IF (TRIM(cdfieldname) == TRIM(ssnd(ji) % clname)) THEN
              id = ssnd(ji) % nid(1, jm)
              mop = OASIS_Out
            END IF
          END IF
        END DO
      END IF
    END DO
    DO ji = 1, nrcv
      IF (srcv(ji) % laction) THEN
        DO jm = 1, ncplmodel
          IF (srcv(ji) % nid(1, jm) /= - 1) THEN
            IF (TRIM(cdfieldname) == TRIM(srcv(ji) % clname)) THEN
              id = srcv(ji) % nid(1, jm)
              mop = OASIS_In
            END IF
          END IF
        END DO
      END IF
    END DO
    IF (id /= - 1) THEN
      CALL oasis_get_freqs(id, 1, itmp, info)
      cpl_freq = itmp(1)
    END IF
  END FUNCTION cpl_freq
  SUBROUTINE cpl_finalize
    DEALLOCATE(exfld)
    IF (nstop == 0) THEN
      CALL oasis_terminate(nerror)
    ELSE
      CALL oasis_abort(ncomp_id, "cpl_finalize", "NEMO ABORT STOP")
    END IF
  END SUBROUTINE cpl_finalize
  SUBROUTINE oasis_init_comp(k1, cd1, k2)
    CHARACTER(LEN = *), INTENT(IN   ) :: cd1
    INTEGER, INTENT(  OUT) :: k1, k2
    k1 = - 1
    k2 = - 1
    WRITE(numout, FMT = *) 'oasis_init_comp: Error you sould not be there...', cd1
  END SUBROUTINE oasis_init_comp
  SUBROUTINE oasis_abort(k1, cd1, cd2)
    INTEGER, INTENT(IN   ) :: k1
    CHARACTER(LEN = *), INTENT(IN   ) :: cd1, cd2
    WRITE(numout, FMT = *) 'oasis_abort: Error you sould not be there...', cd1, cd2
  END SUBROUTINE oasis_abort
  SUBROUTINE oasis_get_localcomm(k1, k2)
    INTEGER, INTENT(  OUT) :: k1, k2
    k1 = - 1
    k2 = - 1
    WRITE(numout, FMT = *) 'oasis_get_localcomm: Error you sould not be there...'
  END SUBROUTINE oasis_get_localcomm
  SUBROUTINE oasis_def_partition(k1, k2, k3, k4)
    INTEGER, INTENT(  OUT) :: k1, k3
    INTEGER, INTENT(IN   ) :: k2(5)
    INTEGER, INTENT(IN   ) :: k4
    k1 = k2(1)
    k3 = k2(5) + k4
    WRITE(numout, FMT = *) 'oasis_def_partition: Error you sould not be there...'
  END SUBROUTINE oasis_def_partition
  SUBROUTINE oasis_def_var(k1, cd1, k2, k3, k4, k5, k6, k7)
    CHARACTER(LEN = *), INTENT(IN   ) :: cd1
    INTEGER, INTENT(IN   ) :: k2, k3(2), k4, k5(2, 2), k6
    INTEGER, INTENT(  OUT) :: k1, k7
    k1 = - 1
    k7 = - 1
    WRITE(numout, FMT = *) 'oasis_def_var: Error you sould not be there...', cd1
  END SUBROUTINE oasis_def_var
  SUBROUTINE oasis_enddef(k1)
    INTEGER, INTENT(  OUT) :: k1
    k1 = - 1
    WRITE(numout, FMT = *) 'oasis_enddef: Error you sould not be there...'
  END SUBROUTINE oasis_enddef
  SUBROUTINE oasis_put(k1, k2, p1, k3)
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN   ) :: p1
    INTEGER, INTENT(IN   ) :: k1, k2
    INTEGER, INTENT(  OUT) :: k3
    k3 = - 1
    WRITE(numout, FMT = *) 'oasis_put: Error you sould not be there...'
  END SUBROUTINE oasis_put
  SUBROUTINE oasis_get(k1, k2, p1, k3)
    REAL(KIND = wp), DIMENSION(:, :), INTENT(  OUT) :: p1
    INTEGER, INTENT(IN   ) :: k1, k2
    INTEGER, INTENT(  OUT) :: k3
    p1(1, 1) = - 1.
    k3 = - 1
    WRITE(numout, FMT = *) 'oasis_get: Error you sould not be there...'
  END SUBROUTINE oasis_get
  SUBROUTINE oasis_get_freqs(k1, k2, k3, k4)
    INTEGER, INTENT(IN   ) :: k1, k2
    INTEGER, DIMENSION(1), INTENT(  OUT) :: k3
    INTEGER, INTENT(  OUT) :: k4
    k3(1) = k1
    k4 = k2
    WRITE(numout, FMT = *) 'oasis_get_freqs: Error you sould not be there...'
  END SUBROUTINE oasis_get_freqs
  SUBROUTINE oasis_terminate(k1)
    INTEGER, INTENT(  OUT) :: k1
    k1 = - 1
    WRITE(numout, FMT = *) 'oasis_terminate: Error you sould not be there...'
  END SUBROUTINE oasis_terminate
END MODULE cpl_oasis3