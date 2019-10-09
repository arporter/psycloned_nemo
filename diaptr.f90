MODULE diaptr
  USE oce
  USE dom_oce
  USE phycst
  USE iom
  USE in_out_manager
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  INTERFACE ptr_sj
    MODULE PROCEDURE ptr_sj_3d, ptr_sj_2d
  END INTERFACE
  PUBLIC :: ptr_sj
  PUBLIC :: ptr_sjk
  PUBLIC :: dia_ptr_init
  PUBLIC :: dia_ptr
  PUBLIC :: dia_ptr_hst
  REAL(KIND = wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:, :) :: htr_adv, htr_ldf, htr_eiv
  REAL(KIND = wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:, :) :: str_adv, str_ldf, str_eiv
  REAL(KIND = wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:, :) :: htr_ove, str_ove
  REAL(KIND = wp), ALLOCATABLE, SAVE, PUBLIC, DIMENSION(:, :) :: htr_btr, str_btr
  LOGICAL, PUBLIC :: ln_diaptr
  LOGICAL, PUBLIC :: ln_subbas
  INTEGER, PUBLIC :: nptr
  REAL(KIND = wp) :: rc_sv = 1.E-6_wp
  REAL(KIND = wp) :: rc_pwatt = 1.E-15_wp
  REAL(KIND = wp) :: rc_ggram = 1.E-6_wp
  CHARACTER(LEN = 3), ALLOCATABLE, SAVE, DIMENSION(:) :: clsubb
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: btmsk
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: btm30
  REAL(KIND = wp), TARGET, ALLOCATABLE, SAVE, DIMENSION(:) :: p_fval1d
  REAL(KIND = wp), TARGET, ALLOCATABLE, SAVE, DIMENSION(:, :) :: p_fval2d
  CONTAINS
  SUBROUTINE dia_ptr(pvtr)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN), OPTIONAL :: pvtr
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp) :: zsfc, zvfc
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z2d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: z3d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zmask
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, jpts) :: zts
    REAL(KIND = wp), DIMENSION(jpj) :: vsum
    REAL(KIND = wp), DIMENSION(jpj, jpts) :: tssum
    REAL(KIND = wp), DIMENSION(jpj, jpk, nptr) :: sjk, r1_sjk
    REAL(KIND = wp), DIMENSION(jpj, jpk, nptr) :: v_msf, sn_jk, tn_jk
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zvn
    CHARACTER(LEN = 12) :: cl1
    TYPE(ProfileData), SAVE :: psy_profile0
    TYPE(ProfileData), SAVE :: psy_profile1
    TYPE(ProfileData), SAVE :: psy_profile2
    TYPE(ProfileData), SAVE :: psy_profile3
    TYPE(ProfileData), SAVE :: psy_profile4
    TYPE(ProfileData), SAVE :: psy_profile5
    TYPE(ProfileData), SAVE :: psy_profile6
    TYPE(ProfileData), SAVE :: psy_profile7
    TYPE(ProfileData), SAVE :: psy_profile8
    TYPE(ProfileData), SAVE :: psy_profile9
    TYPE(ProfileData), SAVE :: psy_profile10
    TYPE(ProfileData), SAVE :: psy_profile11
    TYPE(ProfileData), SAVE :: psy_profile12
    TYPE(ProfileData), SAVE :: psy_profile13
    TYPE(ProfileData), SAVE :: psy_profile14
    TYPE(ProfileData), SAVE :: psy_profile15
    TYPE(ProfileData), SAVE :: psy_profile16
    TYPE(ProfileData), SAVE :: psy_profile17
    TYPE(ProfileData), SAVE :: psy_profile18
    TYPE(ProfileData), SAVE :: psy_profile19
    TYPE(ProfileData), SAVE :: psy_profile20
    TYPE(ProfileData), SAVE :: psy_profile21
    TYPE(ProfileData), SAVE :: psy_profile22
    TYPE(ProfileData), SAVE :: psy_profile23
    TYPE(ProfileData), SAVE :: psy_profile24
    TYPE(ProfileData), SAVE :: psy_profile25
    TYPE(ProfileData), SAVE :: psy_profile26
    TYPE(ProfileData), SAVE :: psy_profile27
    TYPE(ProfileData), SAVE :: psy_profile28
    TYPE(ProfileData), SAVE :: psy_profile29
    TYPE(ProfileData), SAVE :: psy_profile30
    TYPE(ProfileData), SAVE :: psy_profile31
    TYPE(ProfileData), SAVE :: psy_profile32
    TYPE(ProfileData), SAVE :: psy_profile33
    TYPE(ProfileData), SAVE :: psy_profile34
    IF (ln_timing) CALL timing_start('dia_ptr')
    IF (PRESENT(pvtr)) THEN
      IF (iom_use("zomsfglo")) THEN
        CALL ProfileStart('dia_ptr', 'r0', psy_profile0)
        z3d(1, :, :) = ptr_sjk(pvtr(:, :, :))
        CALL ProfileEnd(psy_profile0)
        !$ACC KERNELS
        DO jk = 2, jpkm1
          z3d(1, :, jk) = z3d(1, :, jk - 1) + z3d(1, :, jk)
        END DO
        DO ji = 1, jpi
          z3d(ji, :, :) = z3d(1, :, :)
        END DO
        !$ACC END KERNELS
        CALL ProfileStart('dia_ptr', 'r1', psy_profile1)
        cl1 = TRIM('zomsf' // clsubb(1))
        CALL iom_put(cl1, z3d * rc_sv)
        CALL ProfileEnd(psy_profile1)
        DO jn = 2, nptr
          CALL ProfileStart('dia_ptr', 'r2', psy_profile2)
          z3d(1, :, :) = ptr_sjk(pvtr(:, :, :), btmsk(:, :, jn) * btm30(:, :))
          CALL ProfileEnd(psy_profile2)
          !$ACC KERNELS
          DO jk = 2, jpkm1
            z3d(1, :, jk) = z3d(1, :, jk - 1) + z3d(1, :, jk)
          END DO
          DO ji = 1, jpi
            z3d(ji, :, :) = z3d(1, :, :)
          END DO
          !$ACC END KERNELS
          CALL ProfileStart('dia_ptr', 'r3', psy_profile3)
          cl1 = TRIM('zomsf' // clsubb(jn))
          CALL iom_put(cl1, z3d * rc_sv)
          CALL ProfileEnd(psy_profile3)
        END DO
      END IF
      IF (iom_use("sopstove") .OR. iom_use("sophtove") .OR. iom_use("sopstbtr") .OR. iom_use("sophtbtr")) THEN
        !$ACC KERNELS
        zmask(:, :, :) = 0._wp
        zts(:, :, :, :) = 0._wp
        zvn(:, :, :) = 0._wp
        DO jk = 1, jpkm1
          DO jj = 1, jpjm1
            DO ji = 1, jpi
              zvfc = e1v(ji, jj) * e3v_n(ji, jj, jk)
              zmask(ji, jj, jk) = vmask(ji, jj, jk) * zvfc
              zts(ji, jj, jk, jp_tem) = (tsn(ji, jj, jk, jp_tem) + tsn(ji, jj + 1, jk, jp_tem)) * 0.5 * zvfc
              zts(ji, jj, jk, jp_sal) = (tsn(ji, jj, jk, jp_sal) + tsn(ji, jj + 1, jk, jp_sal)) * 0.5 * zvfc
              zvn(ji, jj, jk) = vn(ji, jj, jk) * zvfc
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      END IF
      IF (iom_use("sopstove") .OR. iom_use("sophtove")) THEN
        CALL ProfileStart('dia_ptr', 'r4', psy_profile4)
        sjk(:, :, 1) = ptr_sjk(zmask(:, :, :), btmsk(:, :, 1))
        CALL ProfileEnd(psy_profile4)
        !$ACC KERNELS
        r1_sjk(:, :, 1) = 0._wp
        WHERE (sjk(:, :, 1) /= 0._wp) r1_sjk(:, :, 1) = 1._wp / sjk(:, :, 1)
        !$ACC END KERNELS
        CALL ProfileStart('dia_ptr', 'r5', psy_profile5)
        tn_jk(:, :, 1) = ptr_sjk(zts(:, :, :, jp_tem)) * r1_sjk(:, :, 1)
        sn_jk(:, :, 1) = ptr_sjk(zts(:, :, :, jp_sal)) * r1_sjk(:, :, 1)
        v_msf(:, :, 1) = ptr_sjk(zvn(:, :, :))
        CALL ProfileEnd(psy_profile5)
        !$ACC KERNELS
        htr_ove(:, 1) = SUM(v_msf(:, :, 1) * tn_jk(:, :, 1), 2)
        str_ove(:, 1) = SUM(v_msf(:, :, 1) * sn_jk(:, :, 1), 2)
        z2d(1, :) = htr_ove(:, 1) * rc_pwatt
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        !$ACC END KERNELS
        CALL ProfileStart('dia_ptr', 'r6', psy_profile6)
        cl1 = 'sophtove'
        CALL iom_put(TRIM(cl1), z2d)
        CALL ProfileEnd(psy_profile6)
        !$ACC KERNELS
        z2d(1, :) = str_ove(:, 1) * rc_ggram
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        !$ACC END KERNELS
        CALL ProfileStart('dia_ptr', 'r7', psy_profile7)
        cl1 = 'sopstove'
        CALL iom_put(TRIM(cl1), z2d)
        CALL ProfileEnd(psy_profile7)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            CALL ProfileStart('dia_ptr', 'r8', psy_profile8)
            sjk(:, :, jn) = ptr_sjk(zmask(:, :, :), btmsk(:, :, jn))
            CALL ProfileEnd(psy_profile8)
            !$ACC KERNELS
            r1_sjk(:, :, jn) = 0._wp
            WHERE (sjk(:, :, jn) /= 0._wp) r1_sjk(:, :, jn) = 1._wp / sjk(:, :, jn)
            !$ACC END KERNELS
            CALL ProfileStart('dia_ptr', 'r9', psy_profile9)
            tn_jk(:, :, jn) = ptr_sjk(zts(:, :, :, jp_tem), btmsk(:, :, jn)) * r1_sjk(:, :, jn)
            sn_jk(:, :, jn) = ptr_sjk(zts(:, :, :, jp_sal), btmsk(:, :, jn)) * r1_sjk(:, :, jn)
            v_msf(:, :, jn) = ptr_sjk(zvn(:, :, :), btmsk(:, :, jn))
            CALL ProfileEnd(psy_profile9)
            !$ACC KERNELS
            htr_ove(:, jn) = SUM(v_msf(:, :, jn) * tn_jk(:, :, jn), 2)
            str_ove(:, jn) = SUM(v_msf(:, :, jn) * sn_jk(:, :, jn), 2)
            z2d(1, :) = htr_ove(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL ProfileStart('dia_ptr', 'r10', psy_profile10)
            cl1 = TRIM('sophtove_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL ProfileEnd(psy_profile10)
            !$ACC KERNELS
            z2d(1, :) = str_ove(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL ProfileStart('dia_ptr', 'r11', psy_profile11)
            cl1 = TRIM('sopstove_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL ProfileEnd(psy_profile11)
          END DO
        END IF
      END IF
      IF (iom_use("sopstbtr") .OR. iom_use("sophtbtr")) THEN
        CALL ProfileStart('dia_ptr', 'r12', psy_profile12)
        sjk(:, 1, 1) = ptr_sj(zmask(:, :, :), btmsk(:, :, 1))
        CALL ProfileEnd(psy_profile12)
        !$ACC KERNELS
        r1_sjk(:, 1, 1) = 0._wp
        WHERE (sjk(:, 1, 1) /= 0._wp) r1_sjk(:, 1, 1) = 1._wp / sjk(:, 1, 1)
        !$ACC END KERNELS
        CALL ProfileStart('dia_ptr', 'r13', psy_profile13)
        vsum = ptr_sj(zvn(:, :, :), btmsk(:, :, 1))
        tssum(:, jp_tem) = ptr_sj(zts(:, :, :, jp_tem), btmsk(:, :, 1))
        tssum(:, jp_sal) = ptr_sj(zts(:, :, :, jp_sal), btmsk(:, :, 1))
        CALL ProfileEnd(psy_profile13)
        !$ACC KERNELS
        htr_btr(:, 1) = vsum * tssum(:, jp_tem) * r1_sjk(:, 1, 1)
        str_btr(:, 1) = vsum * tssum(:, jp_sal) * r1_sjk(:, 1, 1)
        z2d(1, :) = htr_btr(:, 1) * rc_pwatt
        DO ji = 2, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        !$ACC END KERNELS
        CALL ProfileStart('dia_ptr', 'r14', psy_profile14)
        cl1 = 'sophtbtr'
        CALL iom_put(TRIM(cl1), z2d)
        CALL ProfileEnd(psy_profile14)
        !$ACC KERNELS
        z2d(1, :) = str_btr(:, 1) * rc_ggram
        DO ji = 2, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        !$ACC END KERNELS
        CALL ProfileStart('dia_ptr', 'r15', psy_profile15)
        cl1 = 'sopstbtr'
        CALL iom_put(TRIM(cl1), z2d)
        CALL ProfileEnd(psy_profile15)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            CALL ProfileStart('dia_ptr', 'r16', psy_profile16)
            sjk(:, 1, jn) = ptr_sj(zmask(:, :, :), btmsk(:, :, jn))
            CALL ProfileEnd(psy_profile16)
            !$ACC KERNELS
            r1_sjk(:, 1, jn) = 0._wp
            WHERE (sjk(:, 1, jn) /= 0._wp) r1_sjk(:, 1, jn) = 1._wp / sjk(:, 1, jn)
            !$ACC END KERNELS
            CALL ProfileStart('dia_ptr', 'r17', psy_profile17)
            vsum = ptr_sj(zvn(:, :, :), btmsk(:, :, jn))
            tssum(:, jp_tem) = ptr_sj(zts(:, :, :, jp_tem), btmsk(:, :, jn))
            tssum(:, jp_sal) = ptr_sj(zts(:, :, :, jp_sal), btmsk(:, :, jn))
            CALL ProfileEnd(psy_profile17)
            !$ACC KERNELS
            htr_btr(:, jn) = vsum * tssum(:, jp_tem) * r1_sjk(:, 1, jn)
            str_btr(:, jn) = vsum * tssum(:, jp_sal) * r1_sjk(:, 1, jn)
            z2d(1, :) = htr_btr(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL ProfileStart('dia_ptr', 'r18', psy_profile18)
            cl1 = TRIM('sophtbtr_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL ProfileEnd(psy_profile18)
            !$ACC KERNELS
            z2d(1, :) = str_btr(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL ProfileStart('dia_ptr', 'r19', psy_profile19)
            cl1 = TRIM('sopstbtr_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL ProfileEnd(psy_profile19)
          END DO
        END IF
      END IF
    ELSE
      IF (iom_use("zotemglo")) THEN
        !$ACC KERNELS
        DO jk = 1, jpkm1
          DO jj = 1, jpj
            DO ji = 1, jpi
              zsfc = e1t(ji, jj) * e3t_n(ji, jj, jk)
              zmask(ji, jj, jk) = tmask(ji, jj, jk) * zsfc
              zts(ji, jj, jk, jp_tem) = tsn(ji, jj, jk, jp_tem) * zsfc
              zts(ji, jj, jk, jp_sal) = tsn(ji, jj, jk, jp_sal) * zsfc
            END DO
          END DO
        END DO
        !$ACC END KERNELS
        DO jn = 1, nptr
          CALL ProfileStart('dia_ptr', 'r20', psy_profile20)
          zmask(1, :, :) = ptr_sjk(zmask(:, :, :), btmsk(:, :, jn))
          cl1 = TRIM('zosrf' // clsubb(jn))
          CALL iom_put(cl1, zmask)
          z3d(1, :, :) = ptr_sjk(zts(:, :, :, jp_tem), btmsk(:, :, jn)) / MAX(zmask(1, :, :), 10.E-15)
          CALL ProfileEnd(psy_profile20)
          !$ACC KERNELS
          DO ji = 1, jpi
            z3d(ji, :, :) = z3d(1, :, :)
          END DO
          !$ACC END KERNELS
          CALL ProfileStart('dia_ptr', 'r21', psy_profile21)
          cl1 = TRIM('zotem' // clsubb(jn))
          CALL iom_put(cl1, z3d)
          z3d(1, :, :) = ptr_sjk(zts(:, :, :, jp_sal), btmsk(:, :, jn)) / MAX(zmask(1, :, :), 10.E-15)
          CALL ProfileEnd(psy_profile21)
          !$ACC KERNELS
          DO ji = 1, jpi
            z3d(ji, :, :) = z3d(1, :, :)
          END DO
          !$ACC END KERNELS
          CALL ProfileStart('dia_ptr', 'r22', psy_profile22)
          cl1 = TRIM('zosal' // clsubb(jn))
          CALL iom_put(cl1, z3d)
          CALL ProfileEnd(psy_profile22)
        END DO
      END IF
      IF (iom_use("sophtadv") .OR. iom_use("sopstadv")) THEN
        !$ACC KERNELS
        z2d(1, :) = htr_adv(:, 1) * rc_pwatt
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        !$ACC END KERNELS
        CALL ProfileStart('dia_ptr', 'r23', psy_profile23)
        cl1 = 'sophtadv'
        CALL iom_put(TRIM(cl1), z2d)
        CALL ProfileEnd(psy_profile23)
        !$ACC KERNELS
        z2d(1, :) = str_adv(:, 1) * rc_ggram
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        !$ACC END KERNELS
        CALL ProfileStart('dia_ptr', 'r24', psy_profile24)
        cl1 = 'sopstadv'
        CALL iom_put(TRIM(cl1), z2d)
        CALL ProfileEnd(psy_profile24)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            !$ACC KERNELS
            z2d(1, :) = htr_adv(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL ProfileStart('dia_ptr', 'r25', psy_profile25)
            cl1 = TRIM('sophtadv_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL ProfileEnd(psy_profile25)
            !$ACC KERNELS
            z2d(1, :) = str_adv(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL ProfileStart('dia_ptr', 'r26', psy_profile26)
            cl1 = TRIM('sopstadv_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL ProfileEnd(psy_profile26)
          END DO
        END IF
      END IF
      IF (iom_use("sophtldf") .OR. iom_use("sopstldf")) THEN
        !$ACC KERNELS
        z2d(1, :) = htr_ldf(:, 1) * rc_pwatt
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        !$ACC END KERNELS
        CALL ProfileStart('dia_ptr', 'r27', psy_profile27)
        cl1 = 'sophtldf'
        CALL iom_put(TRIM(cl1), z2d)
        CALL ProfileEnd(psy_profile27)
        !$ACC KERNELS
        z2d(1, :) = str_ldf(:, 1) * rc_ggram
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        !$ACC END KERNELS
        CALL ProfileStart('dia_ptr', 'r28', psy_profile28)
        cl1 = 'sopstldf'
        CALL iom_put(TRIM(cl1), z2d)
        CALL ProfileEnd(psy_profile28)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            !$ACC KERNELS
            z2d(1, :) = htr_ldf(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL ProfileStart('dia_ptr', 'r29', psy_profile29)
            cl1 = TRIM('sophtldf_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL ProfileEnd(psy_profile29)
            !$ACC KERNELS
            z2d(1, :) = str_ldf(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL ProfileStart('dia_ptr', 'r30', psy_profile30)
            cl1 = TRIM('sopstldf_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL ProfileEnd(psy_profile30)
          END DO
        END IF
      END IF
      IF (iom_use("sophteiv") .OR. iom_use("sopsteiv")) THEN
        !$ACC KERNELS
        z2d(1, :) = htr_eiv(:, 1) * rc_pwatt
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        !$ACC END KERNELS
        CALL ProfileStart('dia_ptr', 'r31', psy_profile31)
        cl1 = 'sophteiv'
        CALL iom_put(TRIM(cl1), z2d)
        CALL ProfileEnd(psy_profile31)
        !$ACC KERNELS
        z2d(1, :) = str_eiv(:, 1) * rc_ggram
        DO ji = 1, jpi
          z2d(ji, :) = z2d(1, :)
        END DO
        !$ACC END KERNELS
        CALL ProfileStart('dia_ptr', 'r32', psy_profile32)
        cl1 = 'sopsteiv'
        CALL iom_put(TRIM(cl1), z2d)
        CALL ProfileEnd(psy_profile32)
        IF (ln_subbas) THEN
          DO jn = 2, nptr
            !$ACC KERNELS
            z2d(1, :) = htr_eiv(:, jn) * rc_pwatt
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL ProfileStart('dia_ptr', 'r33', psy_profile33)
            cl1 = TRIM('sophteiv_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL ProfileEnd(psy_profile33)
            !$ACC KERNELS
            z2d(1, :) = str_eiv(:, jn) * rc_ggram
            DO ji = 1, jpi
              z2d(ji, :) = z2d(1, :)
            END DO
            !$ACC END KERNELS
            CALL ProfileStart('dia_ptr', 'r34', psy_profile34)
            cl1 = TRIM('sopsteiv_' // clsubb(jn))
            CALL iom_put(cl1, z2d)
            CALL ProfileEnd(psy_profile34)
          END DO
        END IF
      END IF
    END IF
    IF (ln_timing) CALL timing_stop('dia_ptr')
  END SUBROUTINE dia_ptr
  SUBROUTINE dia_ptr_init
    INTEGER :: jn
    INTEGER :: inum, ierr
    INTEGER :: ios
    NAMELIST /namptr/ ln_diaptr, ln_subbas
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, namptr, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namptr in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, namptr, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namptr in configuration namelist', lwp)
    IF (lwm) WRITE(numond, namptr)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dia_ptr_init : poleward transport and msf initialization'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namptr : set ptr parameters'
      WRITE(numout, FMT = *) '      Poleward heat & salt transport (T) or not (F)      ln_diaptr  = ', ln_diaptr
      WRITE(numout, FMT = *) '      Global (F) or glo/Atl/Pac/Ind/Indo-Pac basins      ln_subbas  = ', ln_subbas
    END IF
    IF (ln_diaptr) THEN
      IF (ln_subbas) THEN
        nptr = 5
        ALLOCATE(clsubb(nptr))
        clsubb(1) = 'glo'
        clsubb(2) = 'atl'
        clsubb(3) = 'pac'
        clsubb(4) = 'ind'
        clsubb(5) = 'ipc'
      ELSE
        nptr = 1
        ALLOCATE(clsubb(nptr))
        clsubb(1) = 'glo'
      END IF
      IF (dia_ptr_alloc() /= 0) CALL ctl_stop('STOP', 'dia_ptr_init : unable to allocate arrays')
      rc_pwatt = rc_pwatt * rau0_rcp
      IF (lk_mpp) CALL mpp_ini_znl(numout)
      IF (ln_subbas) THEN
        CALL iom_open('subbasins', inum, ldstop = .FALSE.)
        CALL iom_get(inum, jpdom_data, 'atlmsk', btmsk(:, :, 2))
        CALL iom_get(inum, jpdom_data, 'pacmsk', btmsk(:, :, 3))
        CALL iom_get(inum, jpdom_data, 'indmsk', btmsk(:, :, 4))
        CALL iom_close(inum)
        !$ACC KERNELS
        btmsk(:, :, 5) = MAX(btmsk(:, :, 3), btmsk(:, :, 4))
        WHERE (gphit(:, :) < - 30._wp)
          btm30(:, :) = 0._wp
        ELSEWHERE
          btm30(:, :) = ssmask(:, :)
        END WHERE
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      btmsk(:, :, 1) = tmask_i(:, :)
      DO jn = 1, nptr
        btmsk(:, :, jn) = btmsk(:, :, jn) * tmask_i(:, :)
      END DO
      htr_adv(:, :) = 0._wp
      str_adv(:, :) = 0._wp
      htr_ldf(:, :) = 0._wp
      str_ldf(:, :) = 0._wp
      htr_eiv(:, :) = 0._wp
      str_eiv(:, :) = 0._wp
      htr_ove(:, :) = 0._wp
      str_ove(:, :) = 0._wp
      htr_btr(:, :) = 0._wp
      str_btr(:, :) = 0._wp
      !$ACC END KERNELS
    END IF
  END SUBROUTINE dia_ptr_init
  SUBROUTINE dia_ptr_hst(ktra, cptr, pva)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT(IN ) :: ktra
    CHARACTER(LEN = 3), INTENT(IN) :: cptr
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pva
    INTEGER :: jn
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('dia_ptr_hst', 'r0', psy_profile0)
    IF (cptr == 'adv') THEN
      IF (ktra == jp_tem) htr_adv(:, 1) = ptr_sj(pva(:, :, :))
      IF (ktra == jp_sal) str_adv(:, 1) = ptr_sj(pva(:, :, :))
    END IF
    IF (cptr == 'ldf') THEN
      IF (ktra == jp_tem) htr_ldf(:, 1) = ptr_sj(pva(:, :, :))
      IF (ktra == jp_sal) str_ldf(:, 1) = ptr_sj(pva(:, :, :))
    END IF
    IF (cptr == 'eiv') THEN
      IF (ktra == jp_tem) htr_eiv(:, 1) = ptr_sj(pva(:, :, :))
      IF (ktra == jp_sal) str_eiv(:, 1) = ptr_sj(pva(:, :, :))
    END IF
    IF (ln_subbas) THEN
      IF (cptr == 'adv') THEN
        IF (ktra == jp_tem) THEN
          DO jn = 2, nptr
            htr_adv(:, jn) = ptr_sj(pva(:, :, :), btmsk(:, :, jn))
          END DO
        END IF
        IF (ktra == jp_sal) THEN
          DO jn = 2, nptr
            str_adv(:, jn) = ptr_sj(pva(:, :, :), btmsk(:, :, jn))
          END DO
        END IF
      END IF
      IF (cptr == 'ldf') THEN
        IF (ktra == jp_tem) THEN
          DO jn = 2, nptr
            htr_ldf(:, jn) = ptr_sj(pva(:, :, :), btmsk(:, :, jn))
          END DO
        END IF
        IF (ktra == jp_sal) THEN
          DO jn = 2, nptr
            str_ldf(:, jn) = ptr_sj(pva(:, :, :), btmsk(:, :, jn))
          END DO
        END IF
      END IF
      IF (cptr == 'eiv') THEN
        IF (ktra == jp_tem) THEN
          DO jn = 2, nptr
            htr_eiv(:, jn) = ptr_sj(pva(:, :, :), btmsk(:, :, jn))
          END DO
        END IF
        IF (ktra == jp_sal) THEN
          DO jn = 2, nptr
            str_eiv(:, jn) = ptr_sj(pva(:, :, :), btmsk(:, :, jn))
          END DO
        END IF
      END IF
   END IF
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE dia_ptr_hst
  FUNCTION dia_ptr_alloc()
    INTEGER :: dia_ptr_alloc
    INTEGER, DIMENSION(3) :: ierr
    ierr(:) = 0
    ALLOCATE(btmsk(jpi, jpj, nptr), htr_adv(jpj, nptr), str_adv(jpj, nptr), htr_eiv(jpj, nptr), str_eiv(jpj, nptr), htr_ove(jpj, nptr), str_ove(jpj, nptr), htr_btr(jpj, nptr), str_btr(jpj, nptr), htr_ldf(jpj, nptr), str_ldf(jpj, nptr), STAT = ierr(1))
    ALLOCATE(p_fval1d(jpj), p_fval2d(jpj, jpk), STAT = ierr(2))
    ALLOCATE(btm30(jpi, jpj), STAT = ierr(3))
    dia_ptr_alloc = MAXVAL(ierr)
    IF (lk_mpp) CALL mpp_sum(dia_ptr_alloc)
  END FUNCTION dia_ptr_alloc
  FUNCTION ptr_sj_3d(pva, pmsk) RESULT(p_fval)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj, jpk) :: pva
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj), OPTIONAL :: pmsk
    INTEGER :: ji, jj, jk
    INTEGER :: ijpj
    REAL(KIND = wp), POINTER, DIMENSION(:) :: p_fval
    p_fval => p_fval1d
    ijpj = jpj
    !$ACC KERNELS
    p_fval(:) = 0._wp
    !$ACC END KERNELS
    IF (PRESENT(pmsk)) THEN
      !$ACC KERNELS
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            p_fval(jj) = p_fval(jj) + pva(ji, jj, jk) * tmask_i(ji, jj) * pmsk(ji, jj)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            p_fval(jj) = p_fval(jj) + pva(ji, jj, jk) * tmask_i(ji, jj)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END IF
  END FUNCTION ptr_sj_3d
  FUNCTION ptr_sj_2d(pva, pmsk) RESULT(p_fval)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj) :: pva
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj), OPTIONAL :: pmsk
    INTEGER :: ji, jj
    INTEGER :: ijpj
    REAL(KIND = wp), POINTER, DIMENSION(:) :: p_fval
    p_fval => p_fval1d
    ijpj = jpj
    !$ACC KERNELS
    p_fval(:) = 0._wp
    !$ACC END KERNELS
    IF (PRESENT(pmsk)) THEN
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = nldi, nlei
          p_fval(jj) = p_fval(jj) + pva(ji, jj) * tmask_i(ji, jj) * pmsk(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = nldi, nlei
          p_fval(jj) = p_fval(jj) + pva(ji, jj) * tmask_i(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    END IF
  END FUNCTION ptr_sj_2d
  FUNCTION ptr_sjk(pta, pmsk) RESULT(p_fval)
    IMPLICIT NONE
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj, jpk) :: pta
    REAL(KIND = wp), INTENT(IN), DIMENSION(jpi, jpj), OPTIONAL :: pmsk
    INTEGER :: ji, jj, jk
    REAL(KIND = wp), POINTER, DIMENSION(:, :) :: p_fval
    p_fval => p_fval2d
    !$ACC KERNELS
    p_fval(:, :) = 0._wp
    !$ACC END KERNELS
    IF (PRESENT(pmsk)) THEN
      !$ACC KERNELS
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = nldi, nlei
            p_fval(jj, jk) = p_fval(jj, jk) + pta(ji, jj, jk) * pmsk(ji, jj)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = nldi, nlei
            p_fval(jj, jk) = p_fval(jj, jk) + pta(ji, jj, jk) * tmask_i(ji, jj)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END IF
  END FUNCTION ptr_sjk
END MODULE diaptr
