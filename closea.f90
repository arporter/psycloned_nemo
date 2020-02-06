MODULE closea
  USE oce
  USE dom_oce
  USE phycst
  USE sbc_oce
  USE iom
  USE in_out_manager
  USE lib_fortran, ONLY: glob_sum
  USE lbclnk
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dom_clo
  PUBLIC :: sbc_clo
  PUBLIC :: clo_rnf
  PUBLIC :: clo_bat
  LOGICAL, PUBLIC :: ln_closea
  LOGICAL, PUBLIC :: l_sbc_clo
  LOGICAL, PUBLIC :: l_clo_rnf
  INTEGER, PUBLIC :: jncs
  INTEGER, PUBLIC :: jncsr
  INTEGER, PUBLIC :: jncse
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: closea_mask
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: closea_mask_rnf
  INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:, :) :: closea_mask_empmr
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: surf
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: surfr
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, DIMENSION(:) :: surfe
  CONTAINS
  SUBROUTINE dom_clo
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER :: inum
    INTEGER :: ierr
    INTEGER :: id
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdata_in
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
    CALL ProfileStart('dom_clo', 'r0', psy_profile0)
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) 'dom_clo : read in masks to define closed seas '
    IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
    CALL ProfileEnd(psy_profile0)
    IF (ln_read_cfg) THEN
      CALL ProfileStart('dom_clo', 'r1', psy_profile1)
      CALL iom_open(cn_domcfg, inum)
      id = iom_varid(inum, 'closea_mask', ldstop = .FALSE.)
      CALL ProfileEnd(psy_profile1)
      IF (id > 0) THEN
        CALL ProfileStart('dom_clo', 'r2', psy_profile2)
        l_sbc_clo = .TRUE.
        ALLOCATE(closea_mask(jpi, jpj), STAT = ierr)
        IF (ierr /= 0) CALL ctl_stop('STOP', 'dom_clo: failed to allocate closea_mask array')
        CALL ProfileEnd(psy_profile2)
        !$ACC KERNELS
        zdata_in(:, :) = 0.0
        !$ACC END KERNELS
        CALL iom_get(inum, jpdom_data, 'closea_mask', zdata_in)
        !$ACC KERNELS
        closea_mask(:, :) = NINT(zdata_in(:, :)) * tmask(:, :, 1)
        !$ACC END KERNELS
        CALL ProfileStart('dom_clo', 'r3', psy_profile3)
        jncs = MAXVAL(closea_mask(:, :))
        IF (lk_mpp) CALL mpp_max(jncs)
        IF (jncs > 0) THEN
          IF (lwp) WRITE(numout, FMT = *) 'Number of closed seas : ', jncs
        ELSE
          CALL ctl_stop('Problem with closea_mask field in domain_cfg file. Has no values > 0 so no closed seas defined.')
        END IF
        CALL ProfileEnd(psy_profile3)
      ELSE
        CALL ProfileStart('dom_clo', 'r4', psy_profile4)
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   closea_mask field not found in domain_cfg file.'
        IF (lwp) WRITE(numout, FMT = *) '           No closed seas defined.'
        IF (lwp) WRITE(numout, FMT = *)
        l_sbc_clo = .FALSE.
        jncs = 0
        CALL ProfileEnd(psy_profile4)
      END IF
      CALL ProfileStart('dom_clo', 'r5', psy_profile5)
      l_clo_rnf = .FALSE.
      CALL ProfileEnd(psy_profile5)
      IF (l_sbc_clo) THEN
        CALL ProfileStart('dom_clo', 'r6', psy_profile6)
        id = iom_varid(inum, 'closea_mask_rnf', ldstop = .FALSE.)
        CALL ProfileEnd(psy_profile6)
        IF (id > 0) THEN
          CALL ProfileStart('dom_clo', 'r7', psy_profile7)
          l_clo_rnf = .TRUE.
          ALLOCATE(closea_mask_rnf(jpi, jpj), STAT = ierr)
          IF (ierr /= 0) CALL ctl_stop('STOP', 'dom_clo: failed to allocate closea_mask_rnf array')
          CALL iom_get(inum, jpdom_data, 'closea_mask_rnf', zdata_in)
          CALL ProfileEnd(psy_profile7)
          !$ACC KERNELS
          closea_mask_rnf(:, :) = NINT(zdata_in(:, :)) * tmask(:, :, 1)
          !$ACC END KERNELS
          CALL ProfileStart('dom_clo', 'r8', psy_profile8)
          jncsr = MAXVAL(closea_mask_rnf(:, :))
          IF (lk_mpp) CALL mpp_max(jncsr)
          IF (jncsr > 0) THEN
            IF (lwp) WRITE(numout, FMT = *) 'Number of closed seas rnf mappings : ', jncsr
          ELSE
            CALL ctl_stop('Problem with closea_mask_rnf field in domain_cfg file. Has no values > 0 so no closed seas rnf mappings defined.')
          END IF
          CALL ProfileEnd(psy_profile8)
        ELSE
          CALL ProfileStart('dom_clo', 'r9', psy_profile9)
          IF (lwp) WRITE(numout, FMT = *) 'closea_mask_rnf field not found in domain_cfg file. No closed seas rnf mappings defined.'
          jncsr = 0
          CALL ProfileEnd(psy_profile9)
        END IF
        CALL ProfileStart('dom_clo', 'r10', psy_profile10)
        id = iom_varid(inum, 'closea_mask_empmr', ldstop = .FALSE.)
        CALL ProfileEnd(psy_profile10)
        IF (id > 0) THEN
          CALL ProfileStart('dom_clo', 'r11', psy_profile11)
          l_clo_rnf = .TRUE.
          ALLOCATE(closea_mask_empmr(jpi, jpj), STAT = ierr)
          IF (ierr /= 0) CALL ctl_stop('STOP', 'dom_clo: failed to allocate closea_mask_empmr array')
          CALL iom_get(inum, jpdom_data, 'closea_mask_empmr', zdata_in)
          CALL ProfileEnd(psy_profile11)
          !$ACC KERNELS
          closea_mask_empmr(:, :) = NINT(zdata_in(:, :)) * tmask(:, :, 1)
          !$ACC END KERNELS
          CALL ProfileStart('dom_clo', 'r12', psy_profile12)
          jncse = MAXVAL(closea_mask_empmr(:, :))
          IF (lk_mpp) CALL mpp_max(jncse)
          IF (jncse > 0) THEN
            IF (lwp) WRITE(numout, FMT = *) 'Number of closed seas empmr mappings : ', jncse
          ELSE
            CALL ctl_stop('Problem with closea_mask_empmr field in domain_cfg file. Has no values > 0 so no closed seas empmr mappings defined.')
          END IF
          CALL ProfileEnd(psy_profile12)
        ELSE
          CALL ProfileStart('dom_clo', 'r13', psy_profile13)
          IF (lwp) WRITE(numout, FMT = *) 'closea_mask_empmr field not found in domain_cfg file. No closed seas empmr mappings defined.'
          jncse = 0
          CALL ProfileEnd(psy_profile13)
        END IF
      END IF
      CALL iom_close(inum)
    ELSE
      CALL ProfileStart('dom_clo', 'r14', psy_profile14)
      IF (lwp) WRITE(numout, FMT = *) 'No domain_cfg file so no closed seas defined.'
      l_sbc_clo = .FALSE.
      l_clo_rnf = .FALSE.
      CALL ProfileEnd(psy_profile14)
    END IF
  END SUBROUTINE dom_clo
  SUBROUTINE sbc_clo(kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT(IN ) :: kt
    INTEGER :: ierr
    INTEGER :: jc, jcr, jce
    REAL(KIND = wp), PARAMETER :: rsmall = 1.E-20_wp
    REAL(KIND = wp) :: zfwf_total, zcoef, zcoef1
    REAL(KIND = wp), DIMENSION(jncs) :: zfwf
    REAL(KIND = wp), DIMENSION(jncsr + 1) :: zfwfr
    REAL(KIND = wp), DIMENSION(jncse + 1) :: zfwfe
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztmp2d
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
    IF (ln_timing) CALL timing_start('sbc_clo')
    IF (kt == nit000) THEN
      CALL ProfileStart('sbc_clo', 'r0', psy_profile0)
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'sbc_clo : closed seas '
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
      ALLOCATE(surf(jncs + 1), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('STOP', 'sbc_clo: failed to allocate surf array')
      surf(:) = 0.E0_wp
      ALLOCATE(surfr(jncsr + 1), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('STOP', 'sbc_clo: failed to allocate surfr array')
      surfr(:) = 0.E0_wp
      ALLOCATE(surfe(jncse + 1), STAT = ierr)
      IF (ierr /= 0) CALL ctl_stop('STOP', 'sbc_clo: failed to allocate surfe array')
      surfe(:) = 0.E0_wp
      surf(jncs + 1) = glob_sum(e1e2t(:, :))
      CALL ProfileEnd(psy_profile0)
      DO jc = 1, jncs
        !$ACC KERNELS
        ztmp2d(:, :) = 0.E0_wp
        WHERE (closea_mask(:, :) == jc) ztmp2d(:, :) = e1e2t(:, :) * tmask_i(:, :)
        !$ACC END KERNELS
        CALL ProfileStart('sbc_clo', 'r1', psy_profile1)
        surf(jc) = glob_sum(ztmp2d(:, :))
        CALL ProfileEnd(psy_profile1)
      END DO
      CALL ProfileStart('sbc_clo', 'r2', psy_profile2)
      surf(jncs + 1) = surf(jncs + 1) - SUM(surf(1 : jncs))
      CALL ProfileEnd(psy_profile2)
      IF (jncsr > 0) THEN
        DO jcr = 1, jncsr
          !$ACC KERNELS
          ztmp2d(:, :) = 0.E0_wp
          WHERE (closea_mask_rnf(:, :) == jcr .AND. closea_mask(:, :) == 0) ztmp2d(:, :) = e1e2t(:, :) * tmask_i(:, :)
          !$ACC END KERNELS
          CALL ProfileStart('sbc_clo', 'r3', psy_profile3)
          surfr(jcr) = glob_sum(ztmp2d(:, :))
          CALL ProfileEnd(psy_profile3)
        END DO
      END IF
      IF (jncse > 0) THEN
        DO jce = 1, jncse
          !$ACC KERNELS
          ztmp2d(:, :) = 0.E0_wp
          WHERE (closea_mask_empmr(:, :) == jce .AND. closea_mask(:, :) == 0) ztmp2d(:, :) = e1e2t(:, :) * tmask_i(:, :)
          !$ACC END KERNELS
          CALL ProfileStart('sbc_clo', 'r4', psy_profile4)
          surfe(jce) = glob_sum(ztmp2d(:, :))
          CALL ProfileEnd(psy_profile4)
        END DO
      END IF
      CALL ProfileStart('sbc_clo', 'r5', psy_profile5)
      IF (lwp) WRITE(numout, FMT = *) '     Closed sea surface areas (km2)'
      DO jc = 1, jncs
        IF (lwp) WRITE(numout, FMT = '(1I3,5X,ES12.2)') jc, surf(jc) * 1.0E-6
      END DO
      IF (lwp) WRITE(numout, FMT = '(A,ES12.2)') 'Global surface area excluding closed seas (km2): ', surf(jncs + 1) * 1.0E-6
      IF (jncsr > 0) THEN
        IF (lwp) WRITE(numout, FMT = *) '     Closed sea target rnf surface areas (km2)'
        DO jcr = 1, jncsr
          IF (lwp) WRITE(numout, FMT = '(1I3,5X,ES12.2)') jcr, surfr(jcr) * 1.0E-6
        END DO
      END IF
      IF (jncse > 0) THEN
        IF (lwp) WRITE(numout, FMT = *) '     Closed sea target empmr surface areas (km2)'
        DO jce = 1, jncse
          IF (lwp) WRITE(numout, FMT = '(1I3,5X,ES12.2)') jce, surfe(jce) * 1.0E-6
        END DO
      END IF
      CALL ProfileEnd(psy_profile5)
    END IF
    CALL ProfileStart('sbc_clo', 'r6', psy_profile6)
    zfwf_total = 0._wp
    zfwf(:) = 0.E0_wp
    CALL ProfileEnd(psy_profile6)
    DO jc = 1, jncs
      !$ACC KERNELS
      ztmp2d(:, :) = 0.E0_wp
      WHERE (closea_mask(:, :) == jc) ztmp2d(:, :) = e1e2t(:, :) * (emp(:, :) - rnf(:, :)) * tmask_i(:, :)
      !$ACC END KERNELS
      CALL ProfileStart('sbc_clo', 'r7', psy_profile7)
      zfwf(jc) = glob_sum(ztmp2d(:, :))
      CALL ProfileEnd(psy_profile7)
    END DO
    CALL ProfileStart('sbc_clo', 'r8', psy_profile8)
    zfwf_total = SUM(zfwf)
    zfwfr(:) = 0.E0_wp
    CALL ProfileEnd(psy_profile8)
    IF (jncsr > 0) THEN
      DO jcr = 1, jncsr
        !$ACC KERNELS
        ztmp2d(:, :) = 0.E0_wp
        WHERE (closea_mask_rnf(:, :) == jcr .AND. closea_mask(:, :) > 0) ztmp2d(:, :) = e1e2t(:, :) * (emp(:, :) - rnf(:, :)) * tmask_i(:, :)
        !$ACC END KERNELS
        CALL ProfileStart('sbc_clo', 'r9', psy_profile9)
        zfwfr(jcr) = glob_sum(ztmp2d(:, :))
        CALL ProfileEnd(psy_profile9)
        IF (ABS(zfwfr(jcr) / surf(jncs + 1)) > rsmall) THEN
          IF (zfwfr(jcr) < 0.0) THEN
            !$ACC KERNELS
            zfwf_total = zfwf_total - zfwfr(jcr)
            zcoef = zfwfr(jcr) / surfr(jcr)
            zcoef1 = rcp * zcoef
            WHERE (closea_mask_rnf(:, :) == jcr .AND. closea_mask(:, :) == 0.0)
              emp(:, :) = emp(:, :) + zcoef
              qns(:, :) = qns(:, :) - zcoef1 * sst_m(:, :)
            END WHERE
            !$ACC END KERNELS
          END IF
        END IF
      END DO
    END IF
    CALL ProfileStart('sbc_clo', 'r10', psy_profile10)
    zfwfe(:) = 0.E0_wp
    CALL ProfileEnd(psy_profile10)
    IF (jncse > 0) THEN
      DO jce = 1, jncse
        !$ACC KERNELS
        ztmp2d(:, :) = 0.E0_wp
        WHERE (closea_mask_empmr(:, :) == jce .AND. closea_mask(:, :) > 0) ztmp2d(:, :) = e1e2t(:, :) * (emp(:, :) - rnf(:, :)) * tmask_i(:, :)
        !$ACC END KERNELS
        CALL ProfileStart('sbc_clo', 'r11', psy_profile11)
        zfwfe(jce) = glob_sum(ztmp2d(:, :))
        CALL ProfileEnd(psy_profile11)
        IF (ABS(zfwfe(jce) / surf(jncs + 1)) > rsmall) THEN
          !$ACC KERNELS
          zfwf_total = zfwf_total - zfwfe(jce)
          zcoef = zfwfe(jce) / surfe(jce)
          zcoef1 = rcp * zcoef
          WHERE (closea_mask_empmr(:, :) == jce .AND. closea_mask(:, :) == 0.0)
            emp(:, :) = emp(:, :) + zcoef
            qns(:, :) = qns(:, :) - zcoef1 * sst_m(:, :)
          END WHERE
          !$ACC END KERNELS
        END IF
      END DO
    END IF
    IF (ABS(zfwf_total / surf(jncs + 1)) > rsmall) THEN
      !$ACC KERNELS
      zcoef = zfwf_total / surf(jncs + 1)
      zcoef1 = rcp * zcoef
      WHERE (closea_mask(:, :) == 0)
        emp(:, :) = emp(:, :) + zcoef
        qns(:, :) = qns(:, :) - zcoef1 * sst_m(:, :)
      END WHERE
      !$ACC END KERNELS
    END IF
    DO jc = 1, jncs
      IF (ABS(zfwf(jc) / surf(jncs + 1)) > rsmall) THEN
        !$ACC KERNELS
        zcoef = zfwf(jc) / surf(jc)
        zcoef1 = rcp * zcoef
        WHERE (closea_mask(:, :) == jc)
          emp(:, :) = emp(:, :) - zcoef
          qns(:, :) = qns(:, :) + zcoef1 * sst_m(:, :)
        END WHERE
        !$ACC END KERNELS
      END IF
    END DO
    !$ACC KERNELS
    emp(:, :) = emp(:, :) * tmask(:, :, 1)
    !$ACC END KERNELS
    CALL lbc_lnk(emp, 'T', 1._wp)
  END SUBROUTINE sbc_clo
  SUBROUTINE clo_rnf(p_rnfmsk)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: p_rnfmsk
    !$ACC KERNELS
    IF (jncsr > 0) THEN
      WHERE (closea_mask_rnf(:, :) > 0 .AND. closea_mask(:, :) == 0)
        p_rnfmsk(:, :) = MAX(p_rnfmsk(:, :), 1.0_wp)
      END WHERE
    END IF
    IF (jncse > 0) THEN
      WHERE (closea_mask_empmr(:, :) > 0 .AND. closea_mask(:, :) == 0)
        p_rnfmsk(:, :) = MAX(p_rnfmsk(:, :), 1.0_wp)
      END WHERE
    END IF
    !$ACC END KERNELS
  END SUBROUTINE clo_rnf
  SUBROUTINE clo_bat(k_top, k_bot)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, DIMENSION(:, :), INTENT(INOUT) :: k_top, k_bot
    INTEGER :: inum, id
    INTEGER, DIMENSION(jpi, jpj) :: closea_mask
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdata_in
    TYPE(ProfileData), SAVE :: psy_profile0
    TYPE(ProfileData), SAVE :: psy_profile1
    TYPE(ProfileData), SAVE :: psy_profile2
    TYPE(ProfileData), SAVE :: psy_profile3
    CALL ProfileStart('clo_bat', 'r0', psy_profile0)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'clo_bat : suppression of closed seas'
      WRITE(numout, FMT = *) '~~~~~~~'
    END IF
    CALL ProfileEnd(psy_profile0)
    IF (ln_read_cfg) THEN
      CALL ProfileStart('clo_bat', 'r1', psy_profile1)
      CALL iom_open(cn_domcfg, inum)
      id = iom_varid(inum, 'closea_mask', ldstop = .FALSE.)
      CALL ProfileEnd(psy_profile1)
      IF (id > 0) THEN
        CALL ProfileStart('clo_bat', 'r2', psy_profile2)
        IF (lwp) WRITE(numout, FMT = *) 'Suppressing closed seas in bathymetry based on closea_mask field,'
        CALL iom_get(inum, jpdom_data, 'closea_mask', zdata_in)
        CALL ProfileEnd(psy_profile2)
        !$ACC KERNELS
        closea_mask(:, :) = NINT(zdata_in(:, :))
        WHERE (closea_mask(:, :) > 0)
          k_top(:, :) = 0
          k_bot(:, :) = 0
        END WHERE
        !$ACC END KERNELS
      ELSE
        IF (lwp) WRITE(numout, FMT = *) 'No closea_mask field found in domain_cfg file. No suppression of closed seas.'
      END IF
      CALL iom_close(inum)
    ELSE
      IF (lwp) WRITE(numout, FMT = *) 'No domain_cfg file => no suppression of closed seas.'
    END IF
    CALL ProfileStart('clo_bat', 'r3', psy_profile3)
    l_sbc_clo = .FALSE.
    l_clo_rnf = .FALSE.
    CALL ProfileEnd(psy_profile3)
  END SUBROUTINE clo_bat
END MODULE closea