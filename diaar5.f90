MODULE diaar5
  USE oce
  USE dom_oce
  USE eosbn2
  USE phycst
  USE in_out_manager
  USE zdfddm
  USE zdf_oce
  USE lib_mpp
  USE iom
  USE fldread
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dia_ar5
  PUBLIC :: dia_ar5_alloc
  PUBLIC :: dia_ar5_hst
  REAL(KIND = wp) :: vol0
  REAL(KIND = wp) :: area_tot
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: area
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :) :: thick0
  REAL(KIND = wp), ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: sn0
  LOGICAL :: l_ar5
  CONTAINS
  FUNCTION dia_ar5_alloc()
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER :: dia_ar5_alloc
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('dia_ar5_alloc', 'r0', psy_profile0)
    ALLOCATE(area(jpi, jpj), thick0(jpi, jpj), sn0(jpi, jpj, jpk), STAT = dia_ar5_alloc)
    IF (lk_mpp) CALL mpp_sum(dia_ar5_alloc)
    IF (dia_ar5_alloc /= 0) CALL ctl_warn('dia_ar5_alloc: failed to allocate arrays')
    CALL ProfileEnd(psy_profile0)
  END FUNCTION dia_ar5_alloc
  SUBROUTINE dia_ar5(kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT( IN ) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zvolssh, zvol, zssh_steric, zztmp, zarho, ztemp, zsal, zmass
    REAL(KIND = wp) :: zaw, zbw, zrw
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zarea_ssh, zbotpres
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zpe
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: zrhd, zrhop
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :, :) :: ztsn
    TYPE(ProfileData), SAVE :: psy_profile0
    TYPE(ProfileData), SAVE :: psy_profile1
    TYPE(ProfileData), SAVE :: psy_profile2
    TYPE(ProfileData), SAVE :: psy_profile3
    TYPE(ProfileData), SAVE :: psy_profile4
    TYPE(ProfileData), SAVE :: psy_profile5
    TYPE(ProfileData), SAVE :: psy_profile6
    TYPE(ProfileData), SAVE :: psy_profile7
    TYPE(ProfileData), SAVE :: psy_profile8
    CALL ProfileStart('dia_ar5', 'r0', psy_profile0)
    IF (ln_timing) CALL timing_start('dia_ar5')
    IF (kt == nit000) CALL dia_ar5_init
    CALL ProfileEnd(psy_profile0)
    IF (l_ar5) THEN
      CALL ProfileStart('dia_ar5', 'r1', psy_profile1)
      ALLOCATE(zarea_ssh(jpi, jpj), zbotpres(jpi, jpj))
      ALLOCATE(zrhd(jpi, jpj, jpk), zrhop(jpi, jpj, jpk))
      ALLOCATE(ztsn(jpi, jpj, jpk, jpts))
      CALL ProfileEnd(psy_profile1)
      !$ACC KERNELS
      zarea_ssh(:, :) = area(:, :) * sshn(:, :)
      !$ACC END KERNELS
    END IF
    CALL ProfileStart('dia_ar5', 'r2', psy_profile2)
    IF (iom_use('voltot') .OR. iom_use('sshtot') .OR. iom_use('sshdyn')) THEN
      zvolssh = SUM(zarea_ssh(:, :))
      IF (lk_mpp) CALL mpp_sum(zvolssh)
      zvol = vol0 + zvolssh
      CALL iom_put('voltot', zvol)
      CALL iom_put('sshtot', zvolssh / area_tot)
      CALL iom_put('sshdyn', sshn(:, :) - (zvolssh / area_tot))
    END IF
    CALL ProfileEnd(psy_profile2)
    IF (iom_use('botpres') .OR. iom_use('sshthster') .OR. iom_use('sshsteric')) THEN
      !$ACC KERNELS
      ztsn(:, :, :, jp_tem) = tsn(:, :, :, jp_tem)
      ztsn(:, :, :, jp_sal) = sn0(:, :, :)
      !$ACC END KERNELS
      CALL eos(ztsn, zrhd, gdept_n(:, :, :))
      !$ACC KERNELS
      zbotpres(:, :) = 0._wp
      DO jk = 1, jpkm1
        zbotpres(:, :) = zbotpres(:, :) + e3t_n(:, :, jk) * zrhd(:, :, jk)
      END DO
      !$ACC END KERNELS
      IF (ln_linssh) THEN
        IF (ln_isfcav) THEN
          !$ACC KERNELS
          DO ji = 1, jpi
            DO jj = 1, jpj
              zbotpres(ji, jj) = zbotpres(ji, jj) + sshn(ji, jj) * zrhd(ji, jj, mikt(ji, jj)) + riceload(ji, jj)
            END DO
          END DO
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          zbotpres(:, :) = zbotpres(:, :) + sshn(:, :) * zrhd(:, :, 1)
          !$ACC END KERNELS
        END IF
      END IF
      CALL ProfileStart('dia_ar5', 'r3', psy_profile3)
      zarho = SUM(area(:, :) * zbotpres(:, :))
      IF (lk_mpp) CALL mpp_sum(zarho)
      zssh_steric = - zarho / area_tot
      CALL iom_put('sshthster', zssh_steric)
      CALL eos(tsn, zrhd, zrhop, gdept_n(:, :, :))
      CALL ProfileEnd(psy_profile3)
      !$ACC KERNELS
      zrhop(:, :, jpk) = 0._wp
      !$ACC END KERNELS
      CALL iom_put('rhop', zrhop)
      !$ACC KERNELS
      zbotpres(:, :) = 0._wp
      DO jk = 1, jpkm1
        zbotpres(:, :) = zbotpres(:, :) + e3t_n(:, :, jk) * zrhd(:, :, jk)
      END DO
      !$ACC END KERNELS
      IF (ln_linssh) THEN
        IF (ln_isfcav) THEN
          !$ACC KERNELS
          DO ji = 1, jpi
            DO jj = 1, jpj
              zbotpres(ji, jj) = zbotpres(ji, jj) + sshn(ji, jj) * zrhd(ji, jj, mikt(ji, jj)) + riceload(ji, jj)
            END DO
          END DO
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          zbotpres(:, :) = zbotpres(:, :) + sshn(:, :) * zrhd(:, :, 1)
          !$ACC END KERNELS
        END IF
      END IF
      CALL ProfileStart('dia_ar5', 'r4', psy_profile4)
      zarho = SUM(area(:, :) * zbotpres(:, :))
      IF (lk_mpp) CALL mpp_sum(zarho)
      zssh_steric = - zarho / area_tot
      CALL iom_put('sshsteric', zssh_steric)
      CALL ProfileEnd(psy_profile4)
      !$ACC KERNELS
      zztmp = rau0 * grav * 1.E-4_wp
      zbotpres(:, :) = zztmp * (zbotpres(:, :) + sshn(:, :) + thick0(:, :))
      !$ACC END KERNELS
      CALL iom_put('botpres', zbotpres)
    END IF
    IF (iom_use('masstot') .OR. iom_use('temptot') .OR. iom_use('saltot')) THEN
      !$ACC KERNELS
      ztemp = 0._wp
      zsal = 0._wp
      DO jk = 1, jpkm1
        DO jj = 1, jpj
          DO ji = 1, jpi
            zztmp = area(ji, jj) * e3t_n(ji, jj, jk)
            ztemp = ztemp + zztmp * tsn(ji, jj, jk, jp_tem)
            zsal = zsal + zztmp * tsn(ji, jj, jk, jp_sal)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      IF (ln_linssh) THEN
        IF (ln_isfcav) THEN
          !$ACC KERNELS
          DO ji = 1, jpi
            DO jj = 1, jpj
              ztemp = ztemp + zarea_ssh(ji, jj) * tsn(ji, jj, mikt(ji, jj), jp_tem)
              zsal = zsal + zarea_ssh(ji, jj) * tsn(ji, jj, mikt(ji, jj), jp_sal)
            END DO
          END DO
          !$ACC END KERNELS
        ELSE
          CALL ProfileStart('dia_ar5', 'r5', psy_profile5)
          ztemp = ztemp + SUM(zarea_ssh(:, :) * tsn(:, :, 1, jp_tem))
          zsal = zsal + SUM(zarea_ssh(:, :) * tsn(:, :, 1, jp_sal))
          CALL ProfileEnd(psy_profile5)
        END IF
      END IF
      CALL ProfileStart('dia_ar5', 'r6', psy_profile6)
      IF (lk_mpp) THEN
        CALL mpp_sum(ztemp)
        CALL mpp_sum(zsal)
      END IF
      zmass = rau0 * (zarho + zvol)
      ztemp = ztemp / zvol
      zsal = zsal / zvol
      CALL iom_put('masstot', zmass)
      CALL iom_put('temptot', ztemp)
      CALL iom_put('saltot', zsal)
      CALL ProfileEnd(psy_profile6)
    END IF
    IF (iom_use('tnpeo')) THEN
      CALL ProfileStart('dia_ar5', 'r7', psy_profile7)
      ALLOCATE(zpe(jpi, jpj))
      CALL ProfileEnd(psy_profile7)
      !$ACC KERNELS
      zpe(:, :) = 0._wp
      !$ACC END KERNELS
      IF (ln_zdfddm) THEN
        !$ACC KERNELS
        DO jk = 2, jpk
          DO jj = 1, jpj
            DO ji = 1, jpi
              IF (rn2(ji, jj, jk) > 0._wp) THEN
                zrw = (gdepw_n(ji, jj, jk) - gdept_n(ji, jj, jk)) / (gdept_n(ji, jj, jk - 1) - gdept_n(ji, jj, jk))
                zaw = rab_n(ji, jj, jk, jp_tem) * (1. - zrw) + rab_n(ji, jj, jk - 1, jp_tem) * zrw
                zbw = rab_n(ji, jj, jk, jp_sal) * (1. - zrw) + rab_n(ji, jj, jk - 1, jp_sal) * zrw
                zpe(ji, jj) = zpe(ji, jj) - grav * (avt(ji, jj, jk) * zaw * (tsn(ji, jj, jk - 1, jp_tem) - tsn(ji, jj, jk, jp_tem)) - avs(ji, jj, jk) * zbw * (tsn(ji, jj, jk - 1, jp_sal) - tsn(ji, jj, jk, jp_sal)))
              END IF
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      ELSE
        !$ACC KERNELS
        DO jk = 1, jpk
          DO ji = 1, jpi
            DO jj = 1, jpj
              zpe(ji, jj) = zpe(ji, jj) + avt(ji, jj, jk) * MIN(0._wp, rn2(ji, jj, jk)) * rau0 * e3w_n(ji, jj, jk)
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      END IF
      CALL iom_put('tnpeo', zpe)
      DEALLOCATE(zpe)
    END IF
    CALL ProfileStart('dia_ar5', 'r8', psy_profile8)
    IF (l_ar5) THEN
      DEALLOCATE(zarea_ssh, zbotpres)
      DEALLOCATE(zrhd, zrhop)
      DEALLOCATE(ztsn)
    END IF
    IF (ln_timing) CALL timing_stop('dia_ar5')
    CALL ProfileEnd(psy_profile8)
  END SUBROUTINE dia_ar5
  SUBROUTINE dia_ar5_hst(ktra, cptr, pua, pva)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT(IN ) :: ktra
    CHARACTER(LEN = 3), INTENT(IN) :: cptr
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pua
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN) :: pva
    INTEGER :: ji, jj, jk
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z2d
    TYPE(ProfileData), SAVE :: psy_profile0
    TYPE(ProfileData), SAVE :: psy_profile1
    !$ACC KERNELS
    z2d(:, :) = pua(:, :, 1)
    DO jk = 1, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          z2d(ji, jj) = z2d(ji, jj) + pua(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL ProfileStart('dia_ar5_hst', 'r0', psy_profile0)
    CALL lbc_lnk(z2d, 'U', - 1.)
    IF (cptr == 'adv') THEN
      IF (ktra == jp_tem) CALL iom_put("uadv_heattr", rau0_rcp * z2d)
      IF (ktra == jp_sal) CALL iom_put("uadv_salttr", rau0 * z2d)
    END IF
    IF (cptr == 'ldf') THEN
      IF (ktra == jp_tem) CALL iom_put("udiff_heattr", rau0_rcp * z2d)
      IF (ktra == jp_sal) CALL iom_put("udiff_salttr", rau0 * z2d)
    END IF
    CALL ProfileEnd(psy_profile0)
    !$ACC KERNELS
    z2d(:, :) = pva(:, :, 1)
    DO jk = 1, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          z2d(ji, jj) = z2d(ji, jj) + pva(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL ProfileStart('dia_ar5_hst', 'r1', psy_profile1)
    CALL lbc_lnk(z2d, 'V', - 1.)
    IF (cptr == 'adv') THEN
      IF (ktra == jp_tem) CALL iom_put("vadv_heattr", rau0_rcp * z2d)
      IF (ktra == jp_sal) CALL iom_put("vadv_salttr", rau0 * z2d)
    END IF
    IF (cptr == 'ldf') THEN
      IF (ktra == jp_tem) CALL iom_put("vdiff_heattr", rau0_rcp * z2d)
      IF (ktra == jp_sal) CALL iom_put("vdiff_salttr", rau0 * z2d)
    END IF
    CALL ProfileEnd(psy_profile1)
  END SUBROUTINE dia_ar5_hst
  SUBROUTINE dia_ar5_init
    INTEGER :: inum
    INTEGER :: ik
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zztmp
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :, :) :: zsaldta
    l_ar5 = .FALSE.
    IF (iom_use('voltot') .OR. iom_use('sshtot') .OR. iom_use('sshdyn') .OR. iom_use('masstot') .OR. iom_use('temptot') .OR. iom_use('saltot') .OR. iom_use('botpres') .OR. iom_use('sshthster') .OR. iom_use('sshsteric')) L_ar5 = .TRUE.
    IF (l_ar5) THEN
      IF (dia_ar5_alloc() /= 0) CALL ctl_stop('STOP', 'dia_ar5_init : unable to allocate arrays')
      !$ACC KERNELS
      area(:, :) = e1e2t(:, :) * tmask_i(:, :)
      !$ACC END KERNELS
      area_tot = SUM(area(:, :))
      IF (lk_mpp) CALL mpp_sum(area_tot)
      !$ACC KERNELS
      vol0 = 0._wp
      thick0(:, :) = 0._wp
      !$ACC END KERNELS
      DO jk = 1, jpkm1
        vol0 = vol0 + SUM(area(:, :) * tmask(:, :, jk) * e3t_0(:, :, jk))
        !$ACC KERNELS
        thick0(:, :) = thick0(:, :) + tmask_i(:, :) * tmask(:, :, jk) * e3t_0(:, :, jk)
        !$ACC END KERNELS
      END DO
      IF (lk_mpp) CALL mpp_sum(vol0)
      IF (iom_use('sshthster')) THEN
        ALLOCATE(zsaldta(jpi, jpj, jpj, jpts))
        CALL iom_open('sali_ref_clim_monthly', inum)
        CALL iom_get(inum, jpdom_data, 'vosaline', zsaldta(:, :, :, 1), 1)
        CALL iom_get(inum, jpdom_data, 'vosaline', zsaldta(:, :, :, 2), 12)
        CALL iom_close(inum)
        !$ACC KERNELS
        sn0(:, :, :) = 0.5_wp * (zsaldta(:, :, :, 1) + zsaldta(:, :, :, 2))
        sn0(:, :, :) = sn0(:, :, :) * tmask(:, :, :)
        !$ACC END KERNELS
        IF (ln_zps) THEN
          !$ACC KERNELS
          DO jj = 1, jpj
            DO ji = 1, jpi
              ik = mbkt(ji, jj)
              IF (ik > 1) THEN
                zztmp = (gdept_1d(ik) - gdept_0(ji, jj, ik)) / (gdept_1d(ik) - gdept_1d(ik - 1))
                sn0(ji, jj, ik) = (1._wp - zztmp) * sn0(ji, jj, ik) + zztmp * sn0(ji, jj, ik - 1)
              END IF
            END DO
          END DO
          !$ACC END KERNELS
        END IF
        DEALLOCATE(zsaldta)
      END IF
    END IF
  END SUBROUTINE dia_ar5_init
END MODULE diaar5