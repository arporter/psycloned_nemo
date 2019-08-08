MODULE dynkeg
  USE oce
  USE dom_oce
  USE trd_oce
  USE trddyn
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  USE prtctl
  USE timing
  USE bdy_oce
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_keg
  INTEGER, PARAMETER, PUBLIC :: nkeg_C2 = 0
  INTEGER, PARAMETER, PUBLIC :: nkeg_HW = 1
  REAL(KIND = wp) :: r1_48 = 1._wp / 48._wp
  CONTAINS
  SUBROUTINE dyn_keg(kt, kscheme)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT( IN ) :: kt
    INTEGER, INTENT( IN ) :: kscheme
    INTEGER :: ji, jj, jk, jb
    INTEGER :: ii, ifu, ib_bdy
    INTEGER :: ij, ifv, igrd
    REAL(KIND = wp) :: zu, zv
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zhke
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: ztrdu, ztrdv
    TYPE(ProfileData), SAVE :: psy_profile0
    TYPE(ProfileData), SAVE :: psy_profile1
    TYPE(ProfileData), SAVE :: psy_profile2
    TYPE(ProfileData), SAVE :: psy_profile3
    CALL ProfileStart('dyn_keg', 'r0', psy_profile0)
    IF (ln_timing) CALL timing_start('dyn_keg')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn_keg : kinetic energy gradient trend, scheme number=', kscheme
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~'
    END IF
    CALL ProfileEnd(psy_profile0)
    IF (l_trddyn) THEN
      CALL ProfileStart('dyn_keg', 'r1', psy_profile1)
      ALLOCATE(ztrdu(jpi, jpj, jpk), ztrdv(jpi, jpj, jpk))
      CALL ProfileEnd(psy_profile1)
      !$ACC KERNELS
      ztrdu(:, :, :) = ua(:, :, :)
      ztrdv(:, :, :) = va(:, :, :)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    zhke(:, :, jpk) = 0._wp
    !$ACC END KERNELS
    CALL ProfileStart('dyn_keg', 'r2', psy_profile2)
    IF (ln_bdy) THEN
      DO ib_bdy = 1, nb_bdy
        IF (cn_dyn3d(ib_bdy) /= 'none') THEN
          igrd = 2
          DO jb = 1, idx_bdy(ib_bdy) % nblenrim(igrd)
            DO jk = 1, jpkm1
              ii = idx_bdy(ib_bdy) % nbi(jb, igrd)
              ij = idx_bdy(ib_bdy) % nbj(jb, igrd)
              ifu = NINT(idx_bdy(ib_bdy) % flagu(jb, igrd))
              un(ii - ifu, ij, jk) = un(ii, ij, jk) * umask(ii, ij, jk)
            END DO
          END DO
          igrd = 3
          DO jb = 1, idx_bdy(ib_bdy) % nblenrim(igrd)
            DO jk = 1, jpkm1
              ii = idx_bdy(ib_bdy) % nbi(jb, igrd)
              ij = idx_bdy(ib_bdy) % nbj(jb, igrd)
              ifv = NINT(idx_bdy(ib_bdy) % flagv(jb, igrd))
              vn(ii, ij - ifv, jk) = vn(ii, ij, jk) * vmask(ii, ij, jk)
            END DO
          END DO
        END IF
      END DO
    END IF
    CALL ProfileEnd(psy_profile2)
    SELECT CASE (kscheme)
    CASE (nkeg_C2)
      !$ACC KERNELS
      DO jk = 1, jpkm1
        DO jj = 2, jpj
          DO ji = 2, jpi
            zu = un(ji - 1, jj, jk) * un(ji - 1, jj, jk) + un(ji, jj, jk) * un(ji, jj, jk)
            zv = vn(ji, jj - 1, jk) * vn(ji, jj - 1, jk) + vn(ji, jj, jk) * vn(ji, jj, jk)
            zhke(ji, jj, jk) = 0.25_wp * (zv + zu)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    CASE (nkeg_HW)
      !$ACC KERNELS
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zu = 8._wp * (un(ji - 1, jj, jk) * un(ji - 1, jj, jk) + un(ji, jj, jk) * un(ji, jj, jk)) + (un(ji - 1, jj - 1, jk) + un(ji - 1, jj + 1, jk)) * (un(ji - 1, jj - 1, jk) + un(ji - 1, jj + 1, jk)) + (un(ji, jj - 1, jk) + un(ji, jj + 1, jk)) * (un(ji, jj - 1, jk) + un(ji, jj + 1, jk))
            zv = 8._wp * (vn(ji, jj - 1, jk) * vn(ji, jj - 1, jk) + vn(ji, jj, jk) * vn(ji, jj, jk)) + (vn(ji - 1, jj - 1, jk) + vn(ji + 1, jj - 1, jk)) * (vn(ji - 1, jj - 1, jk) + vn(ji + 1, jj - 1, jk)) + (vn(ji - 1, jj, jk) + vn(ji + 1, jj, jk)) * (vn(ji - 1, jj, jk) + vn(ji + 1, jj, jk))
            zhke(ji, jj, jk) = r1_48 * (zv + zu)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk(zhke, 'T', 1.)
    END SELECT
    IF (ln_bdy) THEN
      !$ACC KERNELS
      un(:, :, :) = un(:, :, :) * umask(:, :, :)
      vn(:, :, :) = vn(:, :, :) * vmask(:, :, :)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    DO jk = 1, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          ua(ji, jj, jk) = ua(ji, jj, jk) - (zhke(ji + 1, jj, jk) - zhke(ji, jj, jk)) / e1u(ji, jj)
          va(ji, jj, jk) = va(ji, jj, jk) - (zhke(ji, jj + 1, jk) - zhke(ji, jj, jk)) / e2v(ji, jj)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    IF (l_trddyn) THEN
      !$ACC KERNELS
      ztrdu(:, :, :) = ua(:, :, :) - ztrdu(:, :, :)
      ztrdv(:, :, :) = va(:, :, :) - ztrdv(:, :, :)
      !$ACC END KERNELS
      CALL trd_dyn(ztrdu, ztrdv, jpdyn_keg, kt)
      DEALLOCATE(ztrdu, ztrdv)
    END IF
    CALL ProfileStart('dyn_keg', 'r3', psy_profile3)
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = ua, clinfo1 = ' keg  - Ua: ', mask1 = umask, tab3d_2 = va, clinfo2 = ' Va: ', mask2 = vmask, clinfo3 = 'dyn')
    IF (ln_timing) CALL timing_stop('dyn_keg')
    CALL ProfileEnd(psy_profile3)
  END SUBROUTINE dyn_keg
END MODULE dynkeg