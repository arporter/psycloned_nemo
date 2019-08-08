MODULE dynspg_exp
  USE oce
  USE dom_oce
  USE sbc_oce
  USE phycst
  USE in_out_manager
  USE lib_mpp
  USE lbclnk
  USE prtctl
  USE iom
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dyn_spg_exp
  CONTAINS
  SUBROUTINE dyn_spg_exp(kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    TYPE(ProfileData), SAVE :: psy_profile0
    TYPE(ProfileData), SAVE :: psy_profile1
    IF (kt == nit000) THEN
      CALL ProfileStart('dyn_spg_exp', 'r0', psy_profile0)
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'dyn_spg_exp : surface pressure gradient trend'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~   (explicit free surface)'
      CALL ProfileEnd(psy_profile0)
      !$ACC KERNELS
      spgu(:, :) = 0._wp
      spgv(:, :) = 0._wp
      !$ACC END KERNELS
      CALL ProfileStart('dyn_spg_exp', 'r1', psy_profile1)
      IF (.NOT. ln_linssh .AND. lwp) WRITE(numout, FMT = *) '      non linear free surface: spg is included in dynhpg'
      CALL ProfileEnd(psy_profile1)
    END IF
    IF (ln_linssh) THEN
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          spgu(ji, jj) = - grav * (sshn(ji + 1, jj) - sshn(ji, jj)) * r1_e1u(ji, jj)
          spgv(ji, jj) = - grav * (sshn(ji, jj + 1) - sshn(ji, jj)) * r1_e2v(ji, jj)
        END DO
      END DO
      DO jk = 1, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            ua(ji, jj, jk) = ua(ji, jj, jk) + spgu(ji, jj)
            va(ji, jj, jk) = va(ji, jj, jk) + spgv(ji, jj)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END IF
  END SUBROUTINE dyn_spg_exp
END MODULE dynspg_exp