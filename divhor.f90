MODULE divhor
  USE oce
  USE dom_oce
  USE sbc_oce, ONLY: ln_rnf, ln_isf
  USE sbcrnf
  USE sbcisf
  USE iscplhsb
  USE iscplini
  USE in_out_manager
  USE lbclnk
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: div_hor
  CONTAINS
  SUBROUTINE div_hor(kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk
    REAL(KIND = wp) :: zraur, zdep
    TYPE(ProfileData), SAVE :: psy_profile0
    TYPE(ProfileData), SAVE :: psy_profile1
    CALL ProfileStart('div_hor', 'r0', psy_profile0)
    IF (ln_timing) CALL timing_start('div_hor')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'div_hor : horizontal velocity divergence '
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~   '
    END IF
    CALL ProfileEnd(psy_profile0)
    !$ACC KERNELS
    DO jk = 1, jpkm1
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          hdivn(ji, jj, jk) = (e2u(ji, jj) * e3u_n(ji, jj, jk) * un(ji, jj, jk) - e2u(ji - 1, jj) * e3u_n(ji - 1, jj, jk) * un(ji - 1, jj, jk) + e1v(ji, jj) * e3v_n(ji, jj, jk) * vn(ji, jj, jk) - e1v(ji, jj - 1) * e3v_n(ji, jj - 1, jk) * vn(ji, jj - 1, jk)) * r1_e1e2t(ji, jj) / e3t_n(ji, jj, jk)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL ProfileStart('div_hor', 'r1', psy_profile1)
    IF (ln_rnf) CALL sbc_rnf_div(hdivn)
    IF (ln_isf) CALL sbc_isf_div(hdivn)
    IF (ln_iscpl .AND. ln_hsb) CALL iscpl_div(hdivn)
    CALL lbc_lnk(hdivn, 'T', 1.)
    IF (ln_timing) CALL timing_stop('div_hor')
    CALL ProfileEnd(psy_profile1)
  END SUBROUTINE div_hor
END MODULE divhor