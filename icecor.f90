MODULE icecor
  USE dom_oce
  USE phycst
  USE ice
  USE ice1D
  USE iceitd
  USE icevar
  USE icectl
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_cor
  CONTAINS
  SUBROUTINE ice_cor(kt, kn)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT(IN) :: kt
    INTEGER, INTENT(IN) :: kn
    INTEGER :: ji, jj, jk, jl
    REAL(KIND = wp) :: zsal, zzc
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zafx
    TYPE(ProfileData), SAVE :: psy_profile0
    TYPE(ProfileData), SAVE :: psy_profile1
    TYPE(ProfileData), SAVE :: psy_profile2
    TYPE(ProfileData), SAVE :: psy_profile3
    TYPE(ProfileData), SAVE :: psy_profile4
    TYPE(ProfileData), SAVE :: psy_profile5
    TYPE(ProfileData), SAVE :: psy_profile6
    CALL ProfileStart('ice_cor', 'r0', psy_profile0)
    IF (ln_timing) CALL timing_start('icecor')
    IF (ln_icediachk) CALL ice_cons_hsm(0, 'icecor', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    IF (kt == nit000 .AND. lwp .AND. kn == 2) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ice_cor:  correct sea ice variables if out of bounds '
      WRITE(numout, FMT = *) '~~~~~~~'
    END IF
    CALL ProfileEnd(psy_profile0)
    !$ACC KERNELS
    WHERE (a_i(:, :, :) >= epsi20)
      h_i(:, :, :) = v_i(:, :, :) / a_i(:, :, :)
    ELSEWHERE
      h_i(:, :, :) = 0._wp
    END WHERE
    WHERE (h_i(:, :, :) < rn_himin) a_i(:, :, :) = a_i(:, :, :) * h_i(:, :, :) / rn_himin
    !$ACC END KERNELS
    CALL ProfileStart('ice_cor', 'r1', psy_profile1)
    at_i(:, :) = SUM(a_i(:, :, :), dim = 3)
    CALL ProfileEnd(psy_profile1)
    !$ACC KERNELS
    DO jl = 1, jpl
      WHERE (at_i(:, :) > rn_amax_2d(:, :)) a_i(:, :, jl) = a_i(:, :, jl) * rn_amax_2d(:, :) / at_i(:, :)
    END DO
    !$ACC END KERNELS
    IF (nn_icesal == 2) THEN
      !$ACC KERNELS
      zzc = rhoi * r1_rdtice
      DO jl = 1, jpl
        DO jj = 1, jpj
          DO ji = 1, jpi
            zsal = sv_i(ji, jj, jl)
            sv_i(ji, jj, jl) = MIN(MAX(rn_simin * v_i(ji, jj, jl), sv_i(ji, jj, jl)), rn_simax * v_i(ji, jj, jl))
            sfx_res(ji, jj) = sfx_res(ji, jj) - (sv_i(ji, jj, jl) - zsal) * zzc
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END IF
    CALL ProfileStart('ice_cor', 'r2', psy_profile2)
    IF (jpl > 1) CALL ice_itd_reb(kt)
    CALL ice_var_zapsmall
    CALL ProfileEnd(psy_profile2)
    IF (kn == 2) THEN
      !CC KERNELS
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          IF (at_i(ji, jj) == 0._wp) THEN
            IF (at_i(ji + 1, jj) == 0._wp) u_ice(ji, jj) = 0._wp
            IF (at_i(ji - 1, jj) == 0._wp) u_ice(ji - 1, jj) = 0._wp
            IF (at_i(ji, jj + 1) == 0._wp) v_ice(ji, jj) = 0._wp
            IF (at_i(ji, jj - 1) == 0._wp) v_ice(ji, jj - 1) = 0._wp
          END IF
        END DO
      END DO
      !CC END KERNELS
      CALL lbc_lnk_multi(u_ice, 'U', - 1., v_ice, 'V', - 1.)
    END IF
    SELECT CASE (kn)
    CASE (1)
      CALL ProfileStart('ice_cor', 'r3', psy_profile3)
      DO jj = 1, jpj
        DO ji = 1, jpi
          diag_heat(ji, jj) = - (SUM(e_i(ji, jj, 1 : nlay_i, :) - e_i_b(ji, jj, 1 : nlay_i, :)) + SUM(e_s(ji, jj, 1 : nlay_s, :) - e_s_b(ji, jj, 1 : nlay_s, :))) * r1_rdtice
          diag_sice(ji, jj) = SUM(sv_i(ji, jj, :) - sv_i_b(ji, jj, :)) * rhoi * r1_rdtice
          diag_vice(ji, jj) = SUM(v_i(ji, jj, :) - v_i_b(ji, jj, :)) * rhoi * r1_rdtice
          diag_vsnw(ji, jj) = SUM(v_s(ji, jj, :) - v_s_b(ji, jj, :)) * rhos * r1_rdtice
        END DO
      END DO
      zafx(:, :) = SUM(a_i(:, :, :) - a_i_b(:, :, :), dim = 3) * r1_rdtice
      CALL ProfileEnd(psy_profile3)
      !$ACC KERNELS
      afx_tot(:, :) = zafx(:, :)
      !$ACC END KERNELS
      IF (iom_use('afxdyn')) CALL iom_put('afxdyn', zafx(:, :))
    CASE (2)
      !$ACC KERNELS
      oa_i(:, :, :) = oa_i(:, :, :) + a_i(:, :, :) * rdt_ice
      !$ACC END KERNELS
      CALL ProfileStart('ice_cor', 'r4', psy_profile4)
      DO jj = 1, jpj
        DO ji = 1, jpi
          diag_heat(ji, jj) = diag_heat(ji, jj) - (SUM(e_i(ji, jj, 1 : nlay_i, :) - e_i_b(ji, jj, 1 : nlay_i, :)) + SUM(e_s(ji, jj, 1 : nlay_s, :) - e_s_b(ji, jj, 1 : nlay_s, :))) * r1_rdtice
          diag_sice(ji, jj) = diag_sice(ji, jj) + SUM(sv_i(ji, jj, :) - sv_i_b(ji, jj, :)) * rhoi * r1_rdtice
          diag_vice(ji, jj) = diag_vice(ji, jj) + SUM(v_i(ji, jj, :) - v_i_b(ji, jj, :)) * rhoi * r1_rdtice
          diag_vsnw(ji, jj) = diag_vsnw(ji, jj) + SUM(v_s(ji, jj, :) - v_s_b(ji, jj, :)) * rhos * r1_rdtice
        END DO
      END DO
      zafx(:, :) = SUM(a_i(:, :, :) - a_i_b(:, :, :), dim = 3) * r1_rdtice
      CALL ProfileEnd(psy_profile4)
      !$ACC KERNELS
      afx_tot(:, :) = afx_tot(:, :) + zafx(:, :)
      !$ACC END KERNELS
      CALL ProfileStart('ice_cor', 'r5', psy_profile5)
      IF (iom_use('afxthd')) CALL iom_put('afxthd', zafx(:, :))
      IF (iom_use('afxtot')) CALL iom_put('afxtot', afx_tot(:, :))
      CALL ProfileEnd(psy_profile5)
    END SELECT
    CALL ProfileStart('ice_cor', 'r6', psy_profile6)
    IF (ln_icediachk) CALL ice_cons_hsm(1, 'icecor', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
    IF (ln_ctl) CALL ice_prt3D('icecor')
    IF (ln_icectl .AND. kn == 2) CALL ice_prt(kt, iiceprt, jiceprt, 2, ' - Final state - ')
    IF (ln_timing) CALL timing_stop('icecor')
    CALL ProfileEnd(psy_profile6)
  END SUBROUTINE ice_cor
END MODULE icecor
