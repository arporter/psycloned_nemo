MODULE icbclv
  USE par_oce
  USE dom_oce
  USE phycst
  USE lib_mpp
  USE lbclnk
  USE icb_oce
  USE icbdia
  USE icbutl
  USE icb_oce
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: icb_clv_flx
  PUBLIC :: icb_clv
  CONTAINS
  SUBROUTINE icb_clv_flx(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    REAL(KIND = wp) :: zcalving_used, zdist, zfact
    INTEGER :: jn, ji, jj
    INTEGER :: imx
    LOGICAL, SAVE :: ll_first_call = .TRUE.
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_clv_flx', 'r0', 0, 0)
    zfact = ((1000._wp) ** 3 / (NINT(rday) * nyear_len(1))) * rn_rho_bergs
    berg_grid % calving(:, :) = src_calving(:, :) * zfact * tmask_i(:, :) * tmask(:, :, 1)
    berg_grid % calving_hflx(:, :) = src_calving_hflx(:, :) * tmask_i(:, :) * tmask(:, :, 1)
    IF (ll_first_call .AND. .NOT. l_restarted_bergs) THEN
      ll_first_call = .FALSE.
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          IF (berg_grid % calving(ji, jj) /= 0._wp) berg_grid % stored_heat(ji, jj) = SUM(berg_grid % stored_ice(ji, jj, :)) * &
&berg_grid % calving_hflx(ji, jj) * e1e2t(ji, jj) / berg_grid % calving(ji, jj)
        END DO
      END DO
    END IF
    DO jj = 1, jpj
      DO ji = 1, jpi
        imx = berg_grid % maxclass(ji, jj)
        zdist = SUM(rn_distribution(1 : nclasses)) / SUM(rn_distribution(1 : imx))
        DO jn = 1, imx
          berg_grid % stored_ice(ji, jj, jn) = berg_grid % stored_ice(ji, jj, jn) + berg_dt * berg_grid % calving(ji, jj) * &
&rn_distribution(jn) * zdist
        END DO
      END DO
    END DO
    zcalving_used = SUM(berg_grid % calving(:, :))
    berg_grid % tmp(:, :) = berg_dt * berg_grid % calving_hflx(:, :) * e1e2t(:, :) * tmask_i(:, :)
    berg_grid % stored_heat(:, :) = berg_grid % stored_heat(:, :) + berg_grid % tmp(:, :)
    CALL icb_dia_income(kt, zcalving_used, berg_grid % tmp)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_clv_flx
  SUBROUTINE icb_clv(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jn
    INTEGER :: icnt, icntmax
    TYPE(iceberg) :: newberg
    TYPE(point) :: newpt
    REAL(KIND = wp) :: zday, zcalved_to_berg, zheat_to_berg
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('icb_clv', 'r0', 0, 0)
    icntmax = 0
    zday = REAL(nday_year, wp) + REAL(nsec_day, wp) / 86400.0_wp
    DO jn = 1, nclasses
      DO jj = nicbdj, nicbej
        DO ji = nicbdi, nicbei
          icnt = 0
          DO WHILE (berg_grid % stored_ice(ji, jj, jn) >= rn_initial_mass(jn) * rn_mass_scaling(jn))
            newpt % lon = glamt(ji, jj)
            newpt % lat = gphit(ji, jj)
            newpt % xi = REAL(mig(ji), wp)
            newpt % yj = REAL(mjg(jj), wp)
            newpt % uvel = 0._wp
            newpt % vvel = 0._wp
            newpt % mass = rn_initial_mass(jn)
            newpt % thickness = rn_initial_thickness(jn)
            newpt % width = first_width(jn)
            newpt % length = first_length(jn)
            newberg % mass_scaling = rn_mass_scaling(jn)
            newpt % mass_of_bits = 0._wp
            newpt % year = nyear
            newpt % day = zday
            newpt % heat_density = berg_grid % stored_heat(ji, jj) / berg_grid % stored_ice(ji, jj, jn)
            CALL icb_utl_incr
            newberg % number(:) = num_bergs(:)
            CALL icb_utl_add(newberg, newpt)
            zcalved_to_berg = rn_initial_mass(jn) * rn_mass_scaling(jn)
            zheat_to_berg = zcalved_to_berg * newpt % heat_density
            berg_grid % stored_heat(ji, jj) = berg_grid % stored_heat(ji, jj) - zheat_to_berg
            berg_grid % stored_ice(ji, jj, jn) = berg_grid % stored_ice(ji, jj, jn) - zcalved_to_berg
            icnt = icnt + 1
            CALL icb_dia_calve(ji, jj, jn, zcalved_to_berg, zheat_to_berg)
          END DO
          icntmax = MAX(icntmax, icnt)
        END DO
      END DO
    END DO
    DO jn = 1, nclasses
      CALL lbc_lnk('icbclv', berg_grid % stored_ice(:, :, jn), 'T', 1._wp)
    END DO
    CALL lbc_lnk('icbclv', berg_grid % stored_heat, 'T', 1._wp)
    IF (nn_verbose_level > 0 .AND. icntmax > 1) WRITE(numicb, FMT = *) 'icb_clv: icnt=', icnt, ' on', narea
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE icb_clv
END MODULE icbclv