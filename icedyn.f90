MODULE icedyn
  USE phycst
  USE dom_oce
  USE ice
  USE icedyn_rhg
  USE icedyn_adv
  USE icedyn_rdgrft
  USE icecor
  USE icevar
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_dyn
  PUBLIC :: ice_dyn_init
  INTEGER :: nice_dyn
  INTEGER, PARAMETER :: np_dynFULL = 1
  INTEGER, PARAMETER :: np_dynRHGADV = 2
  INTEGER, PARAMETER :: np_dynADV = 3
  LOGICAL :: ln_dynFULL
  LOGICAL :: ln_dynRHGADV
  LOGICAL :: ln_dynADV
  REAL(KIND = wp) :: rn_uice
  REAL(KIND = wp) :: rn_vice
  CONTAINS
  SUBROUTINE ice_dyn(kt)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jl
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: zhmax
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data2
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data3
    CALL profile_psy_data0 % PreStart('ice_dyn', 'r0', 0, 0)
    IF (ln_timing) CALL timing_start('icedyn')
    IF (kt == nit000 .AND. lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ice_dyn: sea-ice dynamics'
      WRITE(numout, FMT = *) '~~~~~~~'
    END IF
    CALL profile_psy_data0 % PostEnd
    IF (ln_landfast) THEN
      !$ACC KERNELS
      tau_icebfr(:, :) = 0._wp
      DO jl = 1, jpl
        WHERE (h_i(:, :, jl) > ht_n(:, :) * rn_gamma) tau_icebfr(:, :) = tau_icebfr(:, :) + a_i(:, :, jl) * rn_icebfr
      END DO
      !$ACC END KERNELS
      IF (iom_use('tau_icebfr')) CALL iom_put('tau_icebfr', tau_icebfr)
    END IF
    !$ACC KERNELS
    zhmax(:, :, :) = h_i_b(:, :, :)
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('ice_dyn', 'r1', 0, 0)
    DO jl = 1, jpl
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zhmax(ji, jj, jl) = MAX(epsi20, MAXVAL(h_i_b(ji - 1 : ji + 1, jj - 1 : jj + 1, jl)))
        END DO
      END DO
    END DO
    CALL lbc_lnk(zhmax(:, :, :), 'T', 1.)
    CALL profile_psy_data1 % PostEnd
    SELECT CASE (nice_dyn)
    CASE (np_dynFULL)
      CALL profile_psy_data2 % PreStart('ice_dyn', 'r2', 0, 0)
      CALL ice_dyn_rhg(kt)
      CALL ice_dyn_adv(kt)
      CALL hbig(zhmax)
      CALL ice_dyn_rdgrft(kt)
      CALL ice_cor(kt, 1)
      CALL profile_psy_data2 % PostEnd
    CASE (np_dynRHGADV)
      CALL profile_psy_data3 % PreStart('ice_dyn', 'r3', 0, 0)
      CALL ice_dyn_rhg(kt)
      CALL ice_dyn_adv(kt)
      CALL Hpiling
      CALL profile_psy_data3 % PostEnd
    CASE (np_dynADV)
      !$ACC KERNELS
      u_ice(:, :) = rn_uice * umask(:, :, 1)
      v_ice(:, :) = rn_vice * vmask(:, :, 1)
      !$ACC END KERNELS
      CALL ice_dyn_adv(kt)
    END SELECT
    IF (ln_timing) CALL timing_stop('icedyn')
  END SUBROUTINE ice_dyn
  SUBROUTINE Hbig(phmax)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: phmax
    INTEGER :: ji, jj, jl
    REAL(KIND = wp) :: zh, zdv
    CALL ice_var_zapsmall
    !$ACC KERNELS
    DO jl = 1, jpl
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1, jpj
        DO ji = 1, jpi
          IF (v_i(ji, jj, jl) > 0._wp) THEN
            zh = v_i(ji, jj, jl) / a_i(ji, jj, jl)
            zdv = v_i(ji, jj, jl) - v_i_b(ji, jj, jl)
            IF ((zdv > 0.0 .AND. zh > phmax(ji, jj, jl) .AND. at_i_b(ji, jj) < 0.80) .OR. (zdv <= 0.0 .AND. zh > phmax(ji, jj, jl))) THEN
              a_i(ji, jj, jl) = v_i(ji, jj, jl) / MIN(phmax(ji, jj, jl), hi_max(jpl))
            END IF
          END IF
        END DO
      END DO
    END DO
    WHERE (a_ip(:, :, :) > a_i(:, :, :)) a_ip(:, :, :) = a_i(:, :, :)
    !$ACC END KERNELS
  END SUBROUTINE Hbig
  SUBROUTINE Hpiling
    INTEGER :: jl
    CALL ice_var_zapsmall
    !$ACC KERNELS
    at_i(:, :) = SUM(a_i(:, :, :), dim = 3)
    DO jl = 1, jpl
      WHERE (at_i(:, :) > epsi20)
        a_i(:, :, jl) = a_i(:, :, jl) * (1._wp + MIN(rn_amax_2d(:, :) - at_i(:, :), 0._wp) / at_i(:, :))
      END WHERE
    END DO
    !$ACC END KERNELS
  END SUBROUTINE Hpiling
  SUBROUTINE ice_dyn_init
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ios, ioptio
    NAMELIST /namdyn/ ln_dynFULL, ln_dynRHGADV, ln_dynADV, rn_uice, rn_vice, rn_ishlat, ln_landfast, rn_gamma, rn_icebfr, rn_lfrelax
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    CALL profile_psy_data0 % PreStart('ice_dyn_init', 'r0', 0, 0)
    REWIND(UNIT = numnam_ice_ref)
    READ(numnam_ice_ref, namdyn, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'namdyn in reference namelist', lwp)
    REWIND(UNIT = numnam_ice_cfg)
    READ(numnam_ice_cfg, namdyn, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'namdyn in configuration namelist', lwp)
    IF (lwm) WRITE(numoni, namdyn)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'ice_dyn_init: ice parameters for ice dynamics '
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist namdyn:'
      WRITE(numout, FMT = *) '      Full ice dynamics      (rhg + adv + ridge/raft + corr)  ln_dynFULL   = ', ln_dynFULL
      WRITE(numout, FMT = *) '      No ridge/raft & No cor (rhg + adv)                      ln_dynRHGADV = ', ln_dynRHGADV
      WRITE(numout, FMT = *) '      Advection only         (rn_uvice + adv)                 ln_dynADV    = ', ln_dynADV
      WRITE(numout, FMT = *) '         with prescribed velocity given by   (u,v)_ice = (rn_uice,rn_vice) = (', rn_uice, ',', rn_vice, ')'
      WRITE(numout, FMT = *) '      lateral boundary condition for sea ice dynamics         rn_ishlat    = ', rn_ishlat
      WRITE(numout, FMT = *) '      Landfast: param (T or F)                                ln_landfast  = ', ln_landfast
      WRITE(numout, FMT = *) '         fraction of ocean depth that ice must reach          rn_gamma     = ', rn_gamma
      WRITE(numout, FMT = *) '         maximum bottom stress per unit area of contact       rn_icebfr    = ', rn_icebfr
      WRITE(numout, FMT = *) '         relax time scale (s-1) to reach static friction      rn_lfrelax   = ', rn_lfrelax
      WRITE(numout, FMT = *)
    END IF
    ioptio = 0
    IF (ln_dynfull) THEN
      ioptio = ioptio + 1
      nice_dyn = np_dynfull
    END IF
    IF (ln_dynrhgadv) THEN
      ioptio = ioptio + 1
      nice_dyn = np_dynrhgadv
    END IF
    IF (ln_dynadv) THEN
      ioptio = ioptio + 1
      nice_dyn = np_dynadv
    END IF
    IF (ioptio /= 1) CALL ctl_stop('ice_dyn_init: one and only one ice dynamics option has to be defined ')
    IF (rn_ishlat == 0.) THEN
      IF (lwp) WRITE(numout, FMT = *) '   ===>>>   ice lateral  free-slip'
    ELSE IF (rn_ishlat == 2.) THEN
      IF (lwp) WRITE(numout, FMT = *) '   ===>>>   ice lateral  no-slip'
    ELSE IF (0. < rn_ishlat .AND. rn_ishlat < 2.) THEN
      IF (lwp) WRITE(numout, FMT = *) '   ===>>>   ice lateral  partial-slip'
    ELSE IF (2. < rn_ishlat) THEN
      IF (lwp) WRITE(numout, FMT = *) '   ===>>>   ice lateral  strong-slip'
    END IF
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    IF (.NOT. ln_landfast) tau_icebfr(:, :) = 0._wp
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('ice_dyn_init', 'r1', 0, 0)
    CALL ice_dyn_rdgrft_init
    CALL ice_dyn_rhg_init
    CALL ice_dyn_adv_init
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE ice_dyn_init
END MODULE icedyn