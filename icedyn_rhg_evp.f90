MODULE icedyn_rhg_evp
  USE phycst
  USE dom_oce
  USE sbc_oce, ONLY: ln_ice_embd, nn_fsbc, ssh_m
  USE sbc_ice, ONLY: utau_ice, vtau_ice, snwice_mass, snwice_mass_b
  USE ice
  USE icedyn_rdgrft
  USE bdy_oce, ONLY: ln_bdy
  USE bdyice
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  USE prtctl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_dyn_rhg_evp
  PUBLIC :: rhg_evp_rst
  CONTAINS
  SUBROUTINE ice_dyn_rhg_evp(kt, pstress1_i, pstress2_i, pstress12_i, pshear_i, pdivu_i, pdelta_i)
    INTEGER, INTENT(IN   ) :: kt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pstress1_i, pstress2_i, pstress12_i
    REAL(KIND = wp), DIMENSION(:, :), INTENT(  OUT) :: pshear_i, pdivu_i, pdelta_i
    INTEGER :: ji, jj
    INTEGER :: jter
    REAL(KIND = wp) :: zrhoco
    REAL(KIND = wp) :: zdtevp, z1_dtevp
    REAL(KIND = wp) :: ecc2, z1_ecc2
    REAL(KIND = wp) :: zalph1, z1_alph1, zalph2, z1_alph2
    REAL(KIND = wp) :: zm1, zm2, zm3, zmassU, zmassV
    REAL(KIND = wp) :: zdelta, zp_delf, zds2, zdt, zdt2, zdiv, zdiv2
    REAL(KIND = wp) :: zTauO, zTauB, zTauE, zvel
    REAL(KIND = wp) :: zresm
    REAL(KIND = wp) :: zintb, zintn
    REAL(KIND = wp) :: zfac_x, zfac_y
    REAL(KIND = wp) :: zshear, zdum1, zdum2
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z1_e1t0, z1_e2t0
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zp_delt
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zbeta
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zdt_m
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zaU, zaV
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zmU_t, zmV_t
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zmf
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zTauU_ia, ztauV_ia
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zspgU, zspgV
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: v_oceU, u_oceV, v_iceU, u_iceV
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zfU, zfV
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zds
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zs1, zs2, zs12
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zu_ice, zv_ice, zresr
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zpice
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zCorx, zCory
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztaux_oi, ztauy_oi
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zswitchU, zswitchV
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zmaskU, zmaskV
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zfmask, zwf
    REAL(KIND = wp), PARAMETER :: zepsi = 1.0E-20_wp
    REAL(KIND = wp), PARAMETER :: zmmin = 1._wp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zswi
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zsig1, zsig2, zsig3
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_sig1
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_sig2
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_dssh_dx
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_dssh_dy
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_corstrx
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_corstry
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_intstrx
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_intstry
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_utau_oi
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_vtau_oi
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_xmtrp_ice
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_ymtrp_ice
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_xmtrp_snw
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_ymtrp_snw
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_xatrp
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zdiag_yatrp
    IF (kt == nit000 .AND. lwp) WRITE(numout, FMT = *) '-- ice_dyn_rhg_evp: EVP sea-ice rheology'
    !$ACC KERNELS
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        zfmask(ji, jj) = tmask(ji, jj, 1) * tmask(ji + 1, jj, 1) * tmask(ji, jj + 1, 1) * tmask(ji + 1, jj + 1, 1)
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk(zfmask, 'F', 1._wp)
    !$ACC KERNELS
    zwf(:, :) = zfmask(:, :)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        IF (zfmask(ji, jj) == 0._wp) THEN
          zfmask(ji, jj) = rn_ishlat * MIN(1._wp, MAX(zwf(ji + 1, jj), zwf(ji, jj + 1), zwf(ji - 1, jj), zwf(ji, jj - 1)))
        END IF
      END DO
    END DO
    DO jj = 2, jpjm1
      IF (zfmask(1, jj) == 0._wp) THEN
        zfmask(1, jj) = rn_ishlat * MIN(1._wp, MAX(zwf(2, jj), zwf(1, jj + 1), zwf(1, jj - 1)))
      END IF
      IF (zfmask(jpi, jj) == 0._wp) THEN
        zfmask(jpi, jj) = rn_ishlat * MIN(1._wp, MAX(zwf(jpi, jj + 1), zwf(jpim1, jj), zwf(jpi, jj - 1)))
      END IF
    END DO
    DO ji = 2, jpim1
      IF (zfmask(ji, 1) == 0._wp) THEN
        zfmask(ji, 1) = rn_ishlat * MIN(1._wp, MAX(zwf(ji + 1, 1), zwf(ji, 2), zwf(ji - 1, 1)))
      END IF
      IF (zfmask(ji, jpj) == 0._wp) THEN
        zfmask(ji, jpj) = rn_ishlat * MIN(1._wp, MAX(zwf(ji + 1, jpj), zwf(ji - 1, jpj), zwf(ji, jpjm1)))
      END IF
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk(zfmask, 'F', 1._wp)
    zrhoco = rau0 * rn_cio
    ecc2 = rn_ecc * rn_ecc
    z1_ecc2 = 1._wp / ecc2
    zdtevp = rdt_ice / REAL(nn_nevp)
    z1_dtevp = 1._wp / zdtevp
    IF (.NOT. ln_aEVP) THEN
      zalph1 = (2._wp * rn_relast * rdt_ice) * z1_dtevp
      zalph2 = zalph1 * z1_ecc2
      z1_alph1 = 1._wp / (zalph1 + 1._wp)
      z1_alph2 = 1._wp / (zalph2 + 1._wp)
    END IF
    !$ACC KERNELS
    zs1(:, :) = pstress1_i(:, :)
    zs2(:, :) = pstress2_i(:, :)
    zs12(:, :) = pstress12_i(:, :)
    !$ACC END KERNELS
    CALL ice_strength
    !$ACC KERNELS
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        z1_e1t0(ji, jj) = 1._wp / (e1t(ji + 1, jj) + e1t(ji, jj))
        z1_e2t0(ji, jj) = 1._wp / (e2t(ji, jj + 1) + e2t(ji, jj))
      END DO
    END DO
    !$ACC END KERNELS
    IF (ln_ice_embd) THEN
      !$ACC KERNELS
      zintn = REAL(nn_fsbc - 1) / REAL(nn_fsbc) * 0.5_wp
      zintb = REAL(nn_fsbc + 1) / REAL(nn_fsbc) * 0.5_wp
      zpice(:, :) = ssh_m(:, :) + (zintn * snwice_mass(:, :) + zintb * snwice_mass_b(:, :)) * r1_rau0
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      zpice(:, :) = ssh_m(:, :)
      !$ACC END KERNELS
    END IF
    !$ACC KERNELS
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zaU(ji, jj) = 0.5_wp * (at_i(ji, jj) * e1e2t(ji, jj) + at_i(ji + 1, jj) * e1e2t(ji + 1, jj)) * r1_e1e2u(ji, jj) * umask(ji, jj, 1)
        zaV(ji, jj) = 0.5_wp * (at_i(ji, jj) * e1e2t(ji, jj) + at_i(ji, jj + 1) * e1e2t(ji, jj + 1)) * r1_e1e2v(ji, jj) * vmask(ji, jj, 1)
        zm1 = (rhos * vt_s(ji, jj) + rhoi * vt_i(ji, jj))
        zm2 = (rhos * vt_s(ji + 1, jj) + rhoi * vt_i(ji + 1, jj))
        zm3 = (rhos * vt_s(ji, jj + 1) + rhoi * vt_i(ji, jj + 1))
        zmassU = 0.5_wp * (zm1 * e1e2t(ji, jj) + zm2 * e1e2t(ji + 1, jj)) * r1_e1e2u(ji, jj) * umask(ji, jj, 1)
        zmassV = 0.5_wp * (zm1 * e1e2t(ji, jj) + zm3 * e1e2t(ji, jj + 1)) * r1_e1e2v(ji, jj) * vmask(ji, jj, 1)
        v_oceU(ji, jj) = 0.5_wp * ((v_oce(ji, jj) + v_oce(ji, jj - 1)) * e1t(ji + 1, jj) + (v_oce(ji + 1, jj) + v_oce(ji + 1, jj - 1)) * e1t(ji, jj)) * z1_e1t0(ji, jj) * umask(ji, jj, 1)
        u_oceV(ji, jj) = 0.5_wp * ((u_oce(ji, jj) + u_oce(ji - 1, jj)) * e2t(ji, jj + 1) + (u_oce(ji, jj + 1) + u_oce(ji - 1, jj + 1)) * e2t(ji, jj)) * z1_e2t0(ji, jj) * vmask(ji, jj, 1)
        zmf(ji, jj) = zm1 * ff_t(ji, jj)
        zdt_m(ji, jj) = zdtevp / MAX(zm1, zmmin)
        zmU_t(ji, jj) = zmassU * z1_dtevp
        zmV_t(ji, jj) = zmassV * z1_dtevp
        zTauU_ia(ji, jj) = zaU(ji, jj) * utau_ice(ji, jj)
        zTauV_ia(ji, jj) = zaV(ji, jj) * vtau_ice(ji, jj)
        zspgU(ji, jj) = - zmassU * grav * (zpice(ji + 1, jj) - zpice(ji, jj)) * r1_e1u(ji, jj)
        zspgV(ji, jj) = - zmassV * grav * (zpice(ji, jj + 1) - zpice(ji, jj)) * r1_e2v(ji, jj)
        zmaskU(ji, jj) = 1._wp - MAX(0._wp, SIGN(1._wp, - zmassU))
        zmaskV(ji, jj) = 1._wp - MAX(0._wp, SIGN(1._wp, - zmassV))
        zswitchU(ji, jj) = MAX(0._wp, SIGN(1._wp, zmassU - zmmin))
        zswitchV(ji, jj) = MAX(0._wp, SIGN(1._wp, zmassV - zmmin))
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk_multi(zmf, 'T', 1., zdt_m, 'T', 1.)
    DO jter = 1, nn_nevp
      IF (ln_ctl) THEN
        !$ACC KERNELS
        DO jj = 1, jpjm1
          zu_ice(:, jj) = u_ice(:, jj)
          zv_ice(:, jj) = v_ice(:, jj)
        END DO
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zds(ji, jj) = ((u_ice(ji, jj + 1) * r1_e1u(ji, jj + 1) - u_ice(ji, jj) * r1_e1u(ji, jj)) * e1f(ji, jj) * e1f(ji, jj) + (v_ice(ji + 1, jj) * r1_e2v(ji + 1, jj) - v_ice(ji, jj) * r1_e2v(ji, jj)) * e2f(ji, jj) * e2f(ji, jj)) * r1_e1e2f(ji, jj) * zfmask(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk(zds, 'F', 1.)
      !$ACC KERNELS
      DO jj = 2, jpj
        DO ji = 2, jpi
          zds2 = (zds(ji, jj) * zds(ji, jj) * e1e2f(ji, jj) + zds(ji - 1, jj) * zds(ji - 1, jj) * e1e2f(ji - 1, jj) + zds(ji, jj - 1) * zds(ji, jj - 1) * e1e2f(ji, jj - 1) + zds(ji - 1, jj - 1) * zds(ji - 1, jj - 1) * e1e2f(ji - 1, jj - 1)) * 0.25_wp * r1_e1e2t(ji, jj)
          zdiv = (e2u(ji, jj) * u_ice(ji, jj) - e2u(ji - 1, jj) * u_ice(ji - 1, jj) + e1v(ji, jj) * v_ice(ji, jj) - e1v(ji, jj - 1) * v_ice(ji, jj - 1)) * r1_e1e2t(ji, jj)
          zdiv2 = zdiv * zdiv
          zdt = ((u_ice(ji, jj) * r1_e2u(ji, jj) - u_ice(ji - 1, jj) * r1_e2u(ji - 1, jj)) * e2t(ji, jj) * e2t(ji, jj) - (v_ice(ji, jj) * r1_e1v(ji, jj) - v_ice(ji, jj - 1) * r1_e1v(ji, jj - 1)) * e1t(ji, jj) * e1t(ji, jj)) * r1_e1e2t(ji, jj)
          zdt2 = zdt * zdt
          zdelta = SQRT(zdiv2 + (zdt2 + zds2) * z1_ecc2)
          zp_delt(ji, jj) = strength(ji, jj) / (zdelta + rn_creepl)
          IF (ln_aEVP) THEN
            zalph1 = MAX(50._wp, rpi * SQRT(0.5_wp * zp_delt(ji, jj) * r1_e1e2t(ji, jj) * zdt_m(ji, jj)))
            z1_alph1 = 1._wp / (zalph1 + 1._wp)
            zalph2 = zalph1
            z1_alph2 = z1_alph1
          END IF
          zs1(ji, jj) = (zs1(ji, jj) * zalph1 + zp_delt(ji, jj) * (zdiv - zdelta)) * z1_alph1
          zs2(ji, jj) = (zs2(ji, jj) * zalph2 + zp_delt(ji, jj) * (zdt * z1_ecc2)) * z1_alph2
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk(zp_delt, 'T', 1.)
      !$ACC KERNELS
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          IF (ln_aEVP) THEN
            zalph2 = MAX(50._wp, rpi * SQRT(0.5_wp * zp_delt(ji, jj) * r1_e1e2t(ji, jj) * zdt_m(ji, jj)))
            z1_alph2 = 1._wp / (zalph2 + 1._wp)
            zbeta(ji, jj) = zalph2
          END IF
          zp_delf = 0.25_wp * (zp_delt(ji, jj) + zp_delt(ji + 1, jj) + zp_delt(ji, jj + 1) + zp_delt(ji + 1, jj + 1))
          zs12(ji, jj) = (zs12(ji, jj) * zalph2 + zp_delf * (zds(ji, jj) * z1_ecc2) * 0.5_wp) * z1_alph2
        END DO
      END DO
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zfU(ji, jj) = 0.5_wp * ((zs1(ji + 1, jj) - zs1(ji, jj)) * e2u(ji, jj) + (zs2(ji + 1, jj) * e2t(ji + 1, jj) * e2t(ji + 1, jj) - zs2(ji, jj) * e2t(ji, jj) * e2t(ji, jj)) * r1_e2u(ji, jj) + (zs12(ji, jj) * e1f(ji, jj) * e1f(ji, jj) - zs12(ji, jj - 1) * e1f(ji, jj - 1) * e1f(ji, jj - 1)) * 2._wp * r1_e1u(ji, jj)) * r1_e1e2u(ji, jj)
          zfV(ji, jj) = 0.5_wp * ((zs1(ji, jj + 1) - zs1(ji, jj)) * e1v(ji, jj) - (zs2(ji, jj + 1) * e1t(ji, jj + 1) * e1t(ji, jj + 1) - zs2(ji, jj) * e1t(ji, jj) * e1t(ji, jj)) * r1_e1v(ji, jj) + (zs12(ji, jj) * e2f(ji, jj) * e2f(ji, jj) - zs12(ji - 1, jj) * e2f(ji - 1, jj) * e2f(ji - 1, jj)) * 2._wp * r1_e2v(ji, jj)) * r1_e1e2v(ji, jj)
          u_iceV(ji, jj) = 0.5_wp * ((u_ice(ji, jj) + u_ice(ji - 1, jj)) * e2t(ji, jj + 1) + (u_ice(ji, jj + 1) + u_ice(ji - 1, jj + 1)) * e2t(ji, jj)) * z1_e2t0(ji, jj) * vmask(ji, jj, 1)
          v_iceU(ji, jj) = 0.5_wp * ((v_ice(ji, jj) + v_ice(ji, jj - 1)) * e1t(ji + 1, jj) + (v_ice(ji + 1, jj) + v_ice(ji + 1, jj - 1)) * e1t(ji, jj)) * z1_e1t0(ji, jj) * umask(ji, jj, 1)
        END DO
      END DO
      !$ACC END KERNELS
      IF (MOD(jter, 2) == 0) THEN
        !$ACC KERNELS
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zTauO = zaV(ji, jj) * zrhoco * SQRT((v_ice(ji, jj) - v_oce(ji, jj)) * (v_ice(ji, jj) - v_oce(ji, jj)) + (u_iceV(ji, jj) - u_oceV(ji, jj)) * (u_iceV(ji, jj) - u_oceV(ji, jj)))
            ztauy_oi(ji, jj) = zTauO * (v_oce(ji, jj) - v_ice(ji, jj))
            zvel = MAX(zepsi, SQRT(v_ice(ji, jj) * v_ice(ji, jj) + u_iceV(ji, jj) * u_iceV(ji, jj)))
            zTauB = - tau_icebfr(ji, jj) / zvel
            zCory(ji, jj) = - 0.25_wp * r1_e2v(ji, jj) * (zmf(ji, jj) * (e2u(ji, jj) * u_ice(ji, jj) + e2u(ji - 1, jj) * u_ice(ji - 1, jj)) + zmf(ji, jj + 1) * (e2u(ji, jj + 1) * u_ice(ji, jj + 1) + e2u(ji - 1, jj + 1) * u_ice(ji - 1, jj + 1)))
            zTauE = zfV(ji, jj) + zTauV_ia(ji, jj) + zCory(ji, jj) + zspgV(ji, jj) + ztauy_oi(ji, jj)
            rswitch = 1._wp - MIN(1._wp, ABS(SIGN(1._wp, ztauE - tau_icebfr(ji, jj)) - SIGN(1._wp, zTauE)))
            IF (ln_aEVP) THEN
              v_ice(ji, jj) = ((rswitch * (zmV_t(ji, jj) * (zbeta(ji, jj) * v_ice(ji, jj) + v_ice_b(ji, jj)) + zTauE + zTauO * v_ice(ji, jj)) / MAX(zepsi, zmV_t(ji, jj) * (zbeta(ji, jj) + 1._wp) + zTauO - zTauB) + (1._wp - rswitch) * v_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchV(ji, jj) + v_oce(ji, jj) * (1._wp - zswitchV(ji, jj))) * zmaskV(ji, jj)
            ELSE
              v_ice(ji, jj) = ((rswitch * (zmV_t(ji, jj) * v_ice(ji, jj) + zTauE + zTauO * v_ice(ji, jj)) / MAX(zepsi, zmV_t(ji, jj) + zTauO - zTauB) + (1._wp - rswitch) * v_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchV(ji, jj) + v_oce(ji, jj) * (1._wp - zswitchV(ji, jj))) * zmaskV(ji, jj)
            END IF
          END DO
        END DO
        !$ACC END KERNELS
        CALL lbc_lnk(v_ice, 'V', - 1.)
        IF (ln_bdy) CALL bdy_ice_dyn('V')
        !$ACC KERNELS
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zTauO = zaU(ji, jj) * zrhoco * SQRT((u_ice(ji, jj) - u_oce(ji, jj)) * (u_ice(ji, jj) - u_oce(ji, jj)) + (v_iceU(ji, jj) - v_oceU(ji, jj)) * (v_iceU(ji, jj) - v_oceU(ji, jj)))
            ztaux_oi(ji, jj) = zTauO * (u_oce(ji, jj) - u_ice(ji, jj))
            zvel = MAX(zepsi, SQRT(v_iceU(ji, jj) * v_iceU(ji, jj) + u_ice(ji, jj) * u_ice(ji, jj)))
            zTauB = - tau_icebfr(ji, jj) / zvel
            zCorx(ji, jj) = 0.25_wp * r1_e1u(ji, jj) * (zmf(ji, jj) * (e1v(ji, jj) * v_ice(ji, jj) + e1v(ji, jj - 1) * v_ice(ji, jj - 1)) + zmf(ji + 1, jj) * (e1v(ji + 1, jj) * v_ice(ji + 1, jj) + e1v(ji + 1, jj - 1) * v_ice(ji + 1, jj - 1)))
            zTauE = zfU(ji, jj) + zTauU_ia(ji, jj) + zCorx(ji, jj) + zspgU(ji, jj) + ztaux_oi(ji, jj)
            rswitch = 1._wp - MIN(1._wp, ABS(SIGN(1._wp, ztauE - tau_icebfr(ji, jj)) - SIGN(1._wp, zTauE)))
            IF (ln_aEVP) THEN
              u_ice(ji, jj) = ((rswitch * (zmU_t(ji, jj) * (zbeta(ji, jj) * u_ice(ji, jj) + u_ice_b(ji, jj)) + zTauE + zTauO * u_ice(ji, jj)) / MAX(zepsi, zmU_t(ji, jj) * (zbeta(ji, jj) + 1._wp) + zTauO - zTauB) + (1._wp - rswitch) * u_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchU(ji, jj) + u_oce(ji, jj) * (1._wp - zswitchU(ji, jj))) * zmaskU(ji, jj)
            ELSE
              u_ice(ji, jj) = ((rswitch * (zmU_t(ji, jj) * u_ice(ji, jj) + zTauE + zTauO * u_ice(ji, jj)) / MAX(zepsi, zmU_t(ji, jj) + zTauO - zTauB) + (1._wp - rswitch) * u_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchU(ji, jj) + u_oce(ji, jj) * (1._wp - zswitchU(ji, jj))) * zmaskU(ji, jj)
            END IF
          END DO
        END DO
        !$ACC END KERNELS
        CALL lbc_lnk(u_ice, 'U', - 1.)
        IF (ln_bdy) CALL bdy_ice_dyn('U')
      ELSE
        !$ACC KERNELS
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zTauO = zaU(ji, jj) * zrhoco * SQRT((u_ice(ji, jj) - u_oce(ji, jj)) * (u_ice(ji, jj) - u_oce(ji, jj)) + (v_iceU(ji, jj) - v_oceU(ji, jj)) * (v_iceU(ji, jj) - v_oceU(ji, jj)))
            ztaux_oi(ji, jj) = zTauO * (u_oce(ji, jj) - u_ice(ji, jj))
            zvel = MAX(zepsi, SQRT(v_iceU(ji, jj) * v_iceU(ji, jj) + u_ice(ji, jj) * u_ice(ji, jj)))
            zTauB = - tau_icebfr(ji, jj) / zvel
            zCorx(ji, jj) = 0.25_wp * r1_e1u(ji, jj) * (zmf(ji, jj) * (e1v(ji, jj) * v_ice(ji, jj) + e1v(ji, jj - 1) * v_ice(ji, jj - 1)) + zmf(ji + 1, jj) * (e1v(ji + 1, jj) * v_ice(ji + 1, jj) + e1v(ji + 1, jj - 1) * v_ice(ji + 1, jj - 1)))
            zTauE = zfU(ji, jj) + zTauU_ia(ji, jj) + zCorx(ji, jj) + zspgU(ji, jj) + ztaux_oi(ji, jj)
            rswitch = 1._wp - MIN(1._wp, ABS(SIGN(1._wp, ztauE - tau_icebfr(ji, jj)) - SIGN(1._wp, zTauE)))
            IF (ln_aEVP) THEN
              u_ice(ji, jj) = ((rswitch * (zmU_t(ji, jj) * (zbeta(ji, jj) * u_ice(ji, jj) + u_ice_b(ji, jj)) + zTauE + zTauO * u_ice(ji, jj)) / MAX(zepsi, zmU_t(ji, jj) * (zbeta(ji, jj) + 1._wp) + zTauO - zTauB) + (1._wp - rswitch) * u_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchU(ji, jj) + u_oce(ji, jj) * (1._wp - zswitchU(ji, jj))) * zmaskU(ji, jj)
            ELSE
              u_ice(ji, jj) = ((rswitch * (zmU_t(ji, jj) * u_ice(ji, jj) + zTauE + zTauO * u_ice(ji, jj)) / MAX(zepsi, zmU_t(ji, jj) + zTauO - zTauB) + (1._wp - rswitch) * u_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchU(ji, jj) + u_oce(ji, jj) * (1._wp - zswitchU(ji, jj))) * zmaskU(ji, jj)
            END IF
          END DO
        END DO
        !$ACC END KERNELS
        CALL lbc_lnk(u_ice, 'U', - 1.)
        IF (ln_bdy) CALL bdy_ice_dyn('U')
        !$ACC KERNELS
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zTauO = zaV(ji, jj) * zrhoco * SQRT((v_ice(ji, jj) - v_oce(ji, jj)) * (v_ice(ji, jj) - v_oce(ji, jj)) + (u_iceV(ji, jj) - u_oceV(ji, jj)) * (u_iceV(ji, jj) - u_oceV(ji, jj)))
            ztauy_oi(ji, jj) = zTauO * (v_oce(ji, jj) - v_ice(ji, jj))
            zvel = MAX(zepsi, SQRT(v_ice(ji, jj) * v_ice(ji, jj) + u_iceV(ji, jj) * u_iceV(ji, jj)))
            ztauB = - tau_icebfr(ji, jj) / zvel
            zCory(ji, jj) = - 0.25_wp * r1_e2v(ji, jj) * (zmf(ji, jj) * (e2u(ji, jj) * u_ice(ji, jj) + e2u(ji - 1, jj) * u_ice(ji - 1, jj)) + zmf(ji, jj + 1) * (e2u(ji, jj + 1) * u_ice(ji, jj + 1) + e2u(ji - 1, jj + 1) * u_ice(ji - 1, jj + 1)))
            zTauE = zfV(ji, jj) + zTauV_ia(ji, jj) + zCory(ji, jj) + zspgV(ji, jj) + ztauy_oi(ji, jj)
            rswitch = 1._wp - MIN(1._wp, ABS(SIGN(1._wp, zTauE - tau_icebfr(ji, jj)) - SIGN(1._wp, zTauE)))
            IF (ln_aEVP) THEN
              v_ice(ji, jj) = ((rswitch * (zmV_t(ji, jj) * (zbeta(ji, jj) * v_ice(ji, jj) + v_ice_b(ji, jj)) + zTauE + zTauO * v_ice(ji, jj)) / MAX(zepsi, zmV_t(ji, jj) * (zbeta(ji, jj) + 1._wp) + zTauO - zTauB) + (1._wp - rswitch) * v_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchV(ji, jj) + v_oce(ji, jj) * (1._wp - zswitchV(ji, jj))) * zmaskV(ji, jj)
            ELSE
              v_ice(ji, jj) = ((rswitch * (zmV_t(ji, jj) * v_ice(ji, jj) + zTauE + zTauO * v_ice(ji, jj)) / MAX(zepsi, zmV_t(ji, jj) + zTauO - zTauB) + (1._wp - rswitch) * v_ice(ji, jj) * MAX(0._wp, 1._wp - zdtevp * rn_lfrelax)) * zswitchV(ji, jj) + v_oce(ji, jj) * (1._wp - zswitchV(ji, jj))) * zmaskV(ji, jj)
            END IF
          END DO
        END DO
        !$ACC END KERNELS
        CALL lbc_lnk(v_ice, 'V', - 1.)
        IF (ln_bdy) CALL bdy_ice_dyn('V')
      END IF
      IF (ln_ctl) THEN
        !$ACC KERNELS
        DO jj = 2, jpjm1
          zresr(:, jj) = MAX(ABS(u_ice(:, jj) - zu_ice(:, jj)), ABS(v_ice(:, jj) - zv_ice(:, jj)))
        END DO
        !$ACC END KERNELS
        zresm = MAXVAL(zresr(1 : jpi, 2 : jpjm1))
        IF (lk_mpp) CALL mpp_max(zresm)
      END IF
    END DO
    !$ACC KERNELS
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        zds(ji, jj) = ((u_ice(ji, jj + 1) * r1_e1u(ji, jj + 1) - u_ice(ji, jj) * r1_e1u(ji, jj)) * e1f(ji, jj) * e1f(ji, jj) + (v_ice(ji + 1, jj) * r1_e2v(ji + 1, jj) - v_ice(ji, jj) * r1_e2v(ji, jj)) * e2f(ji, jj) * e2f(ji, jj)) * r1_e1e2f(ji, jj) * zfmask(ji, jj)
      END DO
    END DO
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zdt = ((u_ice(ji, jj) * r1_e2u(ji, jj) - u_ice(ji - 1, jj) * r1_e2u(ji - 1, jj)) * e2t(ji, jj) * e2t(ji, jj) - (v_ice(ji, jj) * r1_e1v(ji, jj) - v_ice(ji, jj - 1) * r1_e1v(ji, jj - 1)) * e1t(ji, jj) * e1t(ji, jj)) * r1_e1e2t(ji, jj)
        zdt2 = zdt * zdt
        zds2 = (zds(ji, jj) * zds(ji, jj) * e1e2f(ji, jj) + zds(ji - 1, jj) * zds(ji - 1, jj) * e1e2f(ji - 1, jj) + zds(ji, jj - 1) * zds(ji, jj - 1) * e1e2f(ji, jj - 1) + zds(ji - 1, jj - 1) * zds(ji - 1, jj - 1) * e1e2f(ji - 1, jj - 1)) * 0.25_wp * r1_e1e2t(ji, jj)
        pshear_i(ji, jj) = SQRT(zdt2 + zds2)
        pdivu_i(ji, jj) = (e2u(ji, jj) * u_ice(ji, jj) - e2u(ji - 1, jj) * u_ice(ji - 1, jj) + e1v(ji, jj) * v_ice(ji, jj) - e1v(ji, jj - 1) * v_ice(ji, jj - 1)) * r1_e1e2t(ji, jj)
        zdelta = SQRT(pdivu_i(ji, jj) * pdivu_i(ji, jj) + (zdt2 + zds2) * z1_ecc2)
        rswitch = 1._wp - MAX(0._wp, SIGN(1._wp, - zdelta))
        pdelta_i(ji, jj) = zdelta + rn_creepl * rswitch
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk_multi(pshear_i, 'T', 1., pdivu_i, 'T', 1., pdelta_i, 'T', 1.)
    CALL lbc_lnk_multi(zs1, 'T', 1., zs2, 'T', 1., zs12, 'F', 1.)
    !$ACC KERNELS
    pstress1_i(:, :) = zs1(:, :)
    pstress2_i(:, :) = zs2(:, :)
    pstress12_i(:, :) = zs12(:, :)
    DO jj = 1, jpj
      DO ji = 1, jpi
        zswi(ji, jj) = MAX(0._wp, SIGN(1._wp, at_i(ji, jj) - epsi06))
      END DO
    END DO
    !$ACC END KERNELS
    IF (iom_use('icediv')) CALL iom_put("icediv", pdivu_i(:, :) * zswi(:, :))
    IF (iom_use('iceshe')) CALL iom_put("iceshe", pshear_i(:, :) * zswi(:, :))
    IF (iom_use('icestr')) CALL iom_put("icestr", strength(:, :) * zswi(:, :))
    IF (iom_use('isig1') .OR. iom_use('isig2') .OR. iom_use('isig3')) THEN
      ALLOCATE(zsig1(jpi, jpj), zsig2(jpi, jpj), zsig3(jpi, jpj))
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zdum1 = (zswi(ji - 1, jj) * pstress12_i(ji - 1, jj) + zswi(ji, jj - 1) * pstress12_i(ji, jj - 1) + zswi(ji, jj) * pstress12_i(ji, jj) + zswi(ji - 1, jj - 1) * pstress12_i(ji - 1, jj - 1)) / MAX(1._wp, zswi(ji - 1, jj) + zswi(ji, jj - 1) + zswi(ji, jj) + zswi(ji - 1, jj - 1))
          zshear = SQRT(pstress2_i(ji, jj) * pstress2_i(ji, jj) + 4._wp * zdum1 * zdum1)
          zdum2 = zswi(ji, jj) / MAX(1._wp, strength(ji, jj))
          zsig1(ji, jj) = 0.5_wp * zdum2 * (pstress1_i(ji, jj))
          zsig2(ji, jj) = 0.5_wp * zdum2 * (zshear)
          zsig3(ji, jj) = zdum2 ** 2 * ((pstress1_i(ji, jj) + strength(ji, jj)) ** 2 + (rn_ecc * zshear) ** 2)
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk_multi(zsig1, 'T', 1., zsig2, 'T', 1., zsig3, 'T', 1.)
      IF (iom_use('isig1')) CALL iom_put("isig1", zsig1)
      IF (iom_use('isig2')) CALL iom_put("isig2", zsig2)
      IF (iom_use('isig3')) CALL iom_put("isig3", zsig3)
      DEALLOCATE(zsig1, zsig2, zsig3)
    END IF
    IF (iom_use('normstr') .OR. iom_use('sheastr') .OR. iom_use('dssh_dx') .OR. iom_use('dssh_dy') .OR. iom_use('corstrx') .OR. iom_use('corstry') .OR. iom_use('intstrx') .OR. iom_use('intstry') .OR. iom_use('utau_oi') .OR. iom_use('vtau_oi') .OR. iom_use('xmtrpice') .OR. iom_use('ymtrpice') .OR. iom_use('xmtrpsnw') .OR. iom_use('ymtrpsnw') .OR. iom_use('xatrp') .OR. iom_use('yatrp')) THEN
      ALLOCATE(zdiag_sig1(jpi, jpj), zdiag_sig2(jpi, jpj), zdiag_dssh_dx(jpi, jpj), zdiag_dssh_dy(jpi, jpj), zdiag_corstrx(jpi, jpj), zdiag_corstry(jpi, jpj), zdiag_intstrx(jpi, jpj), zdiag_intstry(jpi, jpj), zdiag_utau_oi(jpi, jpj), zdiag_vtau_oi(jpi, jpj), zdiag_xmtrp_ice(jpi, jpj), zdiag_ymtrp_ice(jpi, jpj), zdiag_xmtrp_snw(jpi, jpj), zdiag_ymtrp_snw(jpi, jpj), zdiag_xatrp(jpi, jpj), zdiag_yatrp(jpi, jpj))
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          rswitch = MAX(0._wp, SIGN(1._wp, at_i(ji, jj) - epsi06))
          zdiag_sig1(ji, jj) = (zs1(ji, jj) + zs2(ji, jj)) * rswitch
          zdiag_sig2(ji, jj) = SQRT((zs1(ji, jj) - zs2(ji, jj)) ** 2 + 4 * zs12(ji, jj) ** 2) * rswitch
          zdiag_dssh_dx(ji, jj) = zspgU(ji, jj) * rswitch
          zdiag_dssh_dy(ji, jj) = zspgV(ji, jj) * rswitch
          zdiag_corstrx(ji, jj) = zCorx(ji, jj) * rswitch
          zdiag_corstry(ji, jj) = zCory(ji, jj) * rswitch
          zdiag_intstrx(ji, jj) = zfU(ji, jj) * rswitch
          zdiag_intstry(ji, jj) = zfV(ji, jj) * rswitch
          zdiag_utau_oi(ji, jj) = ztaux_oi(ji, jj) * rswitch
          zdiag_vtau_oi(ji, jj) = ztauy_oi(ji, jj) * rswitch
          zfac_x = 0.5 * u_ice(ji, jj) * e2u(ji, jj) * rswitch
          zfac_y = 0.5 * v_ice(ji, jj) * e1v(ji, jj) * rswitch
          zdiag_xmtrp_ice(ji, jj) = rhoi * zfac_x * (vt_i(ji + 1, jj) + vt_i(ji, jj))
          zdiag_ymtrp_ice(ji, jj) = rhoi * zfac_y * (vt_i(ji, jj + 1) + vt_i(ji, jj))
          zdiag_xmtrp_snw(ji, jj) = rhos * zfac_x * (vt_s(ji + 1, jj) + vt_s(ji, jj))
          zdiag_ymtrp_snw(ji, jj) = rhos * zfac_y * (vt_s(ji, jj + 1) + vt_s(ji, jj))
          zdiag_xatrp(ji, jj) = zfac_x * (at_i(ji + 1, jj) + at_i(ji, jj))
          zdiag_yatrp(ji, jj) = zfac_y * (at_i(ji, jj + 1) + at_i(ji, jj))
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk_multi(zdiag_sig1, 'T', 1., zdiag_sig2, 'T', 1., zdiag_dssh_dx, 'U', - 1., zdiag_dssh_dy, 'V', - 1., zdiag_corstrx, 'U', - 1., zdiag_corstry, 'V', - 1., zdiag_intstrx, 'U', - 1., zdiag_intstry, 'V', - 1.)
      CALL lbc_lnk_multi(zdiag_utau_oi, 'U', - 1., zdiag_vtau_oi, 'V', - 1., zdiag_xmtrp_ice, 'U', - 1., zdiag_xmtrp_snw, 'U', - 1., zdiag_xatrp, 'U', - 1., zdiag_ymtrp_ice, 'V', - 1., zdiag_ymtrp_snw, 'V', - 1., zdiag_yatrp, 'V', - 1.)
      IF (iom_use('normstr')) CALL iom_put('normstr', zdiag_sig1(:, :))
      IF (iom_use('sheastr')) CALL iom_put('sheastr', zdiag_sig2(:, :))
      IF (iom_use('dssh_dx')) CALL iom_put('dssh_dx', zdiag_dssh_dx(:, :))
      IF (iom_use('dssh_dy')) CALL iom_put('dssh_dy', zdiag_dssh_dy(:, :))
      IF (iom_use('corstrx')) CALL iom_put('corstrx', zdiag_corstrx(:, :))
      IF (iom_use('corstry')) CALL iom_put('corstry', zdiag_corstry(:, :))
      IF (iom_use('intstrx')) CALL iom_put('intstrx', zdiag_intstrx(:, :))
      IF (iom_use('intstry')) CALL iom_put('intstry', zdiag_intstry(:, :))
      IF (iom_use('utau_oi')) CALL iom_put('utau_oi', zdiag_utau_oi(:, :))
      IF (iom_use('vtau_oi')) CALL iom_put('vtau_oi', zdiag_vtau_oi(:, :))
      IF (iom_use('xmtrpice')) CALL iom_put('xmtrpice', zdiag_xmtrp_ice(:, :))
      IF (iom_use('ymtrpice')) CALL iom_put('ymtrpice', zdiag_ymtrp_ice(:, :))
      IF (iom_use('xmtrpsnw')) CALL iom_put('xmtrpsnw', zdiag_xmtrp_snw(:, :))
      IF (iom_use('ymtrpsnw')) CALL iom_put('ymtrpsnw', zdiag_ymtrp_snw(:, :))
      IF (iom_use('xatrp')) CALL iom_put('xatrp', zdiag_xatrp(:, :))
      IF (iom_use('yatrp')) CALL iom_put('yatrp', zdiag_yatrp(:, :))
      DEALLOCATE(zdiag_sig1, zdiag_sig2, zdiag_dssh_dx, zdiag_dssh_dy, zdiag_corstrx, zdiag_corstry, zdiag_intstrx, zdiag_intstry, zdiag_utau_oi, zdiag_vtau_oi, zdiag_xmtrp_ice, zdiag_ymtrp_ice, zdiag_xmtrp_snw, zdiag_ymtrp_snw, zdiag_xatrp, zdiag_yatrp)
    END IF
  END SUBROUTINE ice_dyn_rhg_evp
  SUBROUTINE rhg_evp_rst(cdrw, kt)
    CHARACTER(LEN = *), INTENT(IN) :: cdrw
    INTEGER, OPTIONAL, INTENT(IN) :: kt
    INTEGER :: iter
    INTEGER :: id1, id2, id3
    IF (TRIM(cdrw) == 'READ') THEN
      IF (ln_rstart) THEN
        id1 = iom_varid(numrir, 'stress1_i', ldstop = .FALSE.)
        id2 = iom_varid(numrir, 'stress2_i', ldstop = .FALSE.)
        id3 = iom_varid(numrir, 'stress12_i', ldstop = .FALSE.)
        IF (MIN(id1, id2, id3) > 0) THEN
          CALL iom_get(numrir, jpdom_autoglo, 'stress1_i', stress1_i)
          CALL iom_get(numrir, jpdom_autoglo, 'stress2_i', stress2_i)
          CALL iom_get(numrir, jpdom_autoglo, 'stress12_i', stress12_i)
        ELSE
          IF (lwp) WRITE(numout, FMT = *)
          IF (lwp) WRITE(numout, FMT = *) '   ==>>>   previous run without rheology, set stresses to 0'
          !$ACC KERNELS
          stress1_i(:, :) = 0._wp
          stress2_i(:, :) = 0._wp
          stress12_i(:, :) = 0._wp
          !$ACC END KERNELS
        END IF
      ELSE
        IF (lwp) WRITE(numout, FMT = *)
        IF (lwp) WRITE(numout, FMT = *) '   ==>>>   start from rest: set stresses to 0'
        !$ACC KERNELS
        stress1_i(:, :) = 0._wp
        stress2_i(:, :) = 0._wp
        stress12_i(:, :) = 0._wp
        !$ACC END KERNELS
      END IF
    ELSE IF (TRIM(cdrw) == 'WRITE') THEN
      IF (lwp) WRITE(numout, FMT = *) '---- rhg-rst ----'
      iter = kt + nn_fsbc - 1
      CALL iom_rstput(iter, nitrst, numriw, 'stress1_i', stress1_i)
      CALL iom_rstput(iter, nitrst, numriw, 'stress2_i', stress2_i)
      CALL iom_rstput(iter, nitrst, numriw, 'stress12_i', stress12_i)
    END IF
  END SUBROUTINE rhg_evp_rst
END MODULE icedyn_rhg_evp