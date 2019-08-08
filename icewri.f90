MODULE icewri
  USE dianam
  USE phycst
  USE dom_oce
  USE sbc_oce
  USE sbc_ice
  USE ice
  USE icevar
  USE ioipsl
  USE in_out_manager
  USE iom
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_wri
  PUBLIC :: ice_wri_state
  CONTAINS
  SUBROUTINE ice_wri(kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT(IN) :: kt
    INTEGER :: ji, jj, jk, jl
    REAL(KIND = wp) :: z2da, z2db, zrho1, zrho2
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: z2d
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zmsk00, zmsk05, zmsk15, zmsksn
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpl) :: zmsk00l, zmsksnl
    REAL(KIND = wp) :: zdiag_area_nh, zdiag_extt_nh, zdiag_volu_nh
    REAL(KIND = wp) :: zdiag_area_sh, zdiag_extt_sh, zdiag_volu_sh
    TYPE(ProfileData), SAVE :: psy_profile0
    TYPE(ProfileData), SAVE :: psy_profile1
    TYPE(ProfileData), SAVE :: psy_profile2
    TYPE(ProfileData), SAVE :: psy_profile3
    TYPE(ProfileData), SAVE :: psy_profile4
    TYPE(ProfileData), SAVE :: psy_profile5
    CALL ProfileStart('ice_wri', 'r0', psy_profile0)
    IF (ln_timing) CALL timing_start('icewri')
    CALL ice_var_bv
    CALL ProfileEnd(psy_profile0)
    !$ACC KERNELS
    DO jj = 1, jpj
      DO ji = 1, jpi
        zmsk00(ji, jj) = MAX(0._wp, SIGN(1._wp, at_i(ji, jj) - epsi06))
        zmsk05(ji, jj) = MAX(0._wp, SIGN(1._wp, at_i(ji, jj) - 0.05_wp))
        zmsk15(ji, jj) = MAX(0._wp, SIGN(1._wp, at_i(ji, jj) - 0.15_wp))
        zmsksn(ji, jj) = MAX(0._wp, SIGN(1._wp, vt_s(ji, jj) - epsi06))
      END DO
    END DO
    DO jl = 1, jpl
      DO jj = 1, jpj
        DO ji = 1, jpi
          zmsk00l(ji, jj, jl) = MAX(0._wp, SIGN(1._wp, a_i(ji, jj, jl) - epsi06))
          zmsksnl(ji, jj, jl) = MAX(0._wp, SIGN(1._wp, v_s(ji, jj, jl) - epsi06))
        END DO
      END DO
    END DO
    zrho1 = (rau0 - rhoi) * r1_rau0
    zrho2 = rhos * r1_rau0
    !$ACC END KERNELS
    CALL ProfileStart('ice_wri', 'r1', psy_profile1)
    IF (iom_use('icemask')) CALL iom_put("icemask", zmsk00)
    IF (iom_use('icemask05')) CALL iom_put("icemask05", zmsk05)
    IF (iom_use('icemask15')) CALL iom_put("icemask15", zmsk15)
    IF (iom_use('icemass')) CALL iom_put("icemass", rhoi * vt_i * zmsk00)
    IF (iom_use('snwmass')) CALL iom_put("snwmass", rhos * vt_s * zmsksn)
    IF (iom_use('icepres')) CALL iom_put("icepres", zmsk00)
    IF (iom_use('iceconc')) CALL iom_put("iceconc", at_i * zmsk00)
    IF (iom_use('icevolu')) CALL iom_put("icevolu", vt_i * zmsk00)
    IF (iom_use('icethic')) CALL iom_put("icethic", hm_i * zmsk00)
    IF (iom_use('snwthic')) CALL iom_put("snwthic", hm_s * zmsk00)
    IF (iom_use('icebrv')) CALL iom_put("icebrv", bvm_i * zmsk00 * 100.)
    IF (iom_use('iceage')) CALL iom_put("iceage", om_i * zmsk00 / rday)
    IF (iom_use('icehnew')) CALL iom_put("icehnew", ht_i_new)
    IF (iom_use('snwvolu')) CALL iom_put("snwvolu", vt_s * zmsksn)
    CALL ProfileEnd(psy_profile1)
    IF (iom_use('icefrb')) THEN
      !$ACC KERNELS
      z2d(:, :) = (zrho1 * hm_i(:, :) - zrho2 * hm_s(:, :))
      !$ACC END KERNELS
      CALL ProfileStart('ice_wri', 'r2', psy_profile2)
      WHERE (z2d < 0._wp) z2d = 0._wp
      CALL iom_put("icefrb", z2d * zmsk00)
      CALL ProfileEnd(psy_profile2)
    END IF
    CALL ProfileStart('ice_wri', 'r3', psy_profile3)
    IF (iom_use('iceapnd')) CALL iom_put("iceapnd", at_ip * zmsk00)
    IF (iom_use('icevpnd')) CALL iom_put("icevpnd", vt_ip * zmsk00)
    IF (iom_use('icesalt')) CALL iom_put("icesalt", sm_i * zmsk00)
    IF (iom_use('icesalm')) CALL iom_put("icesalm", SUM(sv_i, DIM = 3) * rhoi * 1.0E-3 * zmsk00)
    IF (iom_use('icetemp')) CALL iom_put("icetemp", (tm_i - rt0) * zmsk00)
    IF (iom_use('snwtemp')) CALL iom_put("snwtemp", (tm_s - rt0) * zmsksn)
    IF (iom_use('icettop')) CALL iom_put("icettop", (tm_su - rt0) * zmsk00)
    IF (iom_use('icetbot')) CALL iom_put("icetbot", (t_bo - rt0) * zmsk00)
    IF (iom_use('icetsni')) CALL iom_put("icetsni", (tm_si - rt0) * zmsk00)
    IF (iom_use('icehc')) CALL iom_put("icehc", - et_i * zmsk00)
    IF (iom_use('snwhc')) CALL iom_put("snwhc", - et_s * zmsksn)
    IF (iom_use('uice')) CALL iom_put("uice", u_ice)
    IF (iom_use('vice')) CALL iom_put("vice", v_ice)
    IF (iom_use('utau_ai')) CALL iom_put("utau_ai", utau_ice * zmsk00)
    IF (iom_use('vtau_ai')) CALL iom_put("vtau_ai", vtau_ice * zmsk00)
    CALL ProfileEnd(psy_profile3)
    IF (iom_use('icevel')) THEN
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          z2da = (u_ice(ji, jj) + u_ice(ji - 1, jj))
          z2db = (v_ice(ji, jj) + v_ice(ji, jj - 1))
          z2d(ji, jj) = 0.5_wp * SQRT(z2da * z2da + z2db * z2db)
        END DO
      END DO
      !$ACC END KERNELS
      CALL ProfileStart('ice_wri', 'r4', psy_profile4)
      CALL lbc_lnk(z2d, 'T', 1.)
      IF (iom_use('icevel')) CALL iom_put("icevel", z2d)
      CALL ProfileEnd(psy_profile4)
    END IF
    CALL ProfileStart('ice_wri', 'r5', psy_profile5)
    IF (iom_use('icemask_cat')) CALL iom_put("icemask_cat", zmsk00l)
    IF (iom_use('iceconc_cat')) CALL iom_put("iceconc_cat", a_i * zmsk00l)
    IF (iom_use('icethic_cat')) CALL iom_put("icethic_cat", h_i * zmsk00l)
    IF (iom_use('snwthic_cat')) CALL iom_put("snwthic_cat", h_s * zmsksnl)
    IF (iom_use('icesalt_cat')) CALL iom_put("icesalt_cat", s_i * zmsk00l)
    IF (iom_use('iceage_cat')) CALL iom_put("iceage_cat", o_i * zmsk00l / rday)
    IF (iom_use('icetemp_cat')) CALL iom_put("icetemp_cat", (SUM(t_i(:, :, :, :), dim = 3) * r1_nlay_i - rt0) * zmsk00l)
    IF (iom_use('snwtemp_cat')) CALL iom_put("snwtemp_cat", (SUM(t_s(:, :, :, :), dim = 3) * r1_nlay_s - rt0) * zmsksnl)
    IF (iom_use('icettop_cat')) CALL iom_put("icettop_cat", (t_su - rt0) * zmsk00l)
    IF (iom_use('icebrv_cat')) CALL iom_put("icebrv_cat", bv_i * 100. * zmsk00l)
    IF (iom_use('iceapnd_cat')) CALL iom_put("iceapnd_cat", a_ip * zmsk00l)
    IF (iom_use('icehpnd_cat')) CALL iom_put("icehpnd_cat", h_ip * zmsk00l)
    IF (iom_use('iceafpnd_cat')) CALL iom_put("iceafpnd_cat", a_ip_frac * zmsk00l)
    IF (iom_use('dmithd')) CALL iom_put("dmithd", - wfx_bog - wfx_bom - wfx_sum - wfx_sni - wfx_opw - wfx_lam - wfx_res)
    IF (iom_use('dmidyn')) CALL iom_put("dmidyn", - wfx_dyn + rhoi * diag_trp_vi)
    IF (iom_use('dmiopw')) CALL iom_put("dmiopw", - wfx_opw)
    IF (iom_use('dmibog')) CALL iom_put("dmibog", - wfx_bog)
    IF (iom_use('dmisni')) CALL iom_put("dmisni", - wfx_sni)
    IF (iom_use('dmisum')) CALL iom_put("dmisum", - wfx_sum)
    IF (iom_use('dmibom')) CALL iom_put("dmibom", - wfx_bom)
    IF (iom_use('dmtsub')) CALL iom_put("dmtsub", - wfx_sub)
    IF (iom_use('dmssub')) CALL iom_put("dmssub", - wfx_snw_sub)
    IF (iom_use('dmisub')) CALL iom_put("dmisub", - wfx_ice_sub)
    IF (iom_use('dmsspr')) CALL iom_put("dmsspr", - wfx_spr)
    IF (iom_use('dmsssi')) CALL iom_put("dmsssi", wfx_sni * rhos * r1_rhoi)
    IF (iom_use('dmsmel')) CALL iom_put("dmsmel", - wfx_snw_sum)
    IF (iom_use('dmsdyn')) CALL iom_put("dmsdyn", - wfx_snw_dyn + rhos * diag_trp_vs)
    IF (iom_use('NH_icearea') .OR. iom_use('NH_icevolu') .OR. iom_use('NH_iceextt')) THEN
      WHERE (ff_t > 0._wp)
        zmsk00(:, :) = 1.0E-12
      ELSEWHERE
        zmsk00(:, :) = 0.
      END WHERE
      zdiag_area_nh = glob_sum(at_i(:, :) * zmsk00(:, :) * e1e2t(:, :))
      zdiag_volu_nh = glob_sum(vt_i(:, :) * zmsk00(:, :) * e1e2t(:, :))
      WHERE (ff_t > 0._wp .AND. at_i > 0.15)
        zmsk00(:, :) = 1.0E-12
      ELSEWHERE
        zmsk00(:, :) = 0.
      END WHERE
      zdiag_extt_nh = glob_sum(zmsk00(:, :) * e1e2t(:, :))
      IF (iom_use('NH_icearea')) CALL iom_put("NH_icearea", zdiag_area_nh)
      IF (iom_use('NH_icevolu')) CALL iom_put("NH_icevolu", zdiag_volu_nh)
      IF (iom_use('NH_iceextt')) CALL iom_put("NH_iceextt", zdiag_extt_nh)
    END IF
    IF (iom_use('SH_icearea') .OR. iom_use('SH_icevolu') .OR. iom_use('SH_iceextt')) THEN
      WHERE (ff_t < 0._wp)
        zmsk00(:, :) = 1.0E-12
      ELSEWHERE
        zmsk00(:, :) = 0.
      END WHERE
      zdiag_area_sh = glob_sum(at_i(:, :) * zmsk00(:, :) * e1e2t(:, :))
      zdiag_volu_sh = glob_sum(vt_i(:, :) * zmsk00(:, :) * e1e2t(:, :))
      WHERE (ff_t < 0._wp .AND. at_i > 0.15)
        zmsk00(:, :) = 1.0E-12
      ELSEWHERE
        zmsk00(:, :) = 0.
      END WHERE
      zdiag_extt_sh = glob_sum(zmsk00(:, :) * e1e2t(:, :))
      IF (iom_use('SH_icearea')) CALL iom_put("SH_icearea", zdiag_area_sh)
      IF (iom_use('SH_icevolu')) CALL iom_put("SH_icevolu", zdiag_volu_sh)
      IF (iom_use('SH_iceextt')) CALL iom_put("SH_iceextt", zdiag_extt_sh)
    END IF
    IF (ln_timing) CALL timing_stop('icewri')
    CALL ProfileEnd(psy_profile5)
  END SUBROUTINE ice_wri
  SUBROUTINE ice_wri_state(kt, kid, kh_i)
    INTEGER, INTENT( IN ) :: kt
    INTEGER, INTENT( IN ) :: kid, kh_i
    INTEGER :: nz_i, jl
    REAL(KIND = wp), DIMENSION(jpl) :: jcat
    !$ACC KERNELS
    DO jl = 1, jpl
      jcat(jl) = REAL(jl)
    END DO
    !$ACC END KERNELS
    CALL histvert(kid, "ncatice", "Ice Categories", "", jpl, jcat, nz_i, "up")
    CALL histdef(kid, "sithic", "Ice thickness", "m", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "siconc", "Ice concentration", "%", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "sitemp", "Ice temperature", "C", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "sivelu", "i-Ice speed ", "m/s", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "sivelv", "j-Ice speed ", "m/s", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "sistru", "i-Wind stress over ice", "Pa", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "sistrv", "j-Wind stress over ice", "Pa", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "sisflx", "Solar flx over ocean", "W/m2", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "sinflx", "NonSolar flx over ocean", "W/m2", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "snwpre", "Snow precipitation", "kg/m2/s", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "sisali", "Ice salinity", "PSU", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "sivolu", "Ice volume", "m", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "sidive", "Ice divergence", "10-8s-1", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "si_amp", "Melt pond fraction", "%", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "si_vmp", "Melt pond volume", "m", jpi, jpj, kh_i, 1, 1, 1, - 99, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "sithicat", "Ice thickness", "m", jpi, jpj, kh_i, jpl, 1, jpl, nz_i, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "siconcat", "Ice concentration", "%", jpi, jpj, kh_i, jpl, 1, jpl, nz_i, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "sisalcat", "Ice salinity", "", jpi, jpj, kh_i, jpl, 1, jpl, nz_i, 32, "inst(x)", rdt, rdt)
    CALL histdef(kid, "snthicat", "Snw thickness", "m", jpi, jpj, kh_i, jpl, 1, jpl, nz_i, 32, "inst(x)", rdt, rdt)
    CALL histend(kid, snc4set)
    CALL histwrite(kid, "sithic", kt, hm_i, jpi * jpj, (/1/))
    CALL histwrite(kid, "siconc", kt, at_i, jpi * jpj, (/1/))
    CALL histwrite(kid, "sitemp", kt, tm_i - rt0, jpi * jpj, (/1/))
    CALL histwrite(kid, "sivelu", kt, u_ice, jpi * jpj, (/1/))
    CALL histwrite(kid, "sivelv", kt, v_ice, jpi * jpj, (/1/))
    CALL histwrite(kid, "sistru", kt, utau_ice, jpi * jpj, (/1/))
    CALL histwrite(kid, "sistrv", kt, vtau_ice, jpi * jpj, (/1/))
    CALL histwrite(kid, "sisflx", kt, qsr, jpi * jpj, (/1/))
    CALL histwrite(kid, "sinflx", kt, qns, jpi * jpj, (/1/))
    CALL histwrite(kid, "snwpre", kt, sprecip, jpi * jpj, (/1/))
    CALL histwrite(kid, "sisali", kt, sm_i, jpi * jpj, (/1/))
    CALL histwrite(kid, "sivolu", kt, vt_i, jpi * jpj, (/1/))
    CALL histwrite(kid, "sidive", kt, divu_i * 1.0E8, jpi * jpj, (/1/))
    CALL histwrite(kid, "si_amp", kt, at_ip, jpi * jpj, (/1/))
    CALL histwrite(kid, "si_vmp", kt, vt_ip, jpi * jpj, (/1/))
    CALL histwrite(kid, "sithicat", kt, h_i, jpi * jpj * jpl, (/1/))
    CALL histwrite(kid, "siconcat", kt, a_i, jpi * jpj * jpl, (/1/))
    CALL histwrite(kid, "sisalcat", kt, s_i, jpi * jpj * jpl, (/1/))
    CALL histwrite(kid, "snthicat", kt, h_s, jpi * jpj * jpl, (/1/))
  END SUBROUTINE ice_wri_state
END MODULE icewri