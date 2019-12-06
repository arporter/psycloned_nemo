MODULE step
  USE step_oce
  USE iom
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: stp
  CONTAINS
  SUBROUTINE stp(kstp)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT(IN) :: kstp
    INTEGER :: ji, jj, jk, jl
    INTEGER :: indic
    INTEGER :: kcall
    TYPE(ProfileData), SAVE :: psy_profile0
    TYPE(ProfileData), SAVE :: psy_profile1
    TYPE(ProfileData), SAVE :: psy_profile2
    CALL ProfileStart('stp', 'r0', psy_profile0)
    IF (ln_timing) CALL timing_start('stp')
    indic = 0
    IF (kstp == nit000) THEN
      CALL iom_init(cxios_context)
      IF (ln_crs) CALL iom_init(TRIM(cxios_context) // "_crs")
    END IF
    IF (kstp /= nit000) CALL day(kstp)
    CALL iom_setkt(kstp - nit000 + 1, cxios_context)
    IF (ln_crs) CALL iom_setkt(kstp - nit000 + 1, TRIM(cxios_context) // "_crs")
    IF (ln_tide) CALL sbc_tide(kstp)
    IF (ln_apr_dyn) CALL sbc_apr(kstp)
    IF (ln_bdy) CALL bdy_dta(kstp, time_offset = + 1)
    CALL sbc(kstp)
    IF (ln_sto_eos) CALL sto_par(kstp)
    IF (ln_sto_eos) CALL sto_pts(tsn)
    CALL eos_rab(tsb, rab_b)
    CALL eos_rab(tsn, rab_n)
    CALL bn2(tsb, rab_b, rn2b)
    CALL bn2(tsn, rab_n, rn2)
    CALL zdf_phy(kstp)
    IF (l_ldfslp) THEN
      CALL eos(tsb, rhd, gdept_0(:, :, :))
      IF (ln_zps .AND. .NOT. ln_isfcav) CALL zps_hde(kstp, jpts, tsb, gtsu, gtsv, rhd, gru, grv)
      IF (ln_zps .AND. ln_isfcav) CALL zps_hde_isf(kstp, jpts, tsb, gtsu, gtsv, gtui, gtvi, rhd, gru, grv, grui, grvi)
      IF (ln_traldf_triad) THEN
        CALL ldf_slp_triad(kstp)
      ELSE
        CALL ldf_slp(kstp, rhd, rn2b)
      END IF
    END IF
    IF (l_ldftra_time .OR. l_ldfeiv_time) CALL ldf_tra(kstp)
    IF (l_ldfdyn_time) CALL ldf_dyn(kstp)
    CALL ssh_nxt(kstp)
    IF (.NOT. ln_linssh) CALL dom_vvl_sf_nxt(kstp)
    CALL wzv(kstp)
    CALL eos(tsn, rhd, rhop, gdept_n(:, :, :))
    IF (ln_zps .AND. .NOT. ln_isfcav) CALL zps_hde(kstp, jpts, tsn, gtsu, gtsv, rhd, gru, grv)
    IF (ln_zps .AND. ln_isfcav) CALL zps_hde_isf(kstp, jpts, tsn, gtsu, gtsv, gtui, gtvi, rhd, gru, grv, grui, grvi)
    CALL ProfileEnd(psy_profile0)
    !$ACC KERNELS
    DO jk=LBOUND(ua,3),UBOUND(ua,3)
      DO jj=LBOUND(ua,2),UBOUND(ua,2)
        DO ji=LBOUND(ua,1),UBOUND(ua,1)
          ua(ji, jj, jk) = 0._wp
        ENDDO
      ENDDO
    ENDDO
    DO jk=LBOUND(va,3),UBOUND(va,3)
      DO jj=LBOUND(va,2),UBOUND(va,2)
        DO ji=LBOUND(va,1),UBOUND(va,1)
          va(ji, jj, jk)  = 0._wp
        ENDDO
      ENDDO
    ENDDO
    !$ACC END KERNELS
    CALL ProfileStart('stp', 'r1', psy_profile1)
    IF (lk_asminc .AND. ln_asmiau .AND. ln_dyninc) CALL dyn_asm_inc(kstp)
    IF (ln_bdy) CALL bdy_dyn3d_dmp(kstp)
    CALL dyn_adv(kstp)
    CALL dyn_vor(kstp)
    CALL dyn_ldf(kstp)
    IF (ln_zdfosm) CALL dyn_osm(kstp)
    CALL dyn_hpg(kstp)
    CALL dyn_spg(kstp)
    IF (ln_dynspg_ts) THEN
      CALL div_hor(kstp)
      IF (.NOT. ln_linssh) CALL dom_vvl_sf_nxt(kstp, kcall = 2)
      CALL wzv(kstp)
    END IF
    CALL dyn_zdf(kstp)
    IF (ln_diurnal) CALL stp_diurnal(kstp)
    IF (lk_floats) CALL flo_stp(kstp)
    IF (ln_diacfl) CALL dia_cfl(kstp)
    IF (lk_diahth) CALL dia_hth(kstp)
    IF (lk_diadct) CALL dia_dct(kstp)
    CALL dia_ar5(kstp)
    IF (lk_diaharm) CALL dia_harm(kstp)
    CALL dia_wri(kstp)
    IF (ln_crs) CALL crs_fld(kstp)
    CALL ProfileEnd(psy_profile1)
    !$ACC KERNELS
    DO jl=LBOUND(tsa,4),UBOUND(tsa,4)
      DO jk=LBOUND(tsa,3),UBOUND(tsa,3)
        DO jj=LBOUND(tsa,2),UBOUND(tsa,2)
          DO ji=LBOUND(tsa,1),UBOUND(tsa,1)
            tsa(ji, jj, jk, jl) = 0._wp
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !$ACC END KERNELS
    CALL ProfileStart('stp', 'r2', psy_profile2)
    IF (lk_asminc .AND. ln_asmiau .AND. ln_trainc) CALL tra_asm_inc(kstp)
    CALL tra_sbc(kstp)
    IF (ln_traqsr) CALL tra_qsr(kstp)
    IF (ln_trabbc) CALL tra_bbc(kstp)
    IF (ln_trabbl) CALL tra_bbl(kstp)
    IF (ln_tradmp) CALL tra_dmp(kstp)
    IF (ln_bdy) CALL bdy_tra_dmp(kstp)
    CALL tra_adv(kstp)
    IF (ln_zdfosm) CALL tra_osm(kstp)
    IF (lrst_oce .AND. ln_zdfosm) CALL osm_rst(kstp, 'WRITE')
    CALL tra_ldf(kstp)
    IF (ln_diaptr) CALL dia_ptr
    CALL tra_zdf(kstp)
    IF (ln_zdfnpc) CALL tra_npc(kstp)
    CALL tra_nxt(kstp)
    CALL dyn_nxt(kstp)
    CALL ssh_swp(kstp)
    IF (.NOT. ln_linssh) CALL dom_vvl_sf_swp(kstp)
    IF (ln_diahsb) CALL dia_hsb(kstp)
    IF (lrst_oce) CALL rst_write(kstp)
    IF (ln_sto_eos) CALL sto_rst_write(kstp)
    IF (ln_diaobs) CALL dia_obs(kstp)
    CALL stp_ctl(kstp, indic)
    IF (kstp == nit000) THEN
      CALL iom_close(numror)
      IF (lwm) CALL FLUSH(numond)
      IF (lwm .AND. numoni /= - 1) CALL FLUSH(numoni)
    END IF
    IF (lk_oasis) CALL sbc_cpl_snd(kstp)
    IF (ln_timing) CALL timing_stop('stp')
    CALL ProfileEnd(psy_profile2)
  END SUBROUTINE stp
END MODULE step
