MODULE diatmb
  USE oce
  USE dom_oce
  USE in_out_manager
  USE iom
  IMPLICIT NONE
  PRIVATE
  LOGICAL, PUBLIC :: ln_diatmb
  PUBLIC :: dia_tmb_init
  PUBLIC :: dia_tmb
  CONTAINS
  SUBROUTINE dia_tmb_init
    INTEGER :: ios
    NAMELIST /nam_diatmb/ ln_diatmb
    REWIND(UNIT = numnam_ref)
    READ(numnam_ref, nam_diatmb, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'nam_diatmb in reference namelist', lwp)
    REWIND(UNIT = numnam_cfg)
    READ(numnam_cfg, nam_diatmb, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'nam_diatmb in configuration namelist', lwp)
    IF (lwm) WRITE(numond, nam_diatmb)
    IF (lwp) THEN
      WRITE(numout, FMT = *)
      WRITE(numout, FMT = *) 'dia_tmb_init : Output Top, Middle, Bottom Diagnostics'
      WRITE(numout, FMT = *) '~~~~~~~~~~~~'
      WRITE(numout, FMT = *) '   Namelist nam_diatmb : set tmb outputs '
      WRITE(numout, FMT = *) '      Switch for TMB diagnostics (T) or not (F)  ln_diatmb  = ', ln_diatmb
    END IF
  END SUBROUTINE dia_tmb_init
  SUBROUTINE dia_calctmb(pfield, ptmb)
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk), INTENT(IN ) :: pfield
    REAL(KIND = wp), DIMENSION(jpi, jpj, 3), INTENT( OUT) :: ptmb
    INTEGER :: ji, jj
    INTEGER :: itop, imid, ibot
    REAL(KIND = wp) :: zmdi = 1.E+20_wp
    !$ACC KERNELS
    DO jj = 1, jpj
      DO ji = 1, jpi
        itop = mikt(ji, jj)
        ibot = mbkt(ji, jj)
        imid = itop + (ibot - itop + 1) / 2
        ptmb(ji, jj, 1) = pfield(ji, jj, itop) * tmask(ji, jj, itop) + zmdi * (1._wp - tmask(ji, jj, itop))
        ptmb(ji, jj, 2) = pfield(ji, jj, imid) * tmask(ji, jj, imid) + zmdi * (1._wp - tmask(ji, jj, imid))
        ptmb(ji, jj, 3) = pfield(ji, jj, ibot) * tmask(ji, jj, ibot) + zmdi * (1._wp - tmask(ji, jj, ibot))
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE dia_calctmb
  SUBROUTINE dia_tmb
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    REAL(KIND = wp) :: zmdi = 1.E+20
    REAL(KIND = wp), DIMENSION(jpi, jpj, 3) :: zwtmb
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('dia_tmb', 'r0', psy_profile0)
    CALL dia_calctmb(tsn(:, :, :, jp_tem), zwtmb)
    CALL iom_put("sshnmasked", sshn(:, :) * tmask(:, :, 1) + zmdi * (1.0 - tmask(:, :, 1)))
    CALL iom_put("top_temp", zwtmb(:, :, 1))
    CALL iom_put("mid_temp", zwtmb(:, :, 2))
    CALL iom_put("bot_temp", zwtmb(:, :, 3))
    CALL dia_calctmb(tsn(:, :, :, jp_sal), zwtmb)
    CALL iom_put("top_sal", zwtmb(:, :, 1))
    CALL iom_put("mid_sal", zwtmb(:, :, 2))
    CALL iom_put("bot_sal", zwtmb(:, :, 3))
    CALL dia_calctmb(un(:, :, :), zwtmb)
    CALL iom_put("top_u", zwtmb(:, :, 1))
    CALL iom_put("mid_u", zwtmb(:, :, 2))
    CALL iom_put("bot_u", zwtmb(:, :, 3))
    CALL dia_calctmb(vn(:, :, :), zwtmb)
    CALL iom_put("top_v", zwtmb(:, :, 1))
    CALL iom_put("mid_v", zwtmb(:, :, 2))
    CALL iom_put("bot_v", zwtmb(:, :, 3))
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE dia_tmb
END MODULE diatmb