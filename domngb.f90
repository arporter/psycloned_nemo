MODULE domngb
  USE dom_oce
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dom_ngb
  CONTAINS
  SUBROUTINE dom_ngb(plon, plat, kii, kjj, cdgrid, kkk)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(IN) :: plon, plat
    INTEGER, INTENT(OUT) :: kii, kjj
    INTEGER, INTENT(IN), OPTIONAL :: kkk
    CHARACTER(LEN = 1), INTENT(IN) :: cdgrid
    INTEGER :: ik
    INTEGER, DIMENSION(2) :: iloc
    REAL(KIND = wp) :: zlon, zmini
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zglam, zgphi, zmask, zdist
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data1
    !$ACC KERNELS
    zmask(:, :) = 0._wp
    ik = 1
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('dom_ngb', 'r0', 0, 0)
    IF (PRESENT(kkk)) ik = kkk
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    SELECT CASE (cdgrid)
    CASE ('U')
      zglam(:, :) = glamu(:, :)
      zgphi(:, :) = gphiu(:, :)
      zmask(nldi : nlei, nldj : nlej) = umask(nldi : nlei, nldj : nlej, ik)
    CASE ('V')
      zglam(:, :) = glamv(:, :)
      zgphi(:, :) = gphiv(:, :)
      zmask(nldi : nlei, nldj : nlej) = vmask(nldi : nlei, nldj : nlej, ik)
    CASE ('F')
      zglam(:, :) = glamf(:, :)
      zgphi(:, :) = gphif(:, :)
      zmask(nldi : nlei, nldj : nlej) = fmask(nldi : nlei, nldj : nlej, ik)
    CASE DEFAULT
      zglam(:, :) = glamt(:, :)
      zgphi(:, :) = gphit(:, :)
      zmask(nldi : nlei, nldj : nlej) = tmask(nldi : nlei, nldj : nlej, ik)
    END SELECT
    zlon = MOD(plon + 720., 360.)
    zglam(:, :) = MOD(zglam(:, :) + 720., 360.)
    IF (zlon > 270.) zlon = zlon - 360.
    IF (zlon < 90.) WHERE (zglam(:, :) > 180.) zglam(:, :) = zglam(:, :) - 360.
    zglam(:, :) = zglam(:, :) - zlon
    zgphi(:, :) = zgphi(:, :) - plat
    zdist(:, :) = zglam(:, :) * zglam(:, :) + zgphi(:, :) * zgphi(:, :)
    !$ACC END KERNELS
    CALL profile_psy_data1 % PreStart('dom_ngb', 'r1', 0, 0)
    IF (lk_mpp) THEN
      CALL mpp_minloc('domngb', zdist(:, :), zmask, zmini, iloc)
      kii = iloc(1)
      kjj = iloc(2)
    ELSE
      iloc(:) = MINLOC(zdist(:, :), mask = zmask(:, :) == 1.E0)
      kii = iloc(1) + nimpp - 1
      kjj = iloc(2) + njmpp - 1
    END IF
    CALL profile_psy_data1 % PostEnd
  END SUBROUTINE dom_ngb
END MODULE domngb