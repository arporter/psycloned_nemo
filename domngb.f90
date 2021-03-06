MODULE domngb
  USE dom_oce
  USE in_out_manager
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dom_ngb
  CONTAINS
  SUBROUTINE dom_ngb(plon, plat, kii, kjj, cdgrid, kkk)
    REAL(KIND = wp), INTENT(IN   ) :: plon, plat
    INTEGER, INTENT(  OUT) :: kii, kjj
    INTEGER, INTENT(IN   ), OPTIONAL :: kkk
    CHARACTER(LEN = 1), INTENT(IN   ) :: cdgrid
    INTEGER :: ik
    INTEGER, DIMENSION(2) :: iloc
    REAL(KIND = wp) :: zlon, zmini
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zglam, zgphi, zmask, zdist
    !$ACC KERNELS
    zmask(:, :) = 0._wp
    ik = 1
    !$ACC END KERNELS
    IF (PRESENT(kkk)) ik = kkk
    SELECT CASE (cdgrid)
    CASE ('U')
      !$ACC KERNELS
      zglam(:, :) = glamu(:, :)
      zgphi(:, :) = gphiu(:, :)
      zmask(nldi : nlei, nldj : nlej) = umask(nldi : nlei, nldj : nlej, ik)
      !$ACC END KERNELS
    CASE ('V')
      !$ACC KERNELS
      zglam(:, :) = glamv(:, :)
      zgphi(:, :) = gphiv(:, :)
      zmask(nldi : nlei, nldj : nlej) = vmask(nldi : nlei, nldj : nlej, ik)
      !$ACC END KERNELS
    CASE ('F')
      !$ACC KERNELS
      zglam(:, :) = glamf(:, :)
      zgphi(:, :) = gphif(:, :)
      zmask(nldi : nlei, nldj : nlej) = fmask(nldi : nlei, nldj : nlej, ik)
      !$ACC END KERNELS
    CASE DEFAULT
      !$ACC KERNELS
      zglam(:, :) = glamt(:, :)
      zgphi(:, :) = gphit(:, :)
      zmask(nldi : nlei, nldj : nlej) = tmask(nldi : nlei, nldj : nlej, ik)
      !$ACC END KERNELS
    END SELECT
    !$ACC KERNELS
    zlon = MOD(plon + 720., 360.)
    zglam(:, :) = MOD(zglam(:, :) + 720., 360.)
    IF (zlon > 270.) zlon = zlon - 360.
    !$ACC END KERNELS
    IF (zlon < 90.) WHERE (zglam(:, :) > 180.) zglam(:, :) = zglam(:, :) - 360.
    !$ACC KERNELS
    zglam(:, :) = zglam(:, :) - zlon
    zgphi(:, :) = zgphi(:, :) - plat
    zdist(:, :) = zglam(:, :) * zglam(:, :) + zgphi(:, :) * zgphi(:, :)
    !$ACC END KERNELS
    IF (lk_mpp) THEN
      CALL mpp_minloc(zdist(:, :), zmask, zmini, kii, kjj)
    ELSE
      iloc(:) = MINLOC(zdist(:, :), mask = zmask(:, :) == 1.E0)
      kii = iloc(1) + nimpp - 1
      kjj = iloc(2) + njmpp - 1
    END IF
  END SUBROUTINE dom_ngb
END MODULE domngb