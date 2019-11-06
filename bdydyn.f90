MODULE bdydyn
  USE oce
  USE dom_oce
  USE bdy_oce
  USE bdydyn2d
  USE bdydyn3d
  USE lbclnk
  USE in_out_manager
  USE domvvl
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: bdy_dyn
  CONTAINS
  SUBROUTINE bdy_dyn(kt, dyn3d_only)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    INTEGER, INTENT(IN) :: kt
    LOGICAL, INTENT(IN), OPTIONAL :: dyn3d_only
    INTEGER :: jk, ii, ij, ib_bdy, ib, igrd
    LOGICAL :: ll_dyn2d, ll_dyn3d, ll_orlanski
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: pua2d, pva2d
    TYPE(ProfileData), SAVE :: psy_profile0
    TYPE(ProfileData), SAVE :: psy_profile1
    CALL ProfileStart('bdy_dyn', 'r0', psy_profile0)
    ll_dyn2d = .TRUE.
    ll_dyn3d = .TRUE.
    IF (PRESENT(dyn3d_only)) THEN
      IF (dyn3d_only) ll_dyn2d = .FALSE.
    END IF
    ll_orlanski = .FALSE.
    DO ib_bdy = 1, nb_bdy
      IF (cn_dyn2d(ib_bdy) == 'orlanski' .OR. cn_dyn2d(ib_bdy) == 'orlanski_npo' .OR. cn_dyn3d(ib_bdy) == 'orlanski' .OR. cn_dyn3d(ib_bdy) == 'orlanski_npo') ll_orlanski = .TRUE.
    END DO
    CALL ProfileEnd(psy_profile0)
    !$ACC KERNELS
    pua2d(:, :) = 0._wp
    pva2d(:, :) = 0._wp
    !$ACC END KERNELS
    DO jk = 1, jpkm1
      !$ACC KERNELS
      pua2d(:, :) = pua2d(:, :) + e3u_a(:, :, jk) * ua(:, :, jk) * umask(:, :, jk)
      pva2d(:, :) = pva2d(:, :) + e3v_a(:, :, jk) * va(:, :, jk) * vmask(:, :, jk)
      !$ACC END KERNELS
    END DO
    !$ACC KERNELS
    pua2d(:, :) = pua2d(:, :) * r1_hu_a(:, :)
    pva2d(:, :) = pva2d(:, :) * r1_hv_a(:, :)
    !$ACC END KERNELS
    DO jk = 1, jpkm1
      !$ACC KERNELS
      ua(:, :, jk) = (ua(:, :, jk) - pua2d(:, :)) * umask(:, :, jk)
      va(:, :, jk) = (va(:, :, jk) - pva2d(:, :)) * vmask(:, :, jk)
      !$ACC END KERNELS
    END DO
    IF (ll_orlanski) THEN
      DO jk = 1, jpkm1
        !$ACC KERNELS
        ub(:, :, jk) = (ub(:, :, jk) - ub_b(:, :)) * umask(:, :, jk)
        vb(:, :, jk) = (vb(:, :, jk) - vb_b(:, :)) * vmask(:, :, jk)
        !$ACC END KERNELS
      END DO
    END IF
    CALL ProfileStart('bdy_dyn', 'r1', psy_profile1)
    IF (ll_dyn2d) CALL bdy_dyn2d(kt, pua2d, pva2d, ub_b, vb_b, r1_hu_a(:, :), r1_hv_a(:, :), ssha)
    IF (ll_dyn3d) CALL bdy_dyn3d(kt)
    CALL ProfileEnd(psy_profile1)
    DO jk = 1, jpkm1
      !$ACC KERNELS
      ua(:, :, jk) = (ua(:, :, jk) + pua2d(:, :)) * umask(:, :, jk)
      va(:, :, jk) = (va(:, :, jk) + pva2d(:, :)) * vmask(:, :, jk)
      !$ACC END KERNELS
    END DO
    IF (ll_orlanski) THEN
      DO jk = 1, jpkm1
        !$ACC KERNELS
        ub(:, :, jk) = (ub(:, :, jk) + ub_b(:, :)) * umask(:, :, jk)
        vb(:, :, jk) = (vb(:, :, jk) + vb_b(:, :)) * vmask(:, :, jk)
        !$ACC END KERNELS
      END DO
    END IF
  END SUBROUTINE bdy_dyn
END MODULE bdydyn