MODULE crslbclnk
  USE par_kind, ONLY: wp
  USE dom_oce
  USE crs
  USE lbclnk
  USE in_out_manager
  INTERFACE crs_lbc_lnk
    MODULE PROCEDURE crs_lbc_lnk_3d, crs_lbc_lnk_2d
  END INTERFACE
  PUBLIC :: crs_lbc_lnk
  CONTAINS
  SUBROUTINE crs_lbc_lnk_3d(pt3d1, cd_type1, psgn, cd_mpp, pval)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    CHARACTER(LEN = 1), INTENT(IN   ) :: cd_type1
    REAL(KIND = wp), INTENT(IN   ) :: psgn
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs, jpk), INTENT(INOUT) :: pt3d1
    REAL(KIND = wp), OPTIONAL, INTENT(IN   ) :: pval
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN   ) :: cd_mpp
    LOGICAL :: ll_grid_crs
    REAL(KIND = wp) :: zval
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('crs_lbc_lnk_3d', 'r0', psy_profile0)
    ll_grid_crs = (jpi == jpi_crs)
    IF (PRESENT(pval)) THEN
      zval = pval
    ELSE
      zval = 0._wp
    END IF
    IF (.NOT. ll_grid_crs) CALL dom_grid_crs
    IF (PRESENT(cd_mpp)) THEN
      CALL lbc_lnk(pt3d1, cd_type1, psgn, cd_mpp, pval = zval)
    ELSE
      CALL lbc_lnk(pt3d1, cd_type1, psgn, pval = zval)
    END IF
    IF (.NOT. ll_grid_crs) CALL dom_grid_glo
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE crs_lbc_lnk_3d
  SUBROUTINE crs_lbc_lnk_2d(pt2d, cd_type, psgn, cd_mpp, pval)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    CHARACTER(LEN = 1), INTENT(IN   ) :: cd_type
    REAL(KIND = wp), INTENT(IN   ) :: psgn
    REAL(KIND = wp), DIMENSION(jpi_crs, jpj_crs), INTENT(INOUT) :: pt2d
    REAL(KIND = wp), OPTIONAL, INTENT(IN   ) :: pval
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN   ) :: cd_mpp
    LOGICAL :: ll_grid_crs
    REAL(KIND = wp) :: zval
    TYPE(ProfileData), SAVE :: psy_profile0
    CALL ProfileStart('crs_lbc_lnk_2d', 'r0', psy_profile0)
    ll_grid_crs = (jpi == jpi_crs)
    IF (PRESENT(pval)) THEN
      zval = pval
    ELSE
      zval = 0._wp
    END IF
    IF (.NOT. ll_grid_crs) CALL dom_grid_crs
    IF (PRESENT(cd_mpp)) THEN
      CALL lbc_lnk(pt2d, cd_type, psgn, cd_mpp, pval = zval)
    ELSE
      CALL lbc_lnk(pt2d, cd_type, psgn, pval = zval)
    END IF
    IF (.NOT. ll_grid_crs) CALL dom_grid_glo
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE crs_lbc_lnk_2d
END MODULE crslbclnk