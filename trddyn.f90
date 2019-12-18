MODULE trddyn
  USE oce
  USE dom_oce
  USE phycst
  USE sbc_oce
  USE zdf_oce
  USE trd_oce
  USE trdken
  USE trdglo
  USE trdvor
  USE trdmxl
  USE in_out_manager
  USE lbclnk
  USE iom
  USE lib_mpp
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: trd_dyn
  CONTAINS
  SUBROUTINE trd_dyn(putrd, pvtrd, ktrd, kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: putrd, pvtrd
    INTEGER, INTENT(IN   ) :: ktrd
    INTEGER, INTENT(IN   ) :: kt
    TYPE(ProfileData), SAVE :: psy_profile0
    !$ACC KERNELS
    putrd(:, :, :) = putrd(:, :, :) * umask(:, :, :)
    pvtrd(:, :, :) = pvtrd(:, :, :) * vmask(:, :, :)
    !$ACC END KERNELS
    CALL ProfileStart('trd_dyn', 'r0', psy_profile0)
    IF (ln_dyn_trd) CALL trd_dyn_iom(putrd, pvtrd, ktrd, kt)
    IF (ln_glo_trd) CALL trd_glo(putrd, pvtrd, ktrd, 'DYN', kt)
    IF (ln_KE_trd) CALL trd_ken(putrd, pvtrd, ktrd, kt)
    IF (ln_vor_trd) CALL trd_vor(putrd, pvtrd, ktrd, kt)
    CALL ProfileEnd(psy_profile0)
  END SUBROUTINE trd_dyn
  SUBROUTINE trd_dyn_iom(putrd, pvtrd, ktrd, kt)
    USE profile_mod, ONLY: ProfileData, ProfileStart, ProfileEnd
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: putrd, pvtrd
    INTEGER, INTENT(IN   ) :: ktrd
    INTEGER, INTENT(IN   ) :: kt
    INTEGER :: ji, jj, jk
    INTEGER :: ikbu, ikbv
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: z2dx, z2dy
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :, :) :: z3dx, z3dy
    TYPE(ProfileData), SAVE :: psy_profile0
    TYPE(ProfileData), SAVE :: psy_profile1
    TYPE(ProfileData), SAVE :: psy_profile2
    TYPE(ProfileData), SAVE :: psy_profile3
    TYPE(ProfileData), SAVE :: psy_profile4
    TYPE(ProfileData), SAVE :: psy_profile5
    TYPE(ProfileData), SAVE :: psy_profile6
    TYPE(ProfileData), SAVE :: psy_profile7
    TYPE(ProfileData), SAVE :: psy_profile8
    TYPE(ProfileData), SAVE :: psy_profile9
    TYPE(ProfileData), SAVE :: psy_profile10
    SELECT CASE (ktrd)
    CASE (jpdyn_hpg)
      CALL ProfileStart('trd_dyn_iom', 'r0', psy_profile0)
      CALL iom_put("utrd_hpg", putrd)
      CALL iom_put("vtrd_hpg", pvtrd)
      CALL ProfileEnd(psy_profile0)
    CASE (jpdyn_spg)
      CALL ProfileStart('trd_dyn_iom', 'r1', psy_profile1)
      CALL iom_put("utrd_spg", putrd)
      CALL iom_put("vtrd_spg", pvtrd)
      CALL ProfileEnd(psy_profile1)
    CASE (jpdyn_pvo)
      CALL ProfileStart('trd_dyn_iom', 'r2', psy_profile2)
      CALL iom_put("utrd_pvo", putrd)
      CALL iom_put("vtrd_pvo", pvtrd)
      CALL ProfileEnd(psy_profile2)
    CASE (jpdyn_rvo)
      CALL ProfileStart('trd_dyn_iom', 'r3', psy_profile3)
      CALL iom_put("utrd_rvo", putrd)
      CALL iom_put("vtrd_rvo", pvtrd)
      CALL ProfileEnd(psy_profile3)
    CASE (jpdyn_keg)
      CALL ProfileStart('trd_dyn_iom', 'r4', psy_profile4)
      CALL iom_put("utrd_keg", putrd)
      CALL iom_put("vtrd_keg", pvtrd)
      ALLOCATE(z3dx(jpi, jpj, jpk), z3dy(jpi, jpj, jpk))
      CALL ProfileEnd(psy_profile4)
      !$ACC KERNELS
      z3dx(:, :, :) = 0._wp
      z3dy(:, :, :) = 0._wp
      DO jk = 1, jpkm1
        !$ACC LOOP INDEPENDENT COLLAPSE(2)
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            z3dx(ji, jj, jk) = un(ji, jj, jk) * (un(ji + 1, jj, jk) - un(ji - 1, jj, jk)) / (2._wp * e1u(ji, jj))
            z3dy(ji, jj, jk) = vn(ji, jj, jk) * (vn(ji, jj + 1, jk) - vn(ji, jj - 1, jk)) / (2._wp * e2v(ji, jj))
          END DO
        END DO
      END DO
      !$ACC END KERNELS
      CALL ProfileStart('trd_dyn_iom', 'r5', psy_profile5)
      CALL lbc_lnk_multi(z3dx, 'U', - 1., z3dy, 'V', - 1.)
      CALL iom_put("utrd_udx", z3dx)
      CALL iom_put("vtrd_vdy", z3dy)
      DEALLOCATE(z3dx, z3dy)
      CALL ProfileEnd(psy_profile5)
    CASE (jpdyn_zad)
      CALL ProfileStart('trd_dyn_iom', 'r6', psy_profile6)
      CALL iom_put("utrd_zad", putrd)
      CALL iom_put("vtrd_zad", pvtrd)
      CALL ProfileEnd(psy_profile6)
    CASE (jpdyn_ldf)
      CALL ProfileStart('trd_dyn_iom', 'r7', psy_profile7)
      CALL iom_put("utrd_ldf", putrd)
      CALL iom_put("vtrd_ldf", pvtrd)
      CALL ProfileEnd(psy_profile7)
    CASE (jpdyn_zdf)
      CALL ProfileStart('trd_dyn_iom', 'r8', psy_profile8)
      CALL iom_put("utrd_zdf", putrd)
      CALL iom_put("vtrd_zdf", pvtrd)
      ALLOCATE(z2dx(jpi, jpj), z2dy(jpi, jpj))
      CALL ProfileEnd(psy_profile8)
      !$ACC KERNELS
      z2dx(:, :) = (utau_b(:, :) + utau(:, :)) / (e3u_n(:, :, 1) * rau0)
      z2dy(:, :) = (vtau_b(:, :) + vtau(:, :)) / (e3v_n(:, :, 1) * rau0)
      !$ACC END KERNELS
      CALL ProfileStart('trd_dyn_iom', 'r9', psy_profile9)
      CALL iom_put("utrd_tau", z2dx)
      CALL iom_put("vtrd_tau", z2dy)
      DEALLOCATE(z2dx, z2dy)
      CALL ProfileEnd(psy_profile9)
    CASE (jpdyn_bfr)
      CALL ProfileStart('trd_dyn_iom', 'r10', psy_profile10)
      CALL iom_put("utrd_bfr", putrd)
      CALL iom_put("vtrd_bfr", pvtrd)
      CALL ProfileEnd(psy_profile10)
    CASE (jpdyn_atf)
      CALL iom_put("utrd_atf", putrd)
      CALL iom_put("vtrd_atf", pvtrd)
    END SELECT
  END SUBROUTINE trd_dyn_iom
END MODULE trddyn