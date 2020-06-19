MODULE sbc_ice
  USE lib_mpp
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: sbc_ice_alloc
  LOGICAL, PUBLIC, PARAMETER :: lk_si3 = .FALSE.
  LOGICAL, PUBLIC, PARAMETER :: lk_cice = .FALSE.
  REAL(KIND = wp), PUBLIC, PARAMETER :: cldf_ice = 0.81
  INTEGER, PUBLIC, PARAMETER :: jpl = 1
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: u_ice, v_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: tn_ice, alb_ice, qns_ice, dqns_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: a_i
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: emp_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: qsr_ice
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: h_i, h_s
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :, :) :: topmelt, botmelt
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: snwice_mass
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: snwice_mass_b
  REAL(KIND = wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:, :) :: snwice_fmass
  CONTAINS
  INTEGER FUNCTION sbc_ice_alloc()
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: ierr(1)
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    !$ACC KERNELS
    ierr(:) = 0
    !$ACC END KERNELS
    CALL profile_psy_data0 % PreStart('sbc_ice_alloc', 'r0', 0, 0)
    ALLOCATE(snwice_mass(jpi, jpj), snwice_mass_b(jpi, jpj), snwice_fmass(jpi, jpj), STAT = ierr(1))
    sbc_ice_alloc = MAXVAL(ierr)
    IF (lk_mpp) CALL mpp_sum(sbc_ice_alloc)
    IF (sbc_ice_alloc > 0) CALL ctl_warn('sbc_ice_alloc: allocation of arrays failed')
    CALL profile_psy_data0 % PostEnd
  END FUNCTION sbc_ice_alloc
END MODULE sbc_ice