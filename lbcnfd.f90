MODULE lbcnfd
  USE dom_oce
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  INTERFACE lbc_nfd
    MODULE PROCEDURE lbc_nfd_2d, lbc_nfd_3d, lbc_nfd_4d
    MODULE PROCEDURE lbc_nfd_2d_ptr, lbc_nfd_3d_ptr, lbc_nfd_4d_ptr
    MODULE PROCEDURE lbc_nfd_2d_ext
  END INTERFACE
  INTERFACE lbc_nfd_nogather
    MODULE PROCEDURE lbc_nfd_nogather_4d
  END INTERFACE
  TYPE, PUBLIC :: PTR_2D
    REAL(KIND = wp), DIMENSION(:, :), POINTER :: pt2d
  END TYPE PTR_2D
  TYPE, PUBLIC :: PTR_3D
    REAL(KIND = wp), DIMENSION(:, :, :), POINTER :: pt3d
  END TYPE PTR_3D
  TYPE, PUBLIC :: PTR_4D
    REAL(KIND = wp), DIMENSION(:, :, :, :), POINTER :: pt4d
  END TYPE PTR_4D
  PUBLIC :: lbc_nfd
  PUBLIC :: lbc_nfd_nogather
  INTEGER, PUBLIC, PARAMETER :: jpmaxngh = 3
  INTEGER, PUBLIC :: nsndto, nfsloop, nfeloop
  INTEGER, PUBLIC, DIMENSION(jpmaxngh) :: isendto
  CONTAINS
  SUBROUTINE lbc_nfd_2d(ptab, cd_nat, psgn)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :)
    CHARACTER(LEN = 1), INTENT(IN   ) :: cd_nat
    REAL(KIND = wp), INTENT(IN   ) :: psgn
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
    ipk = 1
    ipl = 1
    ipf = 1
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          !$ACC KERNELS
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipj) = psgn * ptab(ijt, ipj - 2)
          END DO
          ptab(1, ipj) = psgn * ptab(3, ipj - 2)
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipjm1) = psgn * ptab(ijt, ipjm1)
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipj) = psgn * ptab(iju, ipj - 2)
          END DO
          ptab(1, ipj) = psgn * ptab(2, ipj - 2)
          ptab(jpiglo, ipj) = psgn * ptab(jpiglo - 1, ipj - 2)
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipjm1) = psgn * ptab(iju, ipjm1)
          END DO
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipj - 1) = psgn * ptab(ijt, ipj - 2)
            ptab(ji, ipj) = psgn * ptab(ijt, ipj - 3)
          END DO
          ptab(1, ipj) = psgn * ptab(3, ipj - 3)
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipj - 1) = psgn * ptab(iju, ipj - 2)
            ptab(ji, ipj) = psgn * ptab(iju, ipj - 3)
          END DO
          ptab(1, ipj) = psgn * ptab(2, ipj - 3)
          ptab(jpiglo, ipj) = psgn * ptab(jpiglo - 1, ipj - 3)
          !$ACC END KERNELS
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          !$ACC KERNELS
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipj) = psgn * ptab(ijt, ipj - 1)
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipj) = psgn * ptab(iju, ipj - 1)
          END DO
          ptab(jpiglo, ipj) = psgn * ptab(1, ipj - 1)
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipj) = psgn * ptab(ijt, ipj - 2)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipjm1) = psgn * ptab(ijt, ipjm1)
          END DO
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipj) = psgn * ptab(iju, ipj - 2)
          END DO
          ptab(jpiglo, ipj) = psgn * ptab(1, ipj - 2)
          DO ji = jpiglo / 2 + 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipjm1) = psgn * ptab(iju, ipjm1)
          END DO
          !$ACC END KERNELS
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat)
        CASE ('T', 'U', 'V', 'W')
          !$ACC KERNELS
          ptab(:, 1) = 0._wp
          ptab(:, ipj) = 0._wp
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          ptab(:, ipj) = 0._wp
          !$ACC END KERNELS
        END SELECT
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_2d
  SUBROUTINE lbc_nfd_2d_ptr(ptab, cd_nat, psgn, kfld)
    INTEGER, INTENT(IN   ) :: kfld
    TYPE(PTR_2D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN   ) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN   ) :: psgn(:)
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
    ipk = 1
    ipl = 1
    ipf = kfld
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(ijt, ipj - 2)
          END DO
          ptab(jf) % pt2d(1, ipj) = psgn(jf) * ptab(jf) % pt2d(3, ipj - 2)
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt2d(ji, ipjm1) = psgn(jf) * ptab(jf) % pt2d(ijt, ipjm1)
          END DO
        CASE ('U')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(iju, ipj - 2)
          END DO
          ptab(jf) % pt2d(1, ipj) = psgn(jf) * ptab(jf) % pt2d(2, ipj - 2)
          ptab(jf) % pt2d(jpiglo, ipj) = psgn(jf) * ptab(jf) % pt2d(jpiglo - 1, ipj - 2)
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt2d(ji, ipjm1) = psgn(jf) * ptab(jf) % pt2d(iju, ipjm1)
          END DO
        CASE ('V')
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt2d(ji, ipj - 1) = psgn(jf) * ptab(jf) % pt2d(ijt, ipj - 2)
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(ijt, ipj - 3)
          END DO
          ptab(jf) % pt2d(1, ipj) = psgn(jf) * ptab(jf) % pt2d(3, ipj - 3)
        CASE ('F')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt2d(ji, ipj - 1) = psgn(jf) * ptab(jf) % pt2d(iju, ipj - 2)
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(iju, ipj - 3)
          END DO
          ptab(jf) % pt2d(1, ipj) = psgn(jf) * ptab(jf) % pt2d(2, ipj - 3)
          ptab(jf) % pt2d(jpiglo, ipj) = psgn(jf) * ptab(jf) % pt2d(jpiglo - 1, ipj - 3)
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(ijt, ipj - 1)
          END DO
        CASE ('U')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(iju, ipj - 1)
          END DO
          ptab(jf) % pt2d(jpiglo, ipj) = psgn(jf) * ptab(jf) % pt2d(1, ipj - 1)
        CASE ('V')
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(ijt, ipj - 2)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt2d(ji, ipjm1) = psgn(jf) * ptab(jf) % pt2d(ijt, ipjm1)
          END DO
        CASE ('F')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt2d(ji, ipj) = psgn(jf) * ptab(jf) % pt2d(iju, ipj - 2)
          END DO
          ptab(jf) % pt2d(jpiglo, ipj) = psgn(jf) * ptab(jf) % pt2d(1, ipj - 2)
          DO ji = jpiglo / 2 + 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt2d(ji, ipjm1) = psgn(jf) * ptab(jf) % pt2d(iju, ipjm1)
          END DO
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'U', 'V', 'W')
          ptab(jf) % pt2d(:, 1) = 0._wp
          ptab(jf) % pt2d(:, ipj) = 0._wp
        CASE ('F')
          ptab(jf) % pt2d(:, ipj) = 0._wp
        END SELECT
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_2d_ptr
  SUBROUTINE lbc_nfd_2d_ext(ptab, cd_nat, psgn, kextj)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, 1 - kextj :)
    CHARACTER(LEN = 1), INTENT(IN   ) :: cd_nat
    REAL(KIND = wp), INTENT(IN   ) :: psgn
    INTEGER, INTENT(IN   ) :: kextj
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
    ipk = 1
    ipl = 1
    ipf = 1
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          !$ACC KERNELS
          DO jh = 0, kextj
            DO ji = 2, jpiglo
              ijt = jpiglo - ji + 2
              ptab(ji, ipj + jh) = psgn * ptab(ijt, ipj - 2 - jh)
            END DO
            ptab(1, ipj + jh) = psgn * ptab(3, ipj - 2 - jh)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipjm1) = psgn * ptab(ijt, ipjm1)
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS
          DO jh = 0, kextj
            DO ji = 2, jpiglo - 1
              iju = jpiglo - ji + 1
              ptab(ji, ipj + jh) = psgn * ptab(iju, ipj - 2 - jh)
            END DO
            ptab(1, ipj + jh) = psgn * ptab(2, ipj - 2 - jh)
            ptab(jpiglo, ipj + jh) = psgn * ptab(jpiglo - 1, ipj - 2 - jh)
          END DO
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipjm1) = psgn * ptab(iju, ipjm1)
          END DO
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS
          DO jh = 0, kextj
            DO ji = 2, jpiglo
              ijt = jpiglo - ji + 2
              ptab(ji, ipj - 1 + jh) = psgn * ptab(ijt, ipj - 2 - jh)
              ptab(ji, ipj + jh) = psgn * ptab(ijt, ipj - 3 - jh)
            END DO
            ptab(1, ipj + jh) = psgn * ptab(3, ipj - 3 - jh)
          END DO
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          DO jh = 0, kextj
            DO ji = 1, jpiglo - 1
              iju = jpiglo - ji + 1
              ptab(ji, ipj - 1 + jh) = psgn * ptab(iju, ipj - 2 - jh)
              ptab(ji, ipj + jh) = psgn * ptab(iju, ipj - 3 - jh)
            END DO
          END DO
          DO jh = 0, kextj
            ptab(1, ipj + jh) = psgn * ptab(2, ipj - 3 - jh)
            ptab(jpiglo, ipj + jh) = psgn * ptab(jpiglo - 1, ipj - 3 - jh)
          END DO
          !$ACC END KERNELS
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          !$ACC KERNELS
          DO jh = 0, kextj
            DO ji = 1, jpiglo
              ijt = jpiglo - ji + 1
              ptab(ji, ipj + jh) = psgn * ptab(ijt, ipj - 1 - jh)
            END DO
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS
          DO jh = 0, kextj
            DO ji = 1, jpiglo - 1
              iju = jpiglo - ji
              ptab(ji, ipj + jh) = psgn * ptab(iju, ipj - 1 - jh)
            END DO
            ptab(jpiglo, ipj + jh) = psgn * ptab(1, ipj - 1 - jh)
          END DO
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS
          DO jh = 0, kextj
            DO ji = 1, jpiglo
              ijt = jpiglo - ji + 1
              ptab(ji, ipj + jh) = psgn * ptab(ijt, ipj - 2 - jh)
            END DO
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipjm1) = psgn * ptab(ijt, ipjm1)
          END DO
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          DO jh = 0, kextj
            DO ji = 1, jpiglo - 1
              iju = jpiglo - ji
              ptab(ji, ipj + jh) = psgn * ptab(iju, ipj - 2 - jh)
            END DO
            ptab(jpiglo, ipj + jh) = psgn * ptab(1, ipj - 2 - jh)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipjm1) = psgn * ptab(iju, ipjm1)
          END DO
          !$ACC END KERNELS
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat)
        CASE ('T', 'U', 'V', 'W')
          !$ACC KERNELS
          ptab(:, 1 : 1 - kextj) = 0._wp
          ptab(:, ipj : ipj + kextj) = 0._wp
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          ptab(:, ipj : ipj + kextj) = 0._wp
          !$ACC END KERNELS
        END SELECT
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_2d_ext
  SUBROUTINE lbc_nfd_3d(ptab, cd_nat, psgn)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :)
    CHARACTER(LEN = 1), INTENT(IN   ) :: cd_nat
    REAL(KIND = wp), INTENT(IN   ) :: psgn
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
    ipk = SIZE(ptab, 3)
    ipl = 1
    ipf = 1
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          !$ACC KERNELS
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipj, :) = psgn * ptab(ijt, ipj - 2, :)
          END DO
          ptab(1, ipj, :) = psgn * ptab(3, ipj - 2, :)
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipjm1, :) = psgn * ptab(ijt, ipjm1, :)
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipj, :) = psgn * ptab(iju, ipj - 2, :)
          END DO
          ptab(1, ipj, :) = psgn * ptab(2, ipj - 2, :)
          ptab(jpiglo, ipj, :) = psgn * ptab(jpiglo - 1, ipj - 2, :)
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipjm1, :) = psgn * ptab(iju, ipjm1, :)
          END DO
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipj - 1, :) = psgn * ptab(ijt, ipj - 2, :)
            ptab(ji, ipj, :) = psgn * ptab(ijt, ipj - 3, :)
          END DO
          ptab(1, ipj, :) = psgn * ptab(3, ipj - 3, :)
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipj - 1, :) = psgn * ptab(iju, ipj - 2, :)
            ptab(ji, ipj, :) = psgn * ptab(iju, ipj - 3, :)
          END DO
          ptab(1, ipj, :) = psgn * ptab(2, ipj - 3, :)
          ptab(jpiglo, ipj, :) = psgn * ptab(jpiglo - 1, ipj - 3, :)
          !$ACC END KERNELS
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          !$ACC KERNELS
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipj, :) = psgn * ptab(ijt, ipj - 1, :)
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipj, :) = psgn * ptab(iju, ipj - 1, :)
          END DO
          ptab(jpiglo, ipj, :) = psgn * ptab(1, ipj - 1, :)
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipj, :) = psgn * ptab(ijt, ipj - 2, :)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipjm1, :) = psgn * ptab(ijt, ipjm1, :)
          END DO
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipj, :) = psgn * ptab(iju, ipj - 2, :)
          END DO
          ptab(jpiglo, ipj, :) = psgn * ptab(1, ipj - 2, :)
          DO ji = jpiglo / 2 + 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipjm1, :) = psgn * ptab(iju, ipjm1, :)
          END DO
          !$ACC END KERNELS
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat)
        CASE ('T', 'U', 'V', 'W')
          !$ACC KERNELS
          ptab(:, 1, :) = 0._wp
          ptab(:, ipj, :) = 0._wp
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          ptab(:, ipj, :) = 0._wp
          !$ACC END KERNELS
        END SELECT
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_3d
  SUBROUTINE lbc_nfd_3d_ptr(ptab, cd_nat, psgn, kfld)
    INTEGER, INTENT(IN   ) :: kfld
    TYPE(PTR_3D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN   ) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN   ) :: psgn(:)
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
    ipk = SIZE(ptab(1) % pt3d, 3)
    ipl = 1
    ipf = kfld
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt3d(ji, ipj, :) = psgn(jf) * ptab(jf) % pt3d(ijt, ipj - 2, :)
          END DO
          ptab(jf) % pt3d(1, ipj, :) = psgn(jf) * ptab(jf) % pt3d(3, ipj - 2, :)
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt3d(ji, ipjm1, :) = psgn(jf) * ptab(jf) % pt3d(ijt, ipjm1, :)
          END DO
        CASE ('U')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt3d(ji, ipj, :) = psgn(jf) * ptab(jf) % pt3d(iju, ipj - 2, :)
          END DO
          ptab(jf) % pt3d(1, ipj, :) = psgn(jf) * ptab(jf) % pt3d(2, ipj - 2, :)
          ptab(jf) % pt3d(jpiglo, ipj, :) = psgn(jf) * ptab(jf) % pt3d(jpiglo - 1, ipj - 2, :)
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt3d(ji, ipjm1, :) = psgn(jf) * ptab(jf) % pt3d(iju, ipjm1, :)
          END DO
        CASE ('V')
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt3d(ji, ipj - 1, :) = psgn(jf) * ptab(jf) % pt3d(ijt, ipj - 2, :)
            ptab(jf) % pt3d(ji, ipj, :) = psgn(jf) * ptab(jf) % pt3d(ijt, ipj - 3, :)
          END DO
          ptab(jf) % pt3d(1, ipj, :) = psgn(jf) * ptab(jf) % pt3d(3, ipj - 3, :)
        CASE ('F')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt3d(ji, ipj - 1, :) = psgn(jf) * ptab(jf) % pt3d(iju, ipj - 2, :)
            ptab(jf) % pt3d(ji, ipj, :) = psgn(jf) * ptab(jf) % pt3d(iju, ipj - 3, :)
          END DO
          ptab(jf) % pt3d(1, ipj, :) = psgn(jf) * ptab(jf) % pt3d(2, ipj - 3, :)
          ptab(jf) % pt3d(jpiglo, ipj, :) = psgn(jf) * ptab(jf) % pt3d(jpiglo - 1, ipj - 3, :)
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt3d(ji, ipj, :) = psgn(jf) * ptab(jf) % pt3d(ijt, ipj - 1, :)
          END DO
        CASE ('U')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt3d(ji, ipj, :) = psgn(jf) * ptab(jf) % pt3d(iju, ipj - 1, :)
          END DO
          ptab(jf) % pt3d(jpiglo, ipj, :) = psgn(jf) * ptab(jf) % pt3d(1, ipj - 1, :)
        CASE ('V')
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt3d(ji, ipj, :) = psgn(jf) * ptab(jf) % pt3d(ijt, ipj - 2, :)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt3d(ji, ipjm1, :) = psgn(jf) * ptab(jf) % pt3d(ijt, ipjm1, :)
          END DO
        CASE ('F')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt3d(ji, ipj, :) = psgn(jf) * ptab(jf) % pt3d(iju, ipj - 2, :)
          END DO
          ptab(jf) % pt3d(jpiglo, ipj, :) = psgn(jf) * ptab(jf) % pt3d(1, ipj - 2, :)
          DO ji = jpiglo / 2 + 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt3d(ji, ipjm1, :) = psgn(jf) * ptab(jf) % pt3d(iju, ipjm1, :)
          END DO
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'U', 'V', 'W')
          ptab(jf) % pt3d(:, 1, :) = 0._wp
          ptab(jf) % pt3d(:, ipj, :) = 0._wp
        CASE ('F')
          ptab(jf) % pt3d(:, ipj, :) = 0._wp
        END SELECT
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_3d_ptr
  SUBROUTINE lbc_nfd_4d(ptab, cd_nat, psgn)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :, :)
    CHARACTER(LEN = 1), INTENT(IN   ) :: cd_nat
    REAL(KIND = wp), INTENT(IN   ) :: psgn
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
    ipk = SIZE(ptab, 3)
    ipl = SIZE(ptab, 4)
    ipf = 1
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          !$ACC KERNELS
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipj, :, :) = psgn * ptab(ijt, ipj - 2, :, :)
          END DO
          ptab(1, ipj, :, :) = psgn * ptab(3, ipj - 2, :, :)
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipjm1, :, :) = psgn * ptab(ijt, ipjm1, :, :)
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipj, :, :) = psgn * ptab(iju, ipj - 2, :, :)
          END DO
          ptab(1, ipj, :, :) = psgn * ptab(2, ipj - 2, :, :)
          ptab(jpiglo, ipj, :, :) = psgn * ptab(jpiglo - 1, ipj - 2, :, :)
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipjm1, :, :) = psgn * ptab(iju, ipjm1, :, :)
          END DO
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(ji, ipj - 1, :, :) = psgn * ptab(ijt, ipj - 2, :, :)
            ptab(ji, ipj, :, :) = psgn * ptab(ijt, ipj - 3, :, :)
          END DO
          ptab(1, ipj, :, :) = psgn * ptab(3, ipj - 3, :, :)
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(ji, ipj - 1, :, :) = psgn * ptab(iju, ipj - 2, :, :)
            ptab(ji, ipj, :, :) = psgn * ptab(iju, ipj - 3, :, :)
          END DO
          ptab(1, ipj, :, :) = psgn * ptab(2, ipj - 3, :, :)
          ptab(jpiglo, ipj, :, :) = psgn * ptab(jpiglo - 1, ipj - 3, :, :)
          !$ACC END KERNELS
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          !$ACC KERNELS
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipj, :, :) = psgn * ptab(ijt, ipj - 1, :, :)
          END DO
          !$ACC END KERNELS
        CASE ('U')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipj, :, :) = psgn * ptab(iju, ipj - 1, :, :)
          END DO
          ptab(jpiglo, ipj, :, :) = psgn * ptab(1, ipj - 1, :, :)
          !$ACC END KERNELS
        CASE ('V')
          !$ACC KERNELS
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipj, :, :) = psgn * ptab(ijt, ipj - 2, :, :)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(ji, ipjm1, :, :) = psgn * ptab(ijt, ipjm1, :, :)
          END DO
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipj, :, :) = psgn * ptab(iju, ipj - 2, :, :)
          END DO
          ptab(jpiglo, ipj, :, :) = psgn * ptab(1, ipj - 2, :, :)
          DO ji = jpiglo / 2 + 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(ji, ipjm1, :, :) = psgn * ptab(iju, ipjm1, :, :)
          END DO
          !$ACC END KERNELS
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat)
        CASE ('T', 'U', 'V', 'W')
          !$ACC KERNELS
          ptab(:, 1, :, :) = 0._wp
          ptab(:, ipj, :, :) = 0._wp
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          ptab(:, ipj, :, :) = 0._wp
          !$ACC END KERNELS
        END SELECT
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_4d
  SUBROUTINE lbc_nfd_4d_ptr(ptab, cd_nat, psgn, kfld)
    INTEGER, INTENT(IN   ) :: kfld
    TYPE(PTR_4D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN   ) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN   ) :: psgn(:)
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ipjm1
    ipk = SIZE(ptab(1) % pt4d, 3)
    ipl = SIZE(ptab(1) % pt4d, 4)
    ipf = kfld
    SELECT CASE (jpni)
    CASE (1)
      ipj = nlcj
    CASE DEFAULT
      ipj = 4
    END SELECT
    ipjm1 = ipj - 1
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipj - 2, :, :)
          END DO
          ptab(jf) % pt4d(1, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(3, ipj - 2, :, :)
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt4d(ji, ipjm1, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipjm1, :, :)
          END DO
        CASE ('U')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipj - 2, :, :)
          END DO
          ptab(jf) % pt4d(1, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(2, ipj - 2, :, :)
          ptab(jf) % pt4d(jpiglo, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(jpiglo - 1, ipj - 2, :, :)
          DO ji = jpiglo / 2, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt4d(ji, ipjm1, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipjm1, :, :)
          END DO
        CASE ('V')
          DO ji = 2, jpiglo
            ijt = jpiglo - ji + 2
            ptab(jf) % pt4d(ji, ipj - 1, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipj - 2, :, :)
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipj - 3, :, :)
          END DO
          ptab(jf) % pt4d(1, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(3, ipj - 3, :, :)
        CASE ('F')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji + 1
            ptab(jf) % pt4d(ji, ipj - 1, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipj - 2, :, :)
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipj - 3, :, :)
          END DO
          ptab(jf) % pt4d(1, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(2, ipj - 3, :, :)
          ptab(jf) % pt4d(jpiglo, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(jpiglo - 1, ipj - 3, :, :)
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'W')
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipj - 1, :, :)
          END DO
        CASE ('U')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipj - 1, :, :)
          END DO
          ptab(jf) % pt4d(jpiglo, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(1, ipj - 1, :, :)
        CASE ('V')
          DO ji = 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipj - 2, :, :)
          END DO
          DO ji = jpiglo / 2 + 1, jpiglo
            ijt = jpiglo - ji + 1
            ptab(jf) % pt4d(ji, ipjm1, :, :) = psgn(jf) * ptab(jf) % pt4d(ijt, ipjm1, :, :)
          END DO
        CASE ('F')
          DO ji = 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt4d(ji, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipj - 2, :, :)
          END DO
          ptab(jf) % pt4d(jpiglo, ipj, :, :) = psgn(jf) * ptab(jf) % pt4d(1, ipj - 2, :, :)
          DO ji = jpiglo / 2 + 1, jpiglo - 1
            iju = jpiglo - ji
            ptab(jf) % pt4d(ji, ipjm1, :, :) = psgn(jf) * ptab(jf) % pt4d(iju, ipjm1, :, :)
          END DO
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat(jf))
        CASE ('T', 'U', 'V', 'W')
          ptab(jf) % pt4d(:, 1, :, :) = 0._wp
          ptab(jf) % pt4d(:, ipj, :, :) = 0._wp
        CASE ('F')
          ptab(jf) % pt4d(:, ipj, :, :) = 0._wp
        END SELECT
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_4d_ptr
  SUBROUTINE lbc_nfd_nogather_4d(ptab, ptab2, cd_nat, psgn, kfld)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :, :)
    REAL(KIND = wp), INTENT(INOUT) :: ptab2(:, :, :, :)
    CHARACTER(LEN = 1), INTENT(IN   ) :: cd_nat
    REAL(KIND = wp), INTENT(IN   ) :: psgn
    INTEGER, OPTIONAL, INTENT(IN   ) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: ijt, iju, ijpj, ijpjm1, ijta, ijua, jia, startloop, endloop
    ipk = SIZE(ptab, 3)
    ipl = SIZE(ptab, 4)
    ipf = 1
    SELECT CASE (jpni)
    CASE (1)
      ijpj = nlcj
    CASE DEFAULT
      ijpj = 4
    END SELECT
    ijpjm1 = ijpj - 1
    DO jf = 1, ipf
      SELECT CASE (npolj)
      CASE (3, 4)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            startloop = 2
          END IF
          !$ACC KERNELS
          DO ji = startloop, nlci
            ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
            ptab(ji, ijpj, :, :) = psgn * ptab2(ijt, ijpj - 2, :, :)
          END DO
          !$ACC END KERNELS
          IF (nimpp == 1) THEN
            !$ACC KERNELS
            ptab(1, ijpj, :, :) = psgn * ptab(3, ijpj - 2, :, :)
            !$ACC END KERNELS
          END IF
          IF (nimpp >= jpiglo / 2 + 1) THEN
            startloop = 1
          ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
            startloop = jpiglo / 2 + 1 - nimpp + 1
          ELSE
            startloop = nlci + 1
          END IF
          IF (startloop <= nlci) THEN
            DO ji = startloop, nlci
              ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
              jia = ji + nimpp - 1
              ijta = jpiglo - jia + 2
              IF (ijta >= startloop + nimpp - 1 .AND. ijta < jia) THEN
                !$ACC KERNELS
                ptab(ji, ijpjm1, :, :) = psgn * ptab(ijta - nimpp + 1, ijpjm1, :, :)
                !$ACC END KERNELS
              ELSE
                !$ACC KERNELS
                ptab(ji, ijpjm1, :, :) = psgn * ptab2(ijt, ijpjm1, :, :)
                !$ACC END KERNELS
              END IF
            END DO
          END IF
        CASE ('U')
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          !$ACC KERNELS
          DO ji = 1, endloop
            iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
            ptab(ji, ijpj, :, :) = psgn * ptab2(iju, ijpj - 2, :, :)
          END DO
          !$ACC END KERNELS
          IF (nimpp .EQ. 1) THEN
            !$ACC KERNELS
            ptab(1, ijpj, :, :) = psgn * ptab(2, ijpj - 2, :, :)
            !$ACC END KERNELS
          END IF
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            !$ACC KERNELS
            ptab(nlci, ijpj, :, :) = psgn * ptab(nlci - 1, ijpj - 2, :, :)
            !$ACC END KERNELS
          END IF
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          IF (nimpp >= jpiglo / 2) THEN
            startloop = 1
          ELSE IF ((nimpp + nlci - 1 >= jpiglo / 2) .AND. (nimpp < jpiglo / 2)) THEN
            startloop = jpiglo / 2 - nimpp + 1
          ELSE
            startloop = endloop + 1
          END IF
          IF (startloop <= endloop) THEN
            DO ji = startloop, endloop
              iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
              jia = ji + nimpp - 1
              ijua = jpiglo - jia + 1
              IF (ijua >= startloop + nimpp - 1 .AND. ijua < jia) THEN
                !$ACC KERNELS
                ptab(ji, ijpjm1, :, :) = psgn * ptab(ijua - nimpp + 1, ijpjm1, :, :)
                !$ACC END KERNELS
              ELSE
                !$ACC KERNELS
                ptab(ji, ijpjm1, :, :) = psgn * ptab2(iju, ijpjm1, :, :)
                !$ACC END KERNELS
              END IF
            END DO
          END IF
        CASE ('V')
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            startloop = 2
          END IF
          !$ACC KERNELS
          DO ji = startloop, nlci
            ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
            ptab(ji, ijpj - 1, :, :) = psgn * ptab2(ijt, ijpj - 2, :, :)
            ptab(ji, ijpj, :, :) = psgn * ptab2(ijt, ijpj - 3, :, :)
          END DO
          !$ACC END KERNELS
          IF (nimpp .EQ. 1) THEN
            !$ACC KERNELS
            ptab(1, ijpj, :, :) = psgn * ptab(3, ijpj - 3, :, :)
            !$ACC END KERNELS
          END IF
        CASE ('F')
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          !$ACC KERNELS
          DO ji = 1, endloop
            iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
            ptab(ji, ijpj - 1, :, :) = psgn * ptab2(iju, ijpj - 2, :, :)
            ptab(ji, ijpj, :, :) = psgn * ptab2(iju, ijpj - 3, :, :)
          END DO
          !$ACC END KERNELS
          IF (nimpp .EQ. 1) THEN
            !$ACC KERNELS
            ptab(1, ijpj, :, :) = psgn * ptab(2, ijpj - 3, :, :)
            ptab(1, ijpj - 1, :, :) = psgn * ptab(2, ijpj - 2, :, :)
            !$ACC END KERNELS
          END IF
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            !$ACC KERNELS
            ptab(nlci, ijpj, :, :) = psgn * ptab(nlci - 1, ijpj - 3, :, :)
            ptab(nlci, ijpj - 1, :, :) = psgn * ptab(nlci - 1, ijpj - 2, :, :)
            !$ACC END KERNELS
          END IF
        CASE ('I')
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            !$ACC KERNELS
            startloop = 3
            ptab(2, ijpj, :, :) = psgn * ptab2(3, ijpjm1, :, :)
            !$ACC END KERNELS
          END IF
          !$ACC KERNELS
          DO ji = startloop, nlci
            iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 5
            ptab(ji, ijpj, :, :) = psgn * ptab2(iju, ijpjm1, :, :)
          END DO
          !$ACC END KERNELS
        END SELECT
      CASE (5, 6)
        SELECT CASE (cd_nat)
        CASE ('T', 'W')
          !$ACC KERNELS
          DO ji = 1, nlci
            ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
            ptab(ji, ijpj, :, :) = psgn * ptab2(ijt, ijpj - 1, :, :)
          END DO
          !$ACC END KERNELS
        CASE ('U')
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          !$ACC KERNELS
          DO ji = 1, endloop
            iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
            ptab(ji, ijpj, :, :) = psgn * ptab2(iju, ijpj - 1, :, :)
          END DO
          !$ACC END KERNELS
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            !$ACC KERNELS
            ptab(nlci, ijpj, :, :) = psgn * ptab2(1, ijpj - 1, :, :)
            !$ACC END KERNELS
          END IF
        CASE ('V')
          !$ACC KERNELS
          DO ji = 1, nlci
            ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
            ptab(ji, ijpj, :, :) = psgn * ptab2(ijt, ijpj - 2, :, :)
          END DO
          !$ACC END KERNELS
          IF (nimpp >= jpiglo / 2 + 1) THEN
            startloop = 1
          ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
            startloop = jpiglo / 2 + 1 - nimpp + 1
          ELSE
            startloop = nlci + 1
          END IF
          IF (startloop <= nlci) THEN
            !$ACC KERNELS
            DO ji = startloop, nlci
              ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 3
              ptab(ji, ijpjm1, :, :) = psgn * ptab2(ijt, ijpjm1, :, :)
            END DO
            !$ACC END KERNELS
          END IF
        CASE ('F')
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          !$ACC KERNELS
          DO ji = 1, endloop
            iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
            ptab(ji, ijpj, :, :) = psgn * ptab2(iju, ijpj - 2, :, :)
          END DO
          !$ACC END KERNELS
          IF ((nimpp + nlci - 1) .EQ. jpiglo) THEN
            !$ACC KERNELS
            ptab(nlci, ijpj, :, :) = psgn * ptab2(1, ijpj - 2, :, :)
            !$ACC END KERNELS
          END IF
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          IF (nimpp >= jpiglo / 2 + 1) THEN
            startloop = 1
          ELSE IF (nimpp + nlci - 1 >= jpiglo / 2 + 1 .AND. nimpp < jpiglo / 2 + 1) THEN
            startloop = jpiglo / 2 + 1 - nimpp + 1
          ELSE
            startloop = endloop + 1
          END IF
          IF (startloop <= endloop) THEN
            !$ACC KERNELS
            DO ji = startloop, endloop
              iju = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 2
              ptab(ji, ijpjm1, :, :) = psgn * ptab2(iju, ijpjm1, :, :)
            END DO
            !$ACC END KERNELS
          END IF
        CASE ('I')
          IF (nimpp /= 1) THEN
            startloop = 1
          ELSE
            startloop = 2
          END IF
          IF (nimpp + nlci - 1 /= jpiglo) THEN
            endloop = nlci
          ELSE
            endloop = nlci - 1
          END IF
          !$ACC KERNELS
          DO ji = startloop, endloop
            ijt = jpiglo - ji - nimpp - nfiimpp(isendto(1), jpnj) + 4
            ptab(ji, ijpj, :, :) = 0.5 * (ptab(ji, ijpjm1, :, :) + psgn * ptab2(ijt, ijpjm1, :, :))
          END DO
          !$ACC END KERNELS
        END SELECT
      CASE DEFAULT
        SELECT CASE (cd_nat)
        CASE ('T', 'U', 'V', 'W')
          !$ACC KERNELS
          ptab(:, 1, :, :) = 0._wp
          ptab(:, ijpj, :, :) = 0._wp
          !$ACC END KERNELS
        CASE ('F')
          !$ACC KERNELS
          ptab(:, ijpj, :, :) = 0._wp
          !$ACC END KERNELS
        CASE ('I')
          !$ACC KERNELS
          ptab(:, 1, :, :) = 0._wp
          ptab(:, ijpj, :, :) = 0._wp
          !$ACC END KERNELS
        END SELECT
      END SELECT
    END DO
  END SUBROUTINE lbc_nfd_nogather_4d
END MODULE lbcnfd