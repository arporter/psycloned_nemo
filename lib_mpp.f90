MODULE lib_mpp
  USE dom_oce
  USE lbcnfd
  USE in_out_manager
  IMPLICIT NONE
  PRIVATE
  INTERFACE mpp_nfd
    MODULE PROCEDURE mpp_nfd_2d, mpp_nfd_3d, mpp_nfd_4d
    MODULE PROCEDURE mpp_nfd_2d_ptr, mpp_nfd_3d_ptr, mpp_nfd_4d_ptr
  END INTERFACE
  PUBLIC :: mpp_lnk_2d, mpp_lnk_3d, mpp_lnk_4d
  PUBLIC :: mpp_lnk_2d_ptr, mpp_lnk_3d_ptr, mpp_lnk_4d_ptr
  PUBLIC :: mpp_nfd_2d, mpp_nfd_3d, mpp_nfd_4d
  PUBLIC :: mpp_nfd_2d_ptr, mpp_nfd_3d_ptr, mpp_nfd_4d_ptr
  PUBLIC :: ctl_stop, ctl_warn, get_unit, ctl_opn, ctl_nam
  PUBLIC :: mynode, mppstop, mppsync, mpp_comm_free
  PUBLIC :: mpp_ini_north
  PUBLIC :: mpp_lnk_2d_icb
  PUBLIC :: mpp_lbc_north_icb
  PUBLIC :: mpp_min, mpp_max, mpp_sum, mpp_minloc, mpp_maxloc
  PUBLIC :: mpp_max_multiple
  PUBLIC :: mppscatter, mppgather
  PUBLIC :: mpp_ini_ice, mpp_ini_znl
  PUBLIC :: mppsize
  PUBLIC :: mppsend, mpprecv
  PUBLIC :: mpp_lnk_bdy_2d, mpp_lnk_bdy_3d, mpp_lnk_bdy_4d
  PUBLIC :: mpprank
  INTERFACE mpp_min
    MODULE PROCEDURE mppmin_a_int, mppmin_int, mppmin_a_real, mppmin_real
  END INTERFACE
  INTERFACE mpp_max
    MODULE PROCEDURE mppmax_a_int, mppmax_int, mppmax_a_real, mppmax_real
  END INTERFACE
  INTERFACE mpp_sum
    MODULE PROCEDURE mppsum_a_int, mppsum_int, mppsum_a_real, mppsum_real, mppsum_realdd, mppsum_a_realdd
  END INTERFACE
  INTERFACE mpp_minloc
    MODULE PROCEDURE mpp_minloc2d, mpp_minloc3d
  END INTERFACE
  INTERFACE mpp_maxloc
    MODULE PROCEDURE mpp_maxloc2d, mpp_maxloc3d
  END INTERFACE
  INTERFACE mpp_max_multiple
    MODULE PROCEDURE mppmax_real_multiple
  END INTERFACE
  INCLUDE 'mpif.h'
  LOGICAL, PUBLIC, PARAMETER :: lk_mpp = .TRUE.
  INTEGER, PARAMETER :: nprocmax = 2 ** 10
  INTEGER :: mppsize
  INTEGER :: mpprank
  INTEGER, PUBLIC :: mpi_comm_oce
  INTEGER :: MPI_SUMDD
  INTEGER, PUBLIC :: ncomm_ice
  INTEGER :: ngrp_iworld
  INTEGER :: ngrp_ice
  INTEGER :: ndim_rank_ice
  INTEGER :: n_ice_root
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nrank_ice
  INTEGER, PUBLIC :: ncomm_znl
  LOGICAL, PUBLIC :: l_znl_root
  INTEGER :: ngrp_znl
  INTEGER :: ndim_rank_znl
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nrank_znl
  INTEGER, PUBLIC :: ngrp_world
  INTEGER, PUBLIC :: ngrp_opa
  INTEGER, PUBLIC :: ngrp_north
  INTEGER, PUBLIC :: ncomm_north
  INTEGER, PUBLIC :: ndim_rank_north
  INTEGER, PUBLIC :: njmppmax
  INTEGER, PUBLIC :: north_root
  INTEGER, PUBLIC, DIMENSION(:), ALLOCATABLE, SAVE :: nrank_north
  CHARACTER(LEN = 1), PUBLIC :: cn_mpi_send
  LOGICAL, PUBLIC :: l_isend = .FALSE.
  INTEGER, PUBLIC :: nn_buffer
  REAL(KIND = wp), DIMENSION(:), ALLOCATABLE, SAVE :: tampon
  LOGICAL, PUBLIC :: ln_nnogather
  LOGICAL, PUBLIC :: l_north_nogather = .FALSE.
  CONTAINS
  FUNCTION mynode(ldtxt, ldname, kumnam_ref, kumnam_cfg, kumond, kstop, localComm)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    USE openacc
    CHARACTER(LEN = *), DIMENSION(:), INTENT(OUT) :: ldtxt
    CHARACTER(LEN = *), INTENT(IN) :: ldname
    INTEGER, INTENT(IN) :: kumnam_ref
    INTEGER, INTENT(IN) :: kumnam_cfg
    INTEGER, INTENT(INOUT) :: kumond
    INTEGER, INTENT(INOUT) :: kstop
    INTEGER, OPTIONAL, INTENT(IN) :: localComm
    INTEGER :: mynode, ierr, code, ji, ii, ios
    LOGICAL :: mpi_was_called
    NAMELIST /nammpp/ cn_mpi_send, nn_buffer, jpni, jpnj, jpnij, ln_nnogather
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mynode', 'r0', 0, 0)
    ii = 1
    WRITE(ldtxt(ii), FMT = *)
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) 'mynode : mpi initialisation'
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) '~~~~~~ '
    ii = ii + 1
    REWIND(UNIT = kumnam_ref)
    READ(kumnam_ref, nammpp, IOSTAT = ios, ERR = 901)
901 IF (ios /= 0) CALL ctl_nam(ios, 'nammpp in reference namelist', lwp)
    REWIND(UNIT = kumnam_cfg)
    READ(kumnam_cfg, nammpp, IOSTAT = ios, ERR = 902)
902 IF (ios > 0) CALL ctl_nam(ios, 'nammpp in configuration namelist', lwp)
    WRITE(ldtxt(ii), FMT = *) '   Namelist nammpp'
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) '      mpi send type          cn_mpi_send = ', cn_mpi_send
    ii = ii + 1
    WRITE(ldtxt(ii), FMT = *) '      size exported buffer   nn_buffer   = ', nn_buffer, ' bytes'
    ii = ii + 1
    IF (jpnij < 1) THEN
      jpnij = jpni * jpnj
    END IF
    IF (jpni < 1 .OR. jpnj < 1) THEN
      WRITE(ldtxt(ii), FMT = *) '      jpni, jpnj and jpnij will be calculated automatically'
      ii = ii + 1
    ELSE
      WRITE(ldtxt(ii), FMT = *) '      processor grid extent in i         jpni = ', jpni
      ii = ii + 1
      WRITE(ldtxt(ii), FMT = *) '      processor grid extent in j         jpnj = ', jpnj
      ii = ii + 1
      WRITE(ldtxt(ii), FMT = *) '      number of local domains           jpnij = ', jpnij
      ii = ii + 1
    END IF
    WRITE(ldtxt(ii), FMT = *) '      avoid use of mpi_allgather at the north fold  ln_nnogather = ', ln_nnogather
    ii = ii + 1
    CALL mpi_initialized(mpi_was_called, code)
    IF (code /= MPI_SUCCESS) THEN
      DO ji = 1, SIZE(ldtxt)
        IF (TRIM(ldtxt(ji)) /= '') WRITE(*, FMT = *) ldtxt(ji)
      END DO
      WRITE(*, cform_err)
      WRITE(*, FMT = *) 'lib_mpp: Error in routine mpi_initialized'
      CALL mpi_abort(mpi_comm_world, code, ierr)
    END IF
    IF (mpi_was_called) THEN
      SELECT CASE (cn_mpi_send)
      CASE ('S')
        WRITE(ldtxt(ii), FMT = *) '           Standard blocking mpi send (send)'
        ii = ii + 1
      CASE ('B')
        WRITE(ldtxt(ii), FMT = *) '           Buffer blocking mpi send (bsend)'
        ii = ii + 1
        IF (Agrif_Root()) CALL mpi_init_oce(ldtxt, ii, ierr)
      CASE ('I')
        WRITE(ldtxt(ii), FMT = *) '           Immediate non-blocking send (isend)'
        ii = ii + 1
        l_isend = .TRUE.
      CASE DEFAULT
        WRITE(ldtxt(ii), cform_err)
        ii = ii + 1
        WRITE(ldtxt(ii), FMT = *) '           bad value for cn_mpi_send = ', cn_mpi_send
        ii = ii + 1
        kstop = kstop + 1
      END SELECT
    ELSE IF (PRESENT(localComm) .AND. .NOT. mpi_was_called) THEN
      WRITE(ldtxt(ii), FMT = *) ' lib_mpp: You cannot provide a local communicator '
      ii = ii + 1
      WRITE(ldtxt(ii), FMT = *) '          without calling MPI_Init before ! '
      ii = ii + 1
      kstop = kstop + 1
    ELSE
      SELECT CASE (cn_mpi_send)
      CASE ('S')
        WRITE(ldtxt(ii), FMT = *) '           Standard blocking mpi send (send)'
        ii = ii + 1
        CALL mpi_init(ierr)
      CASE ('B')
        WRITE(ldtxt(ii), FMT = *) '           Buffer blocking mpi send (bsend)'
        ii = ii + 1
        IF (Agrif_Root()) CALL mpi_init_oce(ldtxt, ii, ierr)
      CASE ('I')
        WRITE(ldtxt(ii), FMT = *) '           Immediate non-blocking send (isend)'
        ii = ii + 1
        l_isend = .TRUE.
        CALL mpi_init(ierr)
      CASE DEFAULT
        WRITE(ldtxt(ii), cform_err)
        ii = ii + 1
        WRITE(ldtxt(ii), FMT = *) '           bad value for cn_mpi_send = ', cn_mpi_send
        ii = ii + 1
        kstop = kstop + 1
      END SELECT
    END IF
    IF (PRESENT(localComm)) THEN
      IF (Agrif_Root()) THEN
        mpi_comm_oce = localComm
      END IF
    ELSE
      CALL mpi_comm_dup(mpi_comm_world, mpi_comm_oce, code)
      IF (code /= MPI_SUCCESS) THEN
        DO ji = 1, SIZE(ldtxt)
          IF (TRIM(ldtxt(ji)) /= '') WRITE(*, FMT = *) ldtxt(ji)
        END DO
        WRITE(*, cform_err)
        WRITE(*, FMT = *) ' lib_mpp: Error in routine mpi_comm_dup'
        CALL mpi_abort(mpi_comm_world, code, ierr)
      END IF
    END IF
    CALL mpi_comm_rank(mpi_comm_oce, mpprank, ierr)
    CALL mpi_comm_size(mpi_comm_oce, mppsize, ierr)

    call acc_set_device_num( mpprank, acc_device_nvidia )
    write(*,*) "ARPDBG: MPI rank ", mpprank, " using GPU device ", mpprank

    mynode = mpprank
    IF (mynode == 0) THEN
      CALL ctl_opn(kumond, TRIM(ldname), 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', - 1, 6, .FALSE., 1)
      WRITE(kumond, nammpp)
    END IF
    CALL MPI_OP_CREATE(DDPDD_MPI, .TRUE., MPI_SUMDD, ierr)
    CALL profile_psy_data0 % PostEnd
  END FUNCTION mynode
  SUBROUTINE mpp_lnk_2d(ptab, cd_nat, psgn, cd_mpp, pval)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = 1
    ipl = 1
    ipf = 1
    ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2), zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), &
&zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    IF (PRESENT(cd_mpp)) THEN
      !$ACC KERNELS
      DO jf = 1, ipf
        DO jl = 1, ipl
          !DO jk = 1, ipk
            DO jj = nlcj + 1, jpj
              ptab(nldi : nlei, jj) = ptab(nldi : nlei, nlej)
              ptab(1 : nldi - 1, jj) = ptab(nldi, nlej)
              ptab(nlei + 1 : nlci, jj) = ptab(nlei, nlej)
            END DO
            DO ji = nlci + 1, jpi
              ptab(ji, nldj : nlej) = ptab(nlei, nldj : nlej)
              ptab(ji, 1 : nldj - 1) = ptab(nlei, nldj)
              ptab(ji, nlej + 1 : jpj) = ptab(nlei, nlej)
            END DO
          !END DO
        END DO
     END DO
     !$ACC END KERNELS
    ELSE
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :) = ptab(jpim1, :)
          ptab(jpi, :) = ptab(2, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :) = zland
          ptab(nlci - nn_hls + 1 : jpi, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1) = ptab(:, jpjm1)
          ptab(:, jpj) = ptab(:, 2)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls) = zland
          ptab(:, nlcj - nn_hls + 1 : jpj) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
    SELECT CASE (nbondi)
    CASE (- 1, 0, 1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpj * ipk * ipl * ipf
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    iihom = nlci - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(iihom + jh, :) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(iihom + jh, :) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :) = zt3we(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (nbondj /= 2) THEN
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END IF
    imigr = nn_hls * jpi * ipk * ipl * ipf
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    ijhom = nlcj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, ijhom + jh) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(:, ijhom + jh) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh) = zt3sn(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (npolj /= 0 .AND. .NOT. PRESENT(cd_mpp)) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat, psgn)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat, psgn)
      END SELECT
    END IF
    DEALLOCATE(zt3ns, zt3sn, zt3ew, zt3we)
  END SUBROUTINE mpp_lnk_2d
  SUBROUTINE mpp_lnk_2d_ptr(ptab, cd_nat, psgn, kfld, cd_mpp, pval)
    INTEGER, INTENT(IN) :: kfld
    TYPE(PTR_2D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = 1
    ipl = 1
    ipf = kfld
    ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2), zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), &
&zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    IF (PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj + 1, jpj
              ptab(jf) % pt2d(nldi : nlei, jj) = ptab(jf) % pt2d(nldi : nlei, nlej)
              ptab(jf) % pt2d(1 : nldi - 1, jj) = ptab(jf) % pt2d(nldi, nlej)
              ptab(jf) % pt2d(nlei + 1 : nlci, jj) = ptab(jf) % pt2d(nlei, nlej)
            END DO
            DO ji = nlci + 1, jpi
              ptab(jf) % pt2d(ji, nldj : nlej) = ptab(jf) % pt2d(nlei, nldj : nlej)
              ptab(jf) % pt2d(ji, 1 : nldj - 1) = ptab(jf) % pt2d(nlei, nldj)
              ptab(jf) % pt2d(ji, nlej + 1 : jpj) = ptab(jf) % pt2d(nlei, nlej)
            END DO
          END DO
        END DO
      END DO
    ELSE
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          ptab(jf) % pt2d(1, :) = ptab(jf) % pt2d(jpim1, :)
          ptab(jf) % pt2d(jpi, :) = ptab(jf) % pt2d(2, :)
        ELSE
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt2d(1 : nn_hls, :) = zland
          ptab(jf) % pt2d(nlci - nn_hls + 1 : jpi, :) = zland
        END IF
        IF (l_Jperio) THEN
          ptab(jf) % pt2d(:, 1) = ptab(jf) % pt2d(:, jpjm1)
          ptab(jf) % pt2d(:, jpj) = ptab(jf) % pt2d(:, 2)
        ELSE
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt2d(:, 1 : nn_hls) = zland
          ptab(jf) % pt2d(:, nlcj - nn_hls + 1 : jpj) = zland
        END IF
      END DO
    END IF
    SELECT CASE (nbondi)
    CASE (- 1, 0, 1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(jf) % pt2d(nn_hls + jh, :)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(jf) % pt2d(iihom + jh, :)
            END DO
          END DO
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpj * ipk * ipl * ipf
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    iihom = nlci - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt2d(iihom + jh, :) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt2d(jh, :) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(jf) % pt2d(iihom + jh, :) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt2d(jh, :) = zt3we(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    END SELECT
    IF (nbondj /= 2) THEN
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(jf) % pt2d(:, ijhom + jh)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(jf) % pt2d(:, nn_hls + jh)
            END DO
          END DO
        END DO
      END DO
    END IF
    imigr = nn_hls * jpi * ipk * ipl * ipf
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    ijhom = nlcj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt2d(:, ijhom + jh) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt2d(:, jh) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(jf) % pt2d(:, ijhom + jh) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt2d(:, jh) = zt3sn(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    END SELECT
    IF (npolj /= 0 .AND. .NOT. PRESENT(cd_mpp)) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat(:), psgn(:), ipf)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat(:), psgn(:), ipf)
      END SELECT
    END IF
    DEALLOCATE(zt3ns, zt3sn, zt3ew, zt3we)
  END SUBROUTINE mpp_lnk_2d_ptr
  SUBROUTINE mpp_lnk_3d(ptab, cd_nat, psgn, cd_mpp, pval)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = SIZE(ptab, 3)
    ipl = 1
    ipf = 1
    ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2), zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), &
&zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    IF (PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            !$ACC KERNELS
            DO jj = nlcj + 1, jpj
              ptab(nldi : nlei, jj, jk) = ptab(nldi : nlei, nlej, jk)
              ptab(1 : nldi - 1, jj, jk) = ptab(nldi, nlej, jk)
              ptab(nlei + 1 : nlci, jj, jk) = ptab(nlei, nlej, jk)
            END DO
            DO ji = nlci + 1, jpi
              ptab(ji, nldj : nlej, jk) = ptab(nlei, nldj : nlej, jk)
              ptab(ji, 1 : nldj - 1, jk) = ptab(nlei, nldj, jk)
              ptab(ji, nlej + 1 : jpj, jk) = ptab(nlei, nlej, jk)
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
    ELSE
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :, :) = ptab(jpim1, :, :)
          ptab(jpi, :, :) = ptab(2, :, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :, :) = zland
          ptab(nlci - nn_hls + 1 : jpi, :, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1, :) = ptab(:, jpjm1, :)
          ptab(:, jpj, :) = ptab(:, 2, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls, :) = zland
          ptab(:, nlcj - nn_hls + 1 : jpj, :) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
    SELECT CASE (nbondi)
    CASE (- 1, 0, 1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :, jk)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpj * ipk * ipl * ipf
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    iihom = nlci - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(iihom + jh, :, jk) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(iihom + jh, :, jk) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk) = zt3we(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (nbondj /= 2) THEN
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh, jk)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END IF
    imigr = nn_hls * jpi * ipk * ipl * ipf
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    ijhom = nlcj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, ijhom + jh, jk) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(:, ijhom + jh, jk) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk) = zt3sn(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (npolj /= 0 .AND. .NOT. PRESENT(cd_mpp)) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat, psgn)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat, psgn)
      END SELECT
    END IF
    DEALLOCATE(zt3ns, zt3sn, zt3ew, zt3we)
  END SUBROUTINE mpp_lnk_3d
  SUBROUTINE mpp_lnk_3d_ptr(ptab, cd_nat, psgn, kfld, cd_mpp, pval)
    INTEGER, INTENT(IN) :: kfld
    TYPE(PTR_3D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = SIZE(ptab(1) % pt3d, 3)
    ipl = 1
    ipf = kfld
    ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2), zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), &
&zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    IF (PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj + 1, jpj
              ptab(jf) % pt3d(nldi : nlei, jj, jk) = ptab(jf) % pt3d(nldi : nlei, nlej, jk)
              ptab(jf) % pt3d(1 : nldi - 1, jj, jk) = ptab(jf) % pt3d(nldi, nlej, jk)
              ptab(jf) % pt3d(nlei + 1 : nlci, jj, jk) = ptab(jf) % pt3d(nlei, nlej, jk)
            END DO
            DO ji = nlci + 1, jpi
              ptab(jf) % pt3d(ji, nldj : nlej, jk) = ptab(jf) % pt3d(nlei, nldj : nlej, jk)
              ptab(jf) % pt3d(ji, 1 : nldj - 1, jk) = ptab(jf) % pt3d(nlei, nldj, jk)
              ptab(jf) % pt3d(ji, nlej + 1 : jpj, jk) = ptab(jf) % pt3d(nlei, nlej, jk)
            END DO
          END DO
        END DO
      END DO
    ELSE
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          ptab(jf) % pt3d(1, :, :) = ptab(jf) % pt3d(jpim1, :, :)
          ptab(jf) % pt3d(jpi, :, :) = ptab(jf) % pt3d(2, :, :)
        ELSE
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt3d(1 : nn_hls, :, :) = zland
          ptab(jf) % pt3d(nlci - nn_hls + 1 : jpi, :, :) = zland
        END IF
        IF (l_Jperio) THEN
          ptab(jf) % pt3d(:, 1, :) = ptab(jf) % pt3d(:, jpjm1, :)
          ptab(jf) % pt3d(:, jpj, :) = ptab(jf) % pt3d(:, 2, :)
        ELSE
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt3d(:, 1 : nn_hls, :) = zland
          ptab(jf) % pt3d(:, nlcj - nn_hls + 1 : jpj, :) = zland
        END IF
      END DO
    END IF
    SELECT CASE (nbondi)
    CASE (- 1, 0, 1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(jf) % pt3d(nn_hls + jh, :, jk)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(jf) % pt3d(iihom + jh, :, jk)
            END DO
          END DO
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpj * ipk * ipl * ipf
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    iihom = nlci - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt3d(iihom + jh, :, jk) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt3d(jh, :, jk) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(jf) % pt3d(iihom + jh, :, jk) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt3d(jh, :, jk) = zt3we(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    END SELECT
    IF (nbondj /= 2) THEN
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(jf) % pt3d(:, ijhom + jh, jk)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(jf) % pt3d(:, nn_hls + jh, jk)
            END DO
          END DO
        END DO
      END DO
    END IF
    imigr = nn_hls * jpi * ipk * ipl * ipf
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    ijhom = nlcj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt3d(:, ijhom + jh, jk) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt3d(:, jh, jk) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(jf) % pt3d(:, ijhom + jh, jk) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt3d(:, jh, jk) = zt3sn(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    END SELECT
    IF (npolj /= 0 .AND. .NOT. PRESENT(cd_mpp)) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat(:), psgn(:), ipf)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat(:), psgn(:), ipf)
      END SELECT
    END IF
    DEALLOCATE(zt3ns, zt3sn, zt3ew, zt3we)
  END SUBROUTINE mpp_lnk_3d_ptr
  SUBROUTINE mpp_lnk_4d(ptab, cd_nat, psgn, cd_mpp, pval)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = SIZE(ptab, 3)
    ipl = SIZE(ptab, 4)
    ipf = 1
    ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2), zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), &
&zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    IF (PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            !$ACC KERNELS
            DO jj = nlcj + 1, jpj
              ptab(nldi : nlei, jj, jk, jl) = ptab(nldi : nlei, nlej, jk, jl)
              ptab(1 : nldi - 1, jj, jk, jl) = ptab(nldi, nlej, jk, jl)
              ptab(nlei + 1 : nlci, jj, jk, jl) = ptab(nlei, nlej, jk, jl)
            END DO
            DO ji = nlci + 1, jpi
              ptab(ji, nldj : nlej, jk, jl) = ptab(nlei, nldj : nlej, jk, jl)
              ptab(ji, 1 : nldj - 1, jk, jl) = ptab(nlei, nldj, jk, jl)
              ptab(ji, nlej + 1 : jpj, jk, jl) = ptab(nlei, nlej, jk, jl)
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
    ELSE
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :, :, :) = ptab(jpim1, :, :, :)
          ptab(jpi, :, :, :) = ptab(2, :, :, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :, :, :) = zland
          ptab(nlci - nn_hls + 1 : jpi, :, :, :) = zland
          !$ACC END KERNELS
        END IF
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1, :, :) = ptab(:, jpjm1, :, :)
          ptab(:, jpj, :, :) = ptab(:, 2, :, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls, :, :) = zland
          ptab(:, nlcj - nn_hls + 1 : jpj, :, :) = zland
          !$ACC END KERNELS
        END IF
      END DO
    END IF
    SELECT CASE (nbondi)
    CASE (- 1, 0, 1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :, jk, jl)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :, jk, jl)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpj * ipk * ipl * ipf
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    iihom = nlci - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(iihom + jh, :, jk, jl) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk, jl) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(iihom + jh, :, jk, jl) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk, jl) = zt3we(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (nbondj /= 2) THEN
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh, jk, jl)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh, jk, jl)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END IF
    imigr = nn_hls * jpi * ipk * ipl * ipf
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    ijhom = nlcj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, ijhom + jh, jk, jl) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk, jl) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(:, ijhom + jh, jk, jl) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk, jl) = zt3sn(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
    END SELECT
    IF (npolj /= 0 .AND. .NOT. PRESENT(cd_mpp)) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat, psgn)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat, psgn)
      END SELECT
    END IF
    DEALLOCATE(zt3ns, zt3sn, zt3ew, zt3we)
  END SUBROUTINE mpp_lnk_4d
  SUBROUTINE mpp_lnk_4d_ptr(ptab, cd_nat, psgn, kfld, cd_mpp, pval)
    INTEGER, INTENT(IN) :: kfld
    TYPE(PTR_4D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    CHARACTER(LEN = 3), OPTIONAL, INTENT(IN) :: cd_mpp
    REAL(KIND = wp), OPTIONAL, INTENT(IN) :: pval
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = SIZE(ptab(1) % pt4d, 3)
    ipl = SIZE(ptab(1) % pt4d, 4)
    ipf = kfld
    ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2), zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), &
&zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    IF (PRESENT(pval)) THEN
      zland = pval
    ELSE
      zland = 0._wp
    END IF
    IF (PRESENT(cd_mpp)) THEN
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj + 1, jpj
              ptab(jf) % pt4d(nldi : nlei, jj, jk, jl) = ptab(jf) % pt4d(nldi : nlei, nlej, jk, jl)
              ptab(jf) % pt4d(1 : nldi - 1, jj, jk, jl) = ptab(jf) % pt4d(nldi, nlej, jk, jl)
              ptab(jf) % pt4d(nlei + 1 : nlci, jj, jk, jl) = ptab(jf) % pt4d(nlei, nlej, jk, jl)
            END DO
            DO ji = nlci + 1, jpi
              ptab(jf) % pt4d(ji, nldj : nlej, jk, jl) = ptab(jf) % pt4d(nlei, nldj : nlej, jk, jl)
              ptab(jf) % pt4d(ji, 1 : nldj - 1, jk, jl) = ptab(jf) % pt4d(nlei, nldj, jk, jl)
              ptab(jf) % pt4d(ji, nlej + 1 : jpj, jk, jl) = ptab(jf) % pt4d(nlei, nlej, jk, jl)
            END DO
          END DO
        END DO
      END DO
    ELSE
      DO jf = 1, ipf
        IF (l_Iperio) THEN
          ptab(jf) % pt4d(1, :, :, :) = ptab(jf) % pt4d(jpim1, :, :, :)
          ptab(jf) % pt4d(jpi, :, :, :) = ptab(jf) % pt4d(2, :, :, :)
        ELSE
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt4d(1 : nn_hls, :, :, :) = zland
          ptab(jf) % pt4d(nlci - nn_hls + 1 : jpi, :, :, :) = zland
        END IF
        IF (l_Jperio) THEN
          ptab(jf) % pt4d(:, 1, :, :) = ptab(jf) % pt4d(:, jpjm1, :, :)
          ptab(jf) % pt4d(:, jpj, :, :) = ptab(jf) % pt4d(:, 2, :, :)
        ELSE
          IF (.NOT. cd_nat(jf) == 'F') ptab(jf) % pt4d(:, 1 : nn_hls, :, :) = zland
          ptab(jf) % pt4d(:, nlcj - nn_hls + 1 : jpj, :, :) = zland
        END IF
      END DO
    END IF
    SELECT CASE (nbondi)
    CASE (- 1, 0, 1)
      iihom = nlci - nreci
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(jf) % pt4d(nn_hls + jh, :, jk, jl)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(jf) % pt4d(iihom + jh, :, jk, jl)
            END DO
          END DO
        END DO
      END DO
    END SELECT
    imigr = nn_hls * jpj * ipk * ipl * ipf
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    iihom = nlci - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt4d(iihom + jh, :, jk, jl) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt4d(jh, :, jk, jl) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(jf) % pt4d(iihom + jh, :, jk, jl) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt4d(jh, :, jk, jl) = zt3we(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    END SELECT
    IF (nbondj /= 2) THEN
      ijhom = nlcj - nrecj
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(jf) % pt4d(:, ijhom + jh, jk, jl)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(jf) % pt4d(:, nn_hls + jh, jk, jl)
            END DO
          END DO
        END DO
      END DO
    END IF
    imigr = nn_hls * jpi * ipk * ipl * ipf
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    ijhom = nlcj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt4d(:, ijhom + jh, jk, jl) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (0)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt4d(:, jh, jk, jl) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(jf) % pt4d(:, ijhom + jh, jk, jl) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    CASE (1)
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jf) % pt4d(:, jh, jk, jl) = zt3sn(:, jh, jk, jl, jf, 2)
            END DO
          END DO
        END DO
      END DO
    END SELECT
    IF (npolj /= 0 .AND. .NOT. PRESENT(cd_mpp)) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat(:), psgn(:), ipf)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat(:), psgn(:), ipf)
      END SELECT
    END IF
    DEALLOCATE(zt3ns, zt3sn, zt3ew, zt3we)
  END SUBROUTINE mpp_lnk_4d_ptr
  SUBROUTINE mpp_nfd_2d(ptab, cd_nat, psgn, kfld)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf, jr
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ierr, itaille, ilci, ildi, ilei, iilb
    INTEGER :: ij, iproc
    INTEGER, DIMENSION(jpmaxngh) :: ml_req_nf
    INTEGER :: ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: ztab, ztabl, ztabr
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: znorthloc, zfoldwk
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: znorthgloio
    ipk = 1
    ipl = 1
    ipf = 1
    ipj = 4
    ALLOCATE(znorthloc(jpimax, 4, ipk, ipl, ipf))
    !$ACC KERNELS
    znorthloc(:, :, :, :, :) = 0._wp
    !$ACC END KERNELS
    DO jf = 1, ipf
      DO jl = 1, ipl
        !$ACC KERNELS
        DO jk = 1, ipk
          DO jj = nlcj - ipj + 1, nlcj
            ij = jj - nlcj + ipj
            znorthloc(1 : jpi, ij, jk, jl, jf) = ptab(1 : jpi, jj)
          END DO
        END DO
        !$ACC END KERNELS
      END DO
    END DO
    itaille = jpimax * ipj * ipk * ipl * ipf
    IF (l_north_nogather) THEN
      ALLOCATE(zfoldwk(jpimax, 4, ipk, ipl, ipf))
      ALLOCATE(ztabl(jpimax, 4, ipk, ipl, ipf), ztabr(jpimax * jpmaxngh, 4, ipk, ipl, ipf))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) THEN
        ztabr(:, :, :, :, :) = 0._wp
        ztabl(:, :, :, :, :) = 0._wp
      END IF
      !$ACC END KERNELS
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = nfsloop, nfeloop
                ztabl(ji, ij, jk, jl, jf) = ptab(ji, jj)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DO jr = 1, nsndto
        IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
          CALL mppsend(5, znorthloc, itaille, nfipproc(isendto(jr), jpnj), ml_req_nf(jr))
        END IF
      END DO
      DO jr = 1, nsndto
        iproc = nfipproc(isendto(jr), jpnj)
        IF (iproc /= - 1) THEN
          iilb = nimppt(iproc + 1)
          ilci = nlcit(iproc + 1)
          ildi = nldit(iproc + 1)
          ilei = nleit(iproc + 1)
          IF (iilb == 1) ildi = 1
          IF (iilb + ilci - 1 == jpiglo) ilei = ilci
          iilb = nfiimpp(isendto(jr), jpnj) - nfiimpp(isendto(1), jpnj)
        END IF
        IF (iproc /= narea - 1 .AND. iproc /= - 1) THEN
          CALL mpprecv(5, zfoldwk, itaille, iproc)
          DO jf = 1, ipf
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                !$ACC LOOP INDEPENDENT COLLAPSE(2)
                DO jj = 1, ipj
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = zfoldwk(ji, jj, jk, jl, jf)
                  END DO
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END DO
        ELSE IF (iproc == narea - 1) THEN
          DO jf = 1, ipf
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                !$ACC LOOP INDEPENDENT COLLAPSE(2)
                DO jj = 1, ipj
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = ptab(ji, nlcj - ipj + jj)
                  END DO
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END DO
        END IF
      END DO
      IF (l_isend) THEN
        DO jr = 1, nsndto
          IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
            CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
          END IF
        END DO
      END IF
      DO jf = 1, ipf
        CALL lbc_nfd_nogather(ztabl(:, :, :, :, jf), ztabr(:, :, :, :, jf), cd_nat, psgn)
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(ji, jj) = ztabl(ji, ij, jk, jl, jf)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DEALLOCATE(zfoldwk)
      DEALLOCATE(ztabl, ztabr)
    ELSE
      ALLOCATE(ztab(jpiglo, 4, ipk, ipl, ipf))
      ALLOCATE(znorthgloio(jpimax, 4, ipk, ipl, ipf, jpni))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztab(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      CALL MPI_ALLGATHER(znorthloc, itaille, MPI_DOUBLE_PRECISION, znorthgloio, itaille, MPI_DOUBLE_PRECISION, ncomm_north, ierr)
      DO jr = 1, ndim_rank_north
        iproc = nrank_north(jr) + 1
        iilb = nimppt(iproc)
        ilci = nlcit(iproc)
        ildi = nldit(iproc)
        ilei = nleit(iproc)
        IF (iilb == 1) ildi = 1
        IF (iilb + ilci - 1 == jpiglo) ilei = ilci
        DO jf = 1, ipf
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              !$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO jj = 1, ipj
                DO ji = ildi, ilei
                  ztab(ji + iilb - 1, jj, jk, jl, jf) = znorthgloio(ji, jj, jk, jl, jf, jr)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      DO jf = 1, ipf
        CALL lbc_nfd(ztab(:, :, :, :, jf), cd_nat, psgn)
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(ji, jj) = ztab(ji + nimpp - 1, ij, jk, jl, jf)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DEALLOCATE(ztab)
      DEALLOCATE(znorthgloio)
    END IF
    DEALLOCATE(znorthloc)
  END SUBROUTINE mpp_nfd_2d
  SUBROUTINE mpp_nfd_2d_ptr(ptab, cd_nat, psgn, kfld)
    TYPE(PTR_2D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf, jr
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ierr, itaille, ilci, ildi, ilei, iilb
    INTEGER :: ij, iproc
    INTEGER, DIMENSION(jpmaxngh) :: ml_req_nf
    INTEGER :: ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: ztab, ztabl, ztabr
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: znorthloc, zfoldwk
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: znorthgloio
    ipk = 1
    ipl = 1
    ipf = kfld
    ipj = 4
    ALLOCATE(znorthloc(jpimax, 4, ipk, ipl, ipf))
    !$ACC KERNELS
    znorthloc(:, :, :, :, :) = 0._wp
    !$ACC END KERNELS
    DO jf = 1, ipf
      DO jl = 1, ipl
        DO jk = 1, ipk
          DO jj = nlcj - ipj + 1, nlcj
            ij = jj - nlcj + ipj
            znorthloc(1 : jpi, ij, jk, jl, jf) = ptab(jf) % pt2d(1 : jpi, jj)
          END DO
        END DO
      END DO
    END DO
    itaille = jpimax * ipj * ipk * ipl * ipf
    IF (l_north_nogather) THEN
      ALLOCATE(zfoldwk(jpimax, 4, ipk, ipl, ipf))
      ALLOCATE(ztabl(jpimax, 4, ipk, ipl, ipf), ztabr(jpimax * jpmaxngh, 4, ipk, ipl, ipf))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) THEN
        ztabr(:, :, :, :, :) = 0._wp
        ztabl(:, :, :, :, :) = 0._wp
      END IF
      !$ACC END KERNELS
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = nfsloop, nfeloop
                ztabl(ji, ij, jk, jl, jf) = ptab(jf) % pt2d(ji, jj)
              END DO
            END DO
          END DO
        END DO
      END DO
      DO jr = 1, nsndto
        IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
          CALL mppsend(5, znorthloc, itaille, nfipproc(isendto(jr), jpnj), ml_req_nf(jr))
        END IF
      END DO
      DO jr = 1, nsndto
        iproc = nfipproc(isendto(jr), jpnj)
        IF (iproc /= - 1) THEN
          iilb = nimppt(iproc + 1)
          ilci = nlcit(iproc + 1)
          ildi = nldit(iproc + 1)
          ilei = nleit(iproc + 1)
          IF (iilb == 1) ildi = 1
          IF (iilb + ilci - 1 == jpiglo) ilei = ilci
          iilb = nfiimpp(isendto(jr), jpnj) - nfiimpp(isendto(1), jpnj)
        END IF
        IF (iproc /= narea - 1 .AND. iproc /= - 1) THEN
          CALL mpprecv(5, zfoldwk, itaille, iproc)
          DO jf = 1, ipf
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                !$ACC LOOP INDEPENDENT COLLAPSE(2)
                DO jj = 1, ipj
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = zfoldwk(ji, jj, jk, jl, jf)
                  END DO
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END DO
        ELSE IF (iproc == narea - 1) THEN
          DO jf = 1, ipf
            DO jl = 1, ipl
              DO jk = 1, ipk
                DO jj = 1, ipj
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = ptab(jf) % pt2d(ji, nlcj - ipj + jj)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
      END DO
      IF (l_isend) THEN
        DO jr = 1, nsndto
          IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
            CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
          END IF
        END DO
      END IF
      DO jf = 1, ipf
        CALL lbc_nfd_nogather(ztabl(:, :, :, :, jf), ztabr(:, :, :, :, jf), cd_nat(jf), psgn(jf))
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(jf) % pt2d(ji, jj) = ztabl(ji, ij, jk, jl, jf)
              END DO
            END DO
          END DO
        END DO
      END DO
      DEALLOCATE(zfoldwk)
      DEALLOCATE(ztabl, ztabr)
    ELSE
      ALLOCATE(ztab(jpiglo, 4, ipk, ipl, ipf))
      ALLOCATE(znorthgloio(jpimax, 4, ipk, ipl, ipf, jpni))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztab(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      CALL MPI_ALLGATHER(znorthloc, itaille, MPI_DOUBLE_PRECISION, znorthgloio, itaille, MPI_DOUBLE_PRECISION, ncomm_north, ierr)
      DO jr = 1, ndim_rank_north
        iproc = nrank_north(jr) + 1
        iilb = nimppt(iproc)
        ilci = nlcit(iproc)
        ildi = nldit(iproc)
        ilei = nleit(iproc)
        IF (iilb == 1) ildi = 1
        IF (iilb + ilci - 1 == jpiglo) ilei = ilci
        DO jf = 1, ipf
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              !$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO jj = 1, ipj
                DO ji = ildi, ilei
                  ztab(ji + iilb - 1, jj, jk, jl, jf) = znorthgloio(ji, jj, jk, jl, jf, jr)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      DO jf = 1, ipf
        CALL lbc_nfd(ztab(:, :, :, :, jf), cd_nat(jf), psgn(jf))
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(jf) % pt2d(ji, jj) = ztab(ji + nimpp - 1, ij, jk, jl, jf)
              END DO
            END DO
          END DO
        END DO
      END DO
      DEALLOCATE(ztab)
      DEALLOCATE(znorthgloio)
    END IF
    DEALLOCATE(znorthloc)
  END SUBROUTINE mpp_nfd_2d_ptr
  SUBROUTINE mpp_nfd_3d(ptab, cd_nat, psgn, kfld)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf, jr
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ierr, itaille, ilci, ildi, ilei, iilb
    INTEGER :: ij, iproc
    INTEGER, DIMENSION(jpmaxngh) :: ml_req_nf
    INTEGER :: ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: ztab, ztabl, ztabr
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: znorthloc, zfoldwk
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: znorthgloio
    ipk = SIZE(ptab, 3)
    ipl = 1
    ipf = 1
    ipj = 4
    ALLOCATE(znorthloc(jpimax, 4, ipk, ipl, ipf))
    !$ACC KERNELS
    znorthloc(:, :, :, :, :) = 0._wp
    !$ACC END KERNELS
    DO jf = 1, ipf
      DO jl = 1, ipl
        !$ACC KERNELS
        DO jk = 1, ipk
          DO jj = nlcj - ipj + 1, nlcj
            ij = jj - nlcj + ipj
            znorthloc(1 : jpi, ij, jk, jl, jf) = ptab(1 : jpi, jj, jk)
          END DO
        END DO
        !$ACC END KERNELS
      END DO
    END DO
    itaille = jpimax * ipj * ipk * ipl * ipf
    IF (l_north_nogather) THEN
      ALLOCATE(zfoldwk(jpimax, 4, ipk, ipl, ipf))
      ALLOCATE(ztabl(jpimax, 4, ipk, ipl, ipf), ztabr(jpimax * jpmaxngh, 4, ipk, ipl, ipf))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) THEN
        ztabr(:, :, :, :, :) = 0._wp
        ztabl(:, :, :, :, :) = 0._wp
      END IF
      !$ACC END KERNELS
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = nfsloop, nfeloop
                ztabl(ji, ij, jk, jl, jf) = ptab(ji, jj, jk)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DO jr = 1, nsndto
        IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
          CALL mppsend(5, znorthloc, itaille, nfipproc(isendto(jr), jpnj), ml_req_nf(jr))
        END IF
      END DO
      DO jr = 1, nsndto
        iproc = nfipproc(isendto(jr), jpnj)
        IF (iproc /= - 1) THEN
          iilb = nimppt(iproc + 1)
          ilci = nlcit(iproc + 1)
          ildi = nldit(iproc + 1)
          ilei = nleit(iproc + 1)
          IF (iilb == 1) ildi = 1
          IF (iilb + ilci - 1 == jpiglo) ilei = ilci
          iilb = nfiimpp(isendto(jr), jpnj) - nfiimpp(isendto(1), jpnj)
        END IF
        IF (iproc /= narea - 1 .AND. iproc /= - 1) THEN
          CALL mpprecv(5, zfoldwk, itaille, iproc)
          DO jf = 1, ipf
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                !$ACC LOOP INDEPENDENT COLLAPSE(2)
                DO jj = 1, ipj
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = zfoldwk(ji, jj, jk, jl, jf)
                  END DO
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END DO
        ELSE IF (iproc == narea - 1) THEN
          DO jf = 1, ipf
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                !$ACC LOOP INDEPENDENT COLLAPSE(2)
                DO jj = 1, ipj
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = ptab(ji, nlcj - ipj + jj, jk)
                  END DO
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END DO
        END IF
      END DO
      IF (l_isend) THEN
        DO jr = 1, nsndto
          IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
            CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
          END IF
        END DO
      END IF
      DO jf = 1, ipf
        CALL lbc_nfd_nogather(ztabl(:, :, :, :, jf), ztabr(:, :, :, :, jf), cd_nat, psgn)
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(ji, jj, jk) = ztabl(ji, ij, jk, jl, jf)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DEALLOCATE(zfoldwk)
      DEALLOCATE(ztabl, ztabr)
    ELSE
      ALLOCATE(ztab(jpiglo, 4, ipk, ipl, ipf))
      ALLOCATE(znorthgloio(jpimax, 4, ipk, ipl, ipf, jpni))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztab(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      CALL MPI_ALLGATHER(znorthloc, itaille, MPI_DOUBLE_PRECISION, znorthgloio, itaille, MPI_DOUBLE_PRECISION, ncomm_north, ierr)
      DO jr = 1, ndim_rank_north
        iproc = nrank_north(jr) + 1
        iilb = nimppt(iproc)
        ilci = nlcit(iproc)
        ildi = nldit(iproc)
        ilei = nleit(iproc)
        IF (iilb == 1) ildi = 1
        IF (iilb + ilci - 1 == jpiglo) ilei = ilci
        DO jf = 1, ipf
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              !$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO jj = 1, ipj
                DO ji = ildi, ilei
                  ztab(ji + iilb - 1, jj, jk, jl, jf) = znorthgloio(ji, jj, jk, jl, jf, jr)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      DO jf = 1, ipf
        CALL lbc_nfd(ztab(:, :, :, :, jf), cd_nat, psgn)
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(ji, jj, jk) = ztab(ji + nimpp - 1, ij, jk, jl, jf)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DEALLOCATE(ztab)
      DEALLOCATE(znorthgloio)
    END IF
    DEALLOCATE(znorthloc)
  END SUBROUTINE mpp_nfd_3d
  SUBROUTINE mpp_nfd_3d_ptr(ptab, cd_nat, psgn, kfld)
    TYPE(PTR_3D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf, jr
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ierr, itaille, ilci, ildi, ilei, iilb
    INTEGER :: ij, iproc
    INTEGER, DIMENSION(jpmaxngh) :: ml_req_nf
    INTEGER :: ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: ztab, ztabl, ztabr
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: znorthloc, zfoldwk
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: znorthgloio
    ipk = SIZE(ptab(1) % pt3d, 3)
    ipl = 1
    ipf = kfld
    ipj = 4
    ALLOCATE(znorthloc(jpimax, 4, ipk, ipl, ipf))
    !$ACC KERNELS
    znorthloc(:, :, :, :, :) = 0._wp
    !$ACC END KERNELS
    DO jf = 1, ipf
      DO jl = 1, ipl
        DO jk = 1, ipk
          DO jj = nlcj - ipj + 1, nlcj
            ij = jj - nlcj + ipj
            znorthloc(1 : jpi, ij, jk, jl, jf) = ptab(jf) % pt3d(1 : jpi, jj, jk)
          END DO
        END DO
      END DO
    END DO
    itaille = jpimax * ipj * ipk * ipl * ipf
    IF (l_north_nogather) THEN
      ALLOCATE(zfoldwk(jpimax, 4, ipk, ipl, ipf))
      ALLOCATE(ztabl(jpimax, 4, ipk, ipl, ipf), ztabr(jpimax * jpmaxngh, 4, ipk, ipl, ipf))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) THEN
        ztabr(:, :, :, :, :) = 0._wp
        ztabl(:, :, :, :, :) = 0._wp
      END IF
      !$ACC END KERNELS
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = nfsloop, nfeloop
                ztabl(ji, ij, jk, jl, jf) = ptab(jf) % pt3d(ji, jj, jk)
              END DO
            END DO
          END DO
        END DO
      END DO
      DO jr = 1, nsndto
        IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
          CALL mppsend(5, znorthloc, itaille, nfipproc(isendto(jr), jpnj), ml_req_nf(jr))
        END IF
      END DO
      DO jr = 1, nsndto
        iproc = nfipproc(isendto(jr), jpnj)
        IF (iproc /= - 1) THEN
          iilb = nimppt(iproc + 1)
          ilci = nlcit(iproc + 1)
          ildi = nldit(iproc + 1)
          ilei = nleit(iproc + 1)
          IF (iilb == 1) ildi = 1
          IF (iilb + ilci - 1 == jpiglo) ilei = ilci
          iilb = nfiimpp(isendto(jr), jpnj) - nfiimpp(isendto(1), jpnj)
        END IF
        IF (iproc /= narea - 1 .AND. iproc /= - 1) THEN
          CALL mpprecv(5, zfoldwk, itaille, iproc)
          DO jf = 1, ipf
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                !$ACC LOOP INDEPENDENT COLLAPSE(2)
                DO jj = 1, ipj
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = zfoldwk(ji, jj, jk, jl, jf)
                  END DO
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END DO
        ELSE IF (iproc == narea - 1) THEN
          DO jf = 1, ipf
            DO jl = 1, ipl
              DO jk = 1, ipk
                DO jj = 1, ipj
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = ptab(jf) % pt3d(ji, nlcj - ipj + jj, jk)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
      END DO
      IF (l_isend) THEN
        DO jr = 1, nsndto
          IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
            CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
          END IF
        END DO
      END IF
      DO jf = 1, ipf
        CALL lbc_nfd_nogather(ztabl(:, :, :, :, jf), ztabr(:, :, :, :, jf), cd_nat(jf), psgn(jf))
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(jf) % pt3d(ji, jj, jk) = ztabl(ji, ij, jk, jl, jf)
              END DO
            END DO
          END DO
        END DO
      END DO
      DEALLOCATE(zfoldwk)
      DEALLOCATE(ztabl, ztabr)
    ELSE
      ALLOCATE(ztab(jpiglo, 4, ipk, ipl, ipf))
      ALLOCATE(znorthgloio(jpimax, 4, ipk, ipl, ipf, jpni))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztab(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      CALL MPI_ALLGATHER(znorthloc, itaille, MPI_DOUBLE_PRECISION, znorthgloio, itaille, MPI_DOUBLE_PRECISION, ncomm_north, ierr)
      DO jr = 1, ndim_rank_north
        iproc = nrank_north(jr) + 1
        iilb = nimppt(iproc)
        ilci = nlcit(iproc)
        ildi = nldit(iproc)
        ilei = nleit(iproc)
        IF (iilb == 1) ildi = 1
        IF (iilb + ilci - 1 == jpiglo) ilei = ilci
        DO jf = 1, ipf
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              !$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO jj = 1, ipj
                DO ji = ildi, ilei
                  ztab(ji + iilb - 1, jj, jk, jl, jf) = znorthgloio(ji, jj, jk, jl, jf, jr)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      DO jf = 1, ipf
        CALL lbc_nfd(ztab(:, :, :, :, jf), cd_nat(jf), psgn(jf))
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(jf) % pt3d(ji, jj, jk) = ztab(ji + nimpp - 1, ij, jk, jl, jf)
              END DO
            END DO
          END DO
        END DO
      END DO
      DEALLOCATE(ztab)
      DEALLOCATE(znorthgloio)
    END IF
    DEALLOCATE(znorthloc)
  END SUBROUTINE mpp_nfd_3d_ptr
  SUBROUTINE mpp_nfd_4d(ptab, cd_nat, psgn, kfld)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf, jr
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ierr, itaille, ilci, ildi, ilei, iilb
    INTEGER :: ij, iproc
    INTEGER, DIMENSION(jpmaxngh) :: ml_req_nf
    INTEGER :: ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: ztab, ztabl, ztabr
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: znorthloc, zfoldwk
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: znorthgloio
    ipk = SIZE(ptab, 3)
    ipl = SIZE(ptab, 4)
    ipf = 1
    ipj = 4
    ALLOCATE(znorthloc(jpimax, 4, ipk, ipl, ipf))
    !$ACC KERNELS
    znorthloc(:, :, :, :, :) = 0._wp
    !$ACC END KERNELS
    DO jf = 1, ipf
      DO jl = 1, ipl
        !$ACC KERNELS
        DO jk = 1, ipk
          DO jj = nlcj - ipj + 1, nlcj
            ij = jj - nlcj + ipj
            znorthloc(1 : jpi, ij, jk, jl, jf) = ptab(1 : jpi, jj, jk, jl)
          END DO
        END DO
        !$ACC END KERNELS
      END DO
    END DO
    itaille = jpimax * ipj * ipk * ipl * ipf
    IF (l_north_nogather) THEN
      ALLOCATE(zfoldwk(jpimax, 4, ipk, ipl, ipf))
      ALLOCATE(ztabl(jpimax, 4, ipk, ipl, ipf), ztabr(jpimax * jpmaxngh, 4, ipk, ipl, ipf))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) THEN
        ztabr(:, :, :, :, :) = 0._wp
        ztabl(:, :, :, :, :) = 0._wp
      END IF
      !$ACC END KERNELS
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = nfsloop, nfeloop
                ztabl(ji, ij, jk, jl, jf) = ptab(ji, jj, jk, jl)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DO jr = 1, nsndto
        IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
          CALL mppsend(5, znorthloc, itaille, nfipproc(isendto(jr), jpnj), ml_req_nf(jr))
        END IF
      END DO
      DO jr = 1, nsndto
        iproc = nfipproc(isendto(jr), jpnj)
        IF (iproc /= - 1) THEN
          iilb = nimppt(iproc + 1)
          ilci = nlcit(iproc + 1)
          ildi = nldit(iproc + 1)
          ilei = nleit(iproc + 1)
          IF (iilb == 1) ildi = 1
          IF (iilb + ilci - 1 == jpiglo) ilei = ilci
          iilb = nfiimpp(isendto(jr), jpnj) - nfiimpp(isendto(1), jpnj)
        END IF
        IF (iproc /= narea - 1 .AND. iproc /= - 1) THEN
          CALL mpprecv(5, zfoldwk, itaille, iproc)
          DO jf = 1, ipf
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                !$ACC LOOP INDEPENDENT COLLAPSE(2)
                DO jj = 1, ipj
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = zfoldwk(ji, jj, jk, jl, jf)
                  END DO
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END DO
        ELSE IF (iproc == narea - 1) THEN
          DO jf = 1, ipf
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                !$ACC LOOP INDEPENDENT COLLAPSE(2)
                DO jj = 1, ipj
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = ptab(ji, nlcj - ipj + jj, jk, jl)
                  END DO
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END DO
        END IF
      END DO
      IF (l_isend) THEN
        DO jr = 1, nsndto
          IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
            CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
          END IF
        END DO
      END IF
      DO jf = 1, ipf
        CALL lbc_nfd_nogather(ztabl(:, :, :, :, jf), ztabr(:, :, :, :, jf), cd_nat, psgn)
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(ji, jj, jk, jl) = ztabl(ji, ij, jk, jl, jf)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DEALLOCATE(zfoldwk)
      DEALLOCATE(ztabl, ztabr)
    ELSE
      ALLOCATE(ztab(jpiglo, 4, ipk, ipl, ipf))
      ALLOCATE(znorthgloio(jpimax, 4, ipk, ipl, ipf, jpni))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztab(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      CALL MPI_ALLGATHER(znorthloc, itaille, MPI_DOUBLE_PRECISION, znorthgloio, itaille, MPI_DOUBLE_PRECISION, ncomm_north, ierr)
      DO jr = 1, ndim_rank_north
        iproc = nrank_north(jr) + 1
        iilb = nimppt(iproc)
        ilci = nlcit(iproc)
        ildi = nldit(iproc)
        ilei = nleit(iproc)
        IF (iilb == 1) ildi = 1
        IF (iilb + ilci - 1 == jpiglo) ilei = ilci
        DO jf = 1, ipf
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              !$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO jj = 1, ipj
                DO ji = ildi, ilei
                  ztab(ji + iilb - 1, jj, jk, jl, jf) = znorthgloio(ji, jj, jk, jl, jf, jr)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      DO jf = 1, ipf
        CALL lbc_nfd(ztab(:, :, :, :, jf), cd_nat, psgn)
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(ji, jj, jk, jl) = ztab(ji + nimpp - 1, ij, jk, jl, jf)
              END DO
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END DO
      DEALLOCATE(ztab)
      DEALLOCATE(znorthgloio)
    END IF
    DEALLOCATE(znorthloc)
  END SUBROUTINE mpp_nfd_4d
  SUBROUTINE mpp_nfd_4d_ptr(ptab, cd_nat, psgn, kfld)
    TYPE(PTR_4D), INTENT(INOUT) :: ptab(:)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat(:)
    REAL(KIND = wp), INTENT(IN) :: psgn(:)
    INTEGER, OPTIONAL, INTENT(IN) :: kfld
    INTEGER :: ji, jj, jk, jl, jh, jf, jr
    INTEGER :: ipi, ipj, ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ierr, itaille, ilci, ildi, ilei, iilb
    INTEGER :: ij, iproc
    INTEGER, DIMENSION(jpmaxngh) :: ml_req_nf
    INTEGER :: ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: ztab, ztabl, ztabr
    REAL(KIND = wp), DIMENSION(:, :, :, :, :), ALLOCATABLE :: znorthloc, zfoldwk
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: znorthgloio
    ipk = SIZE(ptab(1) % pt4d, 3)
    ipl = SIZE(ptab(1) % pt4d, 4)
    ipf = kfld
    ipj = 4
    ALLOCATE(znorthloc(jpimax, 4, ipk, ipl, ipf))
    !$ACC KERNELS
    znorthloc(:, :, :, :, :) = 0._wp
    !$ACC END KERNELS
    DO jf = 1, ipf
      DO jl = 1, ipl
        DO jk = 1, ipk
          DO jj = nlcj - ipj + 1, nlcj
            ij = jj - nlcj + ipj
            znorthloc(1 : jpi, ij, jk, jl, jf) = ptab(jf) % pt4d(1 : jpi, jj, jk, jl)
          END DO
        END DO
      END DO
    END DO
    itaille = jpimax * ipj * ipk * ipl * ipf
    IF (l_north_nogather) THEN
      ALLOCATE(zfoldwk(jpimax, 4, ipk, ipl, ipf))
      ALLOCATE(ztabl(jpimax, 4, ipk, ipl, ipf), ztabr(jpimax * jpmaxngh, 4, ipk, ipl, ipf))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) THEN
        ztabr(:, :, :, :, :) = 0._wp
        ztabl(:, :, :, :, :) = 0._wp
      END IF
      !$ACC END KERNELS
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = nfsloop, nfeloop
                ztabl(ji, ij, jk, jl, jf) = ptab(jf) % pt4d(ji, jj, jk, jl)
              END DO
            END DO
          END DO
        END DO
      END DO
      DO jr = 1, nsndto
        IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
          CALL mppsend(5, znorthloc, itaille, nfipproc(isendto(jr), jpnj), ml_req_nf(jr))
        END IF
      END DO
      DO jr = 1, nsndto
        iproc = nfipproc(isendto(jr), jpnj)
        IF (iproc /= - 1) THEN
          iilb = nimppt(iproc + 1)
          ilci = nlcit(iproc + 1)
          ildi = nldit(iproc + 1)
          ilei = nleit(iproc + 1)
          IF (iilb == 1) ildi = 1
          IF (iilb + ilci - 1 == jpiglo) ilei = ilci
          iilb = nfiimpp(isendto(jr), jpnj) - nfiimpp(isendto(1), jpnj)
        END IF
        IF (iproc /= narea - 1 .AND. iproc /= - 1) THEN
          CALL mpprecv(5, zfoldwk, itaille, iproc)
          DO jf = 1, ipf
            DO jl = 1, ipl
              !$ACC KERNELS
              DO jk = 1, ipk
                !$ACC LOOP INDEPENDENT COLLAPSE(2)
                DO jj = 1, ipj
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = zfoldwk(ji, jj, jk, jl, jf)
                  END DO
                END DO
              END DO
              !$ACC END KERNELS
            END DO
          END DO
        ELSE IF (iproc == narea - 1) THEN
          DO jf = 1, ipf
            DO jl = 1, ipl
              DO jk = 1, ipk
                DO jj = 1, ipj
                  DO ji = ildi, ilei
                    ztabr(iilb + ji, jj, jk, jl, jf) = ptab(jf) % pt4d(ji, nlcj - ipj + jj, jk, jl)
                  END DO
                END DO
              END DO
            END DO
          END DO
        END IF
      END DO
      IF (l_isend) THEN
        DO jr = 1, nsndto
          IF (nfipproc(isendto(jr), jpnj) /= narea - 1 .AND. nfipproc(isendto(jr), jpnj) /= - 1) THEN
            CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
          END IF
        END DO
      END IF
      DO jf = 1, ipf
        CALL lbc_nfd_nogather(ztabl(:, :, :, :, jf), ztabr(:, :, :, :, jf), cd_nat(jf), psgn(jf))
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(jf) % pt4d(ji, jj, jk, jl) = ztabl(ji, ij, jk, jl, jf)
              END DO
            END DO
          END DO
        END DO
      END DO
      DEALLOCATE(zfoldwk)
      DEALLOCATE(ztabl, ztabr)
    ELSE
      ALLOCATE(ztab(jpiglo, 4, ipk, ipl, ipf))
      ALLOCATE(znorthgloio(jpimax, 4, ipk, ipl, ipf, jpni))
      !$ACC KERNELS
      IF (jpni * jpnj /= jpnij) ztab(:, :, :, :, :) = 0._wp
      !$ACC END KERNELS
      CALL MPI_ALLGATHER(znorthloc, itaille, MPI_DOUBLE_PRECISION, znorthgloio, itaille, MPI_DOUBLE_PRECISION, ncomm_north, ierr)
      DO jr = 1, ndim_rank_north
        iproc = nrank_north(jr) + 1
        iilb = nimppt(iproc)
        ilci = nlcit(iproc)
        ildi = nldit(iproc)
        ilei = nleit(iproc)
        IF (iilb == 1) ildi = 1
        IF (iilb + ilci - 1 == jpiglo) ilei = ilci
        DO jf = 1, ipf
          DO jl = 1, ipl
            !$ACC KERNELS
            DO jk = 1, ipk
              !$ACC LOOP INDEPENDENT COLLAPSE(2)
              DO jj = 1, ipj
                DO ji = ildi, ilei
                  ztab(ji + iilb - 1, jj, jk, jl, jf) = znorthgloio(ji, jj, jk, jl, jf, jr)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          END DO
        END DO
      END DO
      DO jf = 1, ipf
        CALL lbc_nfd(ztab(:, :, :, :, jf), cd_nat(jf), psgn(jf))
      END DO
      DO jf = 1, ipf
        DO jl = 1, ipl
          DO jk = 1, ipk
            DO jj = nlcj - ipj + 1, nlcj
              ij = jj - nlcj + ipj
              DO ji = 1, nlci
                ptab(jf) % pt4d(ji, jj, jk, jl) = ztab(ji + nimpp - 1, ij, jk, jl, jf)
              END DO
            END DO
          END DO
        END DO
      END DO
      DEALLOCATE(ztab)
      DEALLOCATE(znorthgloio)
    END IF
    DEALLOCATE(znorthloc)
  END SUBROUTINE mpp_nfd_4d_ptr
  SUBROUTINE mpp_lnk_bdy_2d(ptab, cd_nat, psgn, kb_bdy)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: kb_bdy
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = 1
    ipl = 1
    ipf = 1
    ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2), zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), &
&zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    zland = 0._wp
    DO jf = 1, ipf
      IF (nbondi == 2) THEN
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :) = ptab(jpim1, :)
          ptab(jpi, :) = ptab(2, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :) = zland
          ptab(nlci - nn_hls + 1 : jpi, :) = zland
          !$ACC END KERNELS
        END IF
      ELSE IF (nbondi == - 1) THEN
        !$ACC KERNELS
        IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :) = zland
        !$ACC END KERNELS
      ELSE IF (nbondi == 1) THEN
        ptab(nlci - nn_hls + 1 : jpi, :) = zland
      END IF
      IF (nbondj == 2) THEN
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1) = ptab(:, jpjm1)
          ptab(:, jpj) = ptab(:, 2)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls) = zland
          ptab(:, nlcj - nn_hls + 1 : jpj) = zland
          !$ACC END KERNELS
        END IF
      ELSE IF (nbondj == - 1) THEN
        !$ACC KERNELS
        IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls) = zland
        !$ACC END KERNELS
      ELSE IF (nbondj == 1) THEN
        ptab(:, nlcj - nn_hls + 1 : jpj) = zland
      END IF
    END DO
    DO jf = 1, ipf
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1, 0, 1)
        iihom = nlci - nreci
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
      imigr = nn_hls * jpj * ipk * ipl
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1)
        CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
        CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE (0)
        CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
        CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
        CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
        CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
        IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE (1)
        CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
        CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      iihom = nlci - nn_hls
      SELECT CASE (nbondi_bdy_b(kb_bdy))
      CASE (- 1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(iihom + jh, :) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (0)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(iihom + jh, :) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :) = zt3we(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
    END DO
    DO jf = 1, ipf
      IF (nbondj_bdy(kb_bdy) /= 2) THEN
        ijhom = nlcj - nrecj
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END IF
      imigr = nn_hls * jpi * ipk * ipl
      SELECT CASE (nbondj_bdy(kb_bdy))
      CASE (- 1)
        CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
        CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE (0)
        CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
        CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
        CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
        CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
        IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE (1)
        CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
        CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      ijhom = nlcj - nn_hls
      SELECT CASE (nbondj_bdy_b(kb_bdy))
      CASE (- 1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, ijhom + jh) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (0)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(:, ijhom + jh) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh) = zt3sn(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
    END DO
    IF (npolj /= 0) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat, psgn)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat, psgn)
      END SELECT
    END IF
    DEALLOCATE(zt3ns, zt3sn, zt3ew, zt3we)
  END SUBROUTINE mpp_lnk_bdy_2d
  SUBROUTINE mpp_lnk_bdy_3d(ptab, cd_nat, psgn, kb_bdy)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: kb_bdy
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = SIZE(ptab, 3)
    ipl = 1
    ipf = 1
    ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2), zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), &
&zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    zland = 0._wp
    DO jf = 1, ipf
      IF (nbondi == 2) THEN
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :, :) = ptab(jpim1, :, :)
          ptab(jpi, :, :) = ptab(2, :, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :, :) = zland
          ptab(nlci - nn_hls + 1 : jpi, :, :) = zland
          !$ACC END KERNELS
        END IF
      ELSE IF (nbondi == - 1) THEN
        !$ACC KERNELS
        IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :, :) = zland
        !$ACC END KERNELS
      ELSE IF (nbondi == 1) THEN
        ptab(nlci - nn_hls + 1 : jpi, :, :) = zland
      END IF
      IF (nbondj == 2) THEN
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1, :) = ptab(:, jpjm1, :)
          ptab(:, jpj, :) = ptab(:, 2, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls, :) = zland
          ptab(:, nlcj - nn_hls + 1 : jpj, :) = zland
          !$ACC END KERNELS
        END IF
      ELSE IF (nbondj == - 1) THEN
        !$ACC KERNELS
        IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls, :) = zland
        !$ACC END KERNELS
      ELSE IF (nbondj == 1) THEN
        ptab(:, nlcj - nn_hls + 1 : jpj, :) = zland
      END IF
    END DO
    DO jf = 1, ipf
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1, 0, 1)
        iihom = nlci - nreci
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :, jk)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
      imigr = nn_hls * jpj * ipk * ipl
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1)
        CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
        CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE (0)
        CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
        CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
        CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
        CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
        IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE (1)
        CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
        CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      iihom = nlci - nn_hls
      SELECT CASE (nbondi_bdy_b(kb_bdy))
      CASE (- 1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(iihom + jh, :, jk) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (0)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(iihom + jh, :, jk) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk) = zt3we(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
    END DO
    DO jf = 1, ipf
      IF (nbondj_bdy(kb_bdy) /= 2) THEN
        ijhom = nlcj - nrecj
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh, jk)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh, jk)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END IF
      imigr = nn_hls * jpi * ipk * ipl
      SELECT CASE (nbondj_bdy(kb_bdy))
      CASE (- 1)
        CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
        CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE (0)
        CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
        CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
        CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
        CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
        IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE (1)
        CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
        CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      ijhom = nlcj - nn_hls
      SELECT CASE (nbondj_bdy_b(kb_bdy))
      CASE (- 1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, ijhom + jh, jk) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (0)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(:, ijhom + jh, jk) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk) = zt3sn(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
    END DO
    IF (npolj /= 0) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat, psgn)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat, psgn)
      END SELECT
    END IF
    DEALLOCATE(zt3ns, zt3sn, zt3ew, zt3we)
  END SUBROUTINE mpp_lnk_bdy_3d
  SUBROUTINE mpp_lnk_bdy_4d(ptab, cd_nat, psgn, kb_bdy)
    REAL(KIND = wp), INTENT(INOUT) :: ptab(:, :, :, :)
    CHARACTER(LEN = 1), INTENT(IN) :: cd_nat
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: kb_bdy
    INTEGER :: ji, jj, jk, jl, jh, jf
    INTEGER :: ipk, ipl, ipf
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ml_req1, ml_req2, ml_err
    REAL(KIND = wp) :: zland
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ns, zt3sn
    REAL(KIND = wp), DIMENSION(:, :, :, :, :, :), ALLOCATABLE :: zt3ew, zt3we
    ipk = SIZE(ptab, 3)
    ipl = SIZE(ptab, 4)
    ipf = 1
    ALLOCATE(zt3ns(jpi, nn_hls, ipk, ipl, ipf, 2), zt3sn(jpi, nn_hls, ipk, ipl, ipf, 2), zt3ew(jpj, nn_hls, ipk, ipl, ipf, 2), &
&zt3we(jpj, nn_hls, ipk, ipl, ipf, 2))
    zland = 0._wp
    DO jf = 1, ipf
      IF (nbondi == 2) THEN
        IF (l_Iperio) THEN
          !$ACC KERNELS
          ptab(1, :, :, :) = ptab(jpim1, :, :, :)
          ptab(jpi, :, :, :) = ptab(2, :, :, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :, :, :) = zland
          ptab(nlci - nn_hls + 1 : jpi, :, :, :) = zland
          !$ACC END KERNELS
        END IF
      ELSE IF (nbondi == - 1) THEN
        !$ACC KERNELS
        IF (.NOT. cd_nat == 'F') ptab(1 : nn_hls, :, :, :) = zland
        !$ACC END KERNELS
      ELSE IF (nbondi == 1) THEN
        ptab(nlci - nn_hls + 1 : jpi, :, :, :) = zland
      END IF
      IF (nbondj == 2) THEN
        IF (l_Jperio) THEN
          !$ACC KERNELS
          ptab(:, 1, :, :) = ptab(:, jpjm1, :, :)
          ptab(:, jpj, :, :) = ptab(:, 2, :, :)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls, :, :) = zland
          ptab(:, nlcj - nn_hls + 1 : jpj, :, :) = zland
          !$ACC END KERNELS
        END IF
      ELSE IF (nbondj == - 1) THEN
        !$ACC KERNELS
        IF (.NOT. cd_nat == 'F') ptab(:, 1 : nn_hls, :, :) = zland
        !$ACC END KERNELS
      ELSE IF (nbondj == 1) THEN
        ptab(:, nlcj - nn_hls + 1 : jpj, :, :) = zland
      END IF
    END DO
    DO jf = 1, ipf
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1, 0, 1)
        iihom = nlci - nreci
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3ew(:, jh, jk, jl, jf, 1) = ptab(nn_hls + jh, :, jk, jl)
              zt3we(:, jh, jk, jl, jf, 1) = ptab(iihom + jh, :, jk, jl)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
      imigr = nn_hls * jpj * ipk * ipl
      SELECT CASE (nbondi_bdy(kb_bdy))
      CASE (- 1)
        CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req1)
        CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE (0)
        CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
        CALL mppsend(2, zt3we(1, 1, 1, 1, 1, 1), imigr, noea, ml_req2)
        CALL mpprecv(1, zt3ew(1, 1, 1, 1, 1, 2), imigr, noea)
        CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
        IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE (1)
        CALL mppsend(1, zt3ew(1, 1, 1, 1, 1, 1), imigr, nowe, ml_req1)
        CALL mpprecv(2, zt3we(1, 1, 1, 1, 1, 2), imigr, nowe)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      iihom = nlci - nn_hls
      SELECT CASE (nbondi_bdy_b(kb_bdy))
      CASE (- 1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(iihom + jh, :, jk, jl) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (0)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk, jl) = zt3we(:, jh, jk, jl, jf, 2)
              ptab(iihom + jh, :, jk, jl) = zt3ew(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(jh, :, jk, jl) = zt3we(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
    END DO
    DO jf = 1, ipf
      IF (nbondj_bdy(kb_bdy) /= 2) THEN
        ijhom = nlcj - nrecj
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              zt3sn(:, jh, jk, jl, jf, 1) = ptab(:, ijhom + jh, jk, jl)
              zt3ns(:, jh, jk, jl, jf, 1) = ptab(:, nn_hls + jh, jk, jl)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END IF
      imigr = nn_hls * jpi * ipk * ipl
      SELECT CASE (nbondj_bdy(kb_bdy))
      CASE (- 1)
        CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req1)
        CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE (0)
        CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
        CALL mppsend(4, zt3sn(1, 1, 1, 1, 1, 1), imigr, nono, ml_req2)
        CALL mpprecv(3, zt3ns(1, 1, 1, 1, 1, 2), imigr, nono)
        CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
        IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE (1)
        CALL mppsend(3, zt3ns(1, 1, 1, 1, 1, 1), imigr, noso, ml_req1)
        CALL mpprecv(4, zt3sn(1, 1, 1, 1, 1, 2), imigr, noso)
        IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      ijhom = nlcj - nn_hls
      SELECT CASE (nbondj_bdy_b(kb_bdy))
      CASE (- 1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, ijhom + jh, jk, jl) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (0)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk, jl) = zt3sn(:, jh, jk, jl, jf, 2)
              ptab(:, ijhom + jh, jk, jl) = zt3ns(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      CASE (1)
        DO jl = 1, ipl
          !$ACC KERNELS
          DO jk = 1, ipk
            DO jh = 1, nn_hls
              ptab(:, jh, jk, jl) = zt3sn(:, jh, jk, jl, jf, 2)
            END DO
          END DO
          !$ACC END KERNELS
        END DO
      END SELECT
    END DO
    IF (npolj /= 0) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(ptab, cd_nat, psgn)
      CASE DEFAULT
        CALL mpp_nfd(ptab, cd_nat, psgn)
      END SELECT
    END IF
    DEALLOCATE(zt3ns, zt3sn, zt3ew, zt3we)
  END SUBROUTINE mpp_lnk_bdy_4d
  SUBROUTINE mppsend(ktyp, pmess, kbytes, kdest, md_req)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(INOUT) :: pmess(*)
    INTEGER, INTENT(IN) :: kbytes
    INTEGER, INTENT(IN) :: kdest
    INTEGER, INTENT(IN) :: ktyp
    INTEGER, INTENT(IN) :: md_req
    INTEGER :: iflag
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppsend', 'r0', 0, 0)
    SELECT CASE (cn_mpi_send)
    CASE ('S')
      CALL mpi_send(pmess, kbytes, mpi_double_precision, kdest, ktyp, mpi_comm_oce, iflag)
    CASE ('B')
      CALL mpi_bsend(pmess, kbytes, mpi_double_precision, kdest, ktyp, mpi_comm_oce, iflag)
    CASE ('I')
      CALL mpi_isend(pmess, kbytes, mpi_double_precision, kdest, ktyp, mpi_comm_oce, md_req, iflag)
    END SELECT
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mppsend
  SUBROUTINE mpprecv(ktyp, pmess, kbytes, ksource)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(INOUT) :: pmess(*)
    INTEGER, INTENT(IN) :: kbytes
    INTEGER, INTENT(IN) :: ktyp
    INTEGER, OPTIONAL, INTENT(IN) :: ksource
    INTEGER :: istatus(mpi_status_size)
    INTEGER :: iflag
    INTEGER :: use_source
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mpprecv', 'r0', 0, 0)
    use_source = mpi_any_source
    IF (PRESENT(ksource)) use_source = ksource
    CALL mpi_recv(pmess, kbytes, mpi_double_precision, use_source, ktyp, mpi_comm_oce, istatus, iflag)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mpprecv
  SUBROUTINE mppgather(ptab, kp, pio)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: ptab
    INTEGER, INTENT(IN) :: kp
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpnij), INTENT(OUT) :: pio
    INTEGER :: itaille, ierror
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppgather', 'r0', 0, 0)
    itaille = jpi * jpj
    CALL mpi_gather(ptab, itaille, mpi_double_precision, pio, itaille, mpi_double_precision, kp, mpi_comm_oce, ierror)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mppgather
  SUBROUTINE mppscatter(pio, kp, ptab)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpnij) :: pio
    INTEGER :: kp
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ptab
    INTEGER :: itaille, ierror
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppscatter', 'r0', 0, 0)
    itaille = jpi * jpj
    CALL mpi_scatter(pio, itaille, mpi_double_precision, ptab, itaille, mpi_double_precision, kp, mpi_comm_oce, ierror)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mppscatter
  SUBROUTINE mppmax_a_int(ktab, kdim, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kdim
    INTEGER, INTENT(INOUT), DIMENSION(kdim) :: ktab
    INTEGER, INTENT(IN), OPTIONAL :: kcom
    INTEGER :: ierror, ilocalcomm
    INTEGER, DIMENSION(kdim) :: iwork
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmax_a_int', 'r0', 0, 0)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    CALL mpi_allreduce(ktab, iwork, kdim, mpi_integer, mpi_max, ilocalcomm, ierror)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    ktab(:) = iwork(:)
    !$ACC END KERNELS
  END SUBROUTINE mppmax_a_int
  SUBROUTINE mppmax_int(ktab, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(INOUT) :: ktab
    INTEGER, INTENT(IN), OPTIONAL :: kcom
    INTEGER :: ierror, iwork, ilocalcomm
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmax_int', 'r0', 0, 0)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    CALL mpi_allreduce(ktab, iwork, 1, mpi_integer, mpi_max, ilocalcomm, ierror)
    ktab = iwork
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mppmax_int
  SUBROUTINE mppmax_a_real(ptab, kdim, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(kdim), INTENT(INOUT) :: ptab
    INTEGER, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ierror, ilocalcomm
    REAL(KIND = wp), DIMENSION(kdim) :: zwork
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmax_a_real', 'r0', 0, 0)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    CALL mpi_allreduce(ptab, zwork, kdim, mpi_double_precision, mpi_max, ilocalcomm, ierror)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    ptab(:) = zwork(:)
    !$ACC END KERNELS
  END SUBROUTINE mppmax_a_real
  SUBROUTINE mppmax_real(ptab, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(INOUT) :: ptab
    INTEGER, INTENT(IN), OPTIONAL :: kcom
    INTEGER :: ierror, ilocalcomm
    REAL(KIND = wp) :: zwork
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmax_real', 'r0', 0, 0)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    CALL mpi_allreduce(ptab, zwork, 1, mpi_double_precision, mpi_max, ilocalcomm, ierror)
    ptab = zwork
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mppmax_real
  SUBROUTINE mppmin_a_int(ktab, kdim, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kdim
    INTEGER, INTENT(INOUT), DIMENSION(kdim) :: ktab
    INTEGER, INTENT(IN), OPTIONAL :: kcom
    INTEGER :: ierror, ilocalcomm
    INTEGER, DIMENSION(kdim) :: iwork
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmin_a_int', 'r0', 0, 0)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    CALL mpi_allreduce(ktab, iwork, kdim, mpi_integer, mpi_min, ilocalcomm, ierror)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    ktab(:) = iwork(:)
    !$ACC END KERNELS
  END SUBROUTINE mppmin_a_int
  SUBROUTINE mppmin_int(ktab, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(INOUT) :: ktab
    INTEGER, INTENT(IN), OPTIONAL :: kcom
    INTEGER :: ierror, iwork, ilocalcomm
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmin_int', 'r0', 0, 0)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    CALL mpi_allreduce(ktab, iwork, 1, mpi_integer, mpi_min, ilocalcomm, ierror)
    ktab = iwork
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mppmin_int
  SUBROUTINE mppmin_a_real(ptab, kdim, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: kdim
    REAL(KIND = wp), INTENT(INOUT), DIMENSION(kdim) :: ptab
    INTEGER, INTENT(IN), OPTIONAL :: kcom
    INTEGER :: ierror, ilocalcomm
    REAL(KIND = wp), DIMENSION(kdim) :: zwork
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmin_a_real', 'r0', 0, 0)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    CALL mpi_allreduce(ptab, zwork, kdim, mpi_double_precision, mpi_min, ilocalcomm, ierror)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    ptab(:) = zwork(:)
    !$ACC END KERNELS
  END SUBROUTINE mppmin_a_real
  SUBROUTINE mppmin_real(ptab, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), INTENT(INOUT) :: ptab
    INTEGER, INTENT(IN), OPTIONAL :: kcom
    INTEGER :: ierror, ilocalcomm
    REAL(KIND = wp) :: zwork
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmin_real', 'r0', 0, 0)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    CALL mpi_allreduce(ptab, zwork, 1, mpi_double_precision, mpi_min, ilocalcomm, ierror)
    ptab = zwork
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mppmin_real
  SUBROUTINE mppsum_a_int(ktab, kdim)
    INTEGER, INTENT(IN) :: kdim
    INTEGER, INTENT(INOUT), DIMENSION(kdim) :: ktab
    INTEGER :: ierror
    INTEGER, DIMENSION(kdim) :: iwork
    CALL mpi_allreduce(ktab, iwork, kdim, mpi_integer, mpi_sum, mpi_comm_oce, ierror)
    !$ACC KERNELS
    ktab(:) = iwork(:)
    !$ACC END KERNELS
  END SUBROUTINE mppsum_a_int
  SUBROUTINE mppsum_int(ktab)
    INTEGER, INTENT(INOUT) :: ktab
    INTEGER :: ierror, iwork
    CALL mpi_allreduce(ktab, iwork, 1, mpi_integer, mpi_sum, mpi_comm_oce, ierror)
    ktab = iwork
  END SUBROUTINE mppsum_int
  SUBROUTINE mppsum_a_real(ptab, kdim, kcom)
    INTEGER, INTENT(IN) :: kdim
    REAL(KIND = wp), DIMENSION(kdim), INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ierror, ilocalcomm
    REAL(KIND = wp) :: zwork(kdim)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    CALL mpi_allreduce(ptab, zwork, kdim, mpi_double_precision, mpi_sum, ilocalcomm, ierror)
    !$ACC KERNELS
    ptab(:) = zwork(:)
    !$ACC END KERNELS
  END SUBROUTINE mppsum_a_real
  SUBROUTINE mppsum_real(ptab, kcom)
    REAL(KIND = wp), INTENT(INOUT) :: ptab
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ierror, ilocalcomm
    REAL(KIND = wp) :: zwork
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    CALL mpi_allreduce(ptab, zwork, 1, mpi_double_precision, mpi_sum, ilocalcomm, ierror)
    ptab = zwork
  END SUBROUTINE mppsum_real
  SUBROUTINE mppsum_realdd(ytab, kcom)
    COMPLEX(KIND = wp), INTENT(INOUT) :: ytab
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ierror, ilocalcomm
    COMPLEX(KIND = wp) :: zwork
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    CALL MPI_ALLREDUCE(ytab, zwork, 1, MPI_DOUBLE_COMPLEX, MPI_SUMDD, ilocalcomm, ierror)
    ytab = zwork
  END SUBROUTINE mppsum_realdd
  SUBROUTINE mppsum_a_realdd(ytab, kdim, kcom)
    INTEGER, INTENT(IN) :: kdim
    COMPLEX(KIND = wp), DIMENSION(kdim), INTENT(INOUT) :: ytab
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ierror, ilocalcomm
    COMPLEX(KIND = wp), DIMENSION(kdim) :: zwork
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    CALL MPI_ALLREDUCE(ytab, zwork, kdim, MPI_DOUBLE_COMPLEX, MPI_SUMDD, ilocalcomm, ierror)
    !$ACC KERNELS
    ytab(:) = zwork(:)
    !$ACC END KERNELS
  END SUBROUTINE mppsum_a_realdd
  SUBROUTINE mppmax_real_multiple(pt1d, kdim, kcom)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    REAL(KIND = wp), DIMENSION(kdim), INTENT(INOUT) :: pt1d
    INTEGER, INTENT(IN) :: kdim
    INTEGER, OPTIONAL, INTENT(IN) :: kcom
    INTEGER :: ierror, ilocalcomm
    REAL(KIND = wp), DIMENSION(kdim) :: zwork
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppmax_real_multiple', 'r0', 0, 0)
    ilocalcomm = mpi_comm_oce
    IF (PRESENT(kcom)) ilocalcomm = kcom
    CALL mpi_allreduce(pt1d, zwork, kdim, mpi_double_precision, mpi_max, ilocalcomm, ierror)
    CALL profile_psy_data0 % PostEnd
    !$ACC KERNELS
    pt1d(:) = zwork(:)
    !$ACC END KERNELS
  END SUBROUTINE mppmax_real_multiple
  SUBROUTINE mpp_minloc2d(ptab, pmask, pmin, ki, kj)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: ptab
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pmask
    REAL(KIND = wp), INTENT(OUT) :: pmin
    INTEGER, INTENT(OUT) :: ki, kj
    INTEGER :: ierror
    INTEGER, DIMENSION(2) :: ilocs
    REAL(KIND = wp) :: zmin
    REAL(KIND = wp), DIMENSION(2, 1) :: zain, zaout
    zmin = MINVAL(ptab(:, :), mask = pmask == 1._wp)
    ilocs = MINLOC(ptab(:, :), mask = pmask == 1._wp)
    !$ACC KERNELS
    ki = ilocs(1) + nimpp - 1
    kj = ilocs(2) + njmpp - 1
    zain(1, :) = zmin
    zain(2, :) = ki + 10000. * kj
    !$ACC END KERNELS
    CALL MPI_ALLREDUCE(zain, zaout, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_OCE, ierror)
    pmin = zaout(1, 1)
    kj = INT(zaout(2, 1) / 10000.)
    ki = INT(zaout(2, 1) - 10000. * kj)
  END SUBROUTINE mpp_minloc2d
  SUBROUTINE mpp_minloc3d(ptab, pmask, pmin, ki, kj, kk)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: ptab
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pmask
    REAL(KIND = wp), INTENT(OUT) :: pmin
    INTEGER, INTENT(OUT) :: ki, kj, kk
    INTEGER :: ierror
    REAL(KIND = wp) :: zmin
    INTEGER, DIMENSION(3) :: ilocs
    REAL(KIND = wp), DIMENSION(2, 1) :: zain, zaout
    zmin = MINVAL(ptab(:, :, :), mask = pmask == 1._wp)
    ilocs = MINLOC(ptab(:, :, :), mask = pmask == 1._wp)
    !$ACC KERNELS
    ki = ilocs(1) + nimpp - 1
    kj = ilocs(2) + njmpp - 1
    kk = ilocs(3)
    zain(1, :) = zmin
    zain(2, :) = ki + 10000. * kj + 100000000. * kk
    !$ACC END KERNELS
    CALL MPI_ALLREDUCE(zain, zaout, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, MPI_COMM_OCE, ierror)
    pmin = zaout(1, 1)
    kk = INT(zaout(2, 1) / 100000000.)
    kj = INT(zaout(2, 1) - kk * 100000000.) / 10000
    ki = INT(zaout(2, 1) - kk * 100000000. - kj * 10000.)
  END SUBROUTINE mpp_minloc3d
  SUBROUTINE mpp_maxloc2d(ptab, pmask, pmax, ki, kj)
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: ptab
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN) :: pmask
    REAL(KIND = wp), INTENT(OUT) :: pmax
    INTEGER, INTENT(OUT) :: ki, kj
    INTEGER :: ierror
    INTEGER, DIMENSION(2) :: ilocs
    REAL(KIND = wp) :: zmax
    REAL(KIND = wp), DIMENSION(2, 1) :: zain, zaout
    zmax = MAXVAL(ptab(:, :), mask = pmask == 1._wp)
    ilocs = MAXLOC(ptab(:, :), mask = pmask == 1._wp)
    !$ACC KERNELS
    ki = ilocs(1) + nimpp - 1
    kj = ilocs(2) + njmpp - 1
    zain(1, :) = zmax
    zain(2, :) = ki + 10000. * kj
    !$ACC END KERNELS
    CALL MPI_ALLREDUCE(zain, zaout, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_OCE, ierror)
    pmax = zaout(1, 1)
    kj = INT(zaout(2, 1) / 10000.)
    ki = INT(zaout(2, 1) - 10000. * kj)
  END SUBROUTINE mpp_maxloc2d
  SUBROUTINE mpp_maxloc3d(ptab, pmask, pmax, ki, kj, kk)
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: ptab
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(IN) :: pmask
    REAL(KIND = wp), INTENT(OUT) :: pmax
    INTEGER, INTENT(OUT) :: ki, kj, kk
    INTEGER :: ierror
    REAL(KIND = wp) :: zmax
    REAL(KIND = wp), DIMENSION(2, 1) :: zain, zaout
    INTEGER, DIMENSION(3) :: ilocs
    zmax = MAXVAL(ptab(:, :, :), mask = pmask == 1._wp)
    ilocs = MAXLOC(ptab(:, :, :), mask = pmask == 1._wp)
    !$ACC KERNELS
    ki = ilocs(1) + nimpp - 1
    kj = ilocs(2) + njmpp - 1
    kk = ilocs(3)
    zain(1, :) = zmax
    zain(2, :) = ki + 10000. * kj + 100000000. * kk
    !$ACC END KERNELS
    CALL MPI_ALLREDUCE(zain, zaout, 1, MPI_2DOUBLE_PRECISION, MPI_MAXLOC, MPI_COMM_OCE, ierror)
    pmax = zaout(1, 1)
    kk = INT(zaout(2, 1) / 100000000.)
    kj = INT(zaout(2, 1) - kk * 100000000.) / 10000
    ki = INT(zaout(2, 1) - kk * 100000000. - kj * 10000.)
  END SUBROUTINE mpp_maxloc3d
  SUBROUTINE mppsync
    INTEGER :: ierror
    CALL mpi_barrier(mpi_comm_oce, ierror)
  END SUBROUTINE mppsync
  SUBROUTINE mppstop
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER :: info
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('mppstop', 'r0', 0, 0)
    CALL mppsync
    CALL mpi_finalize(info)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE mppstop
  SUBROUTINE mpp_comm_free(kcom)
    INTEGER, INTENT(IN) :: kcom
    INTEGER :: ierr
    CALL MPI_COMM_FREE(kcom, ierr)
  END SUBROUTINE mpp_comm_free
  SUBROUTINE mpp_ini_ice(pindic, kumout)
    INTEGER, INTENT(IN) :: pindic
    INTEGER, INTENT(IN) :: kumout
    INTEGER :: jjproc
    INTEGER :: ii, ierr
    INTEGER, ALLOCATABLE, DIMENSION(:) :: kice
    INTEGER, ALLOCATABLE, DIMENSION(:) :: zwork
    ALLOCATE(kice(jpnij), zwork(jpnij), STAT = ierr)
    IF (ierr /= 0) THEN
      WRITE(kumout, cform_err)
      WRITE(kumout, FMT = *) 'mpp_ini_ice : failed to allocate 2, 1D arrays (jpnij in length)'
      CALL mppstop
    END IF
    !$ACC KERNELS
    kice = 0
    DO jjproc = 1, jpnij
      IF (jjproc == narea .AND. pindic .GT. 0) kice(jjproc) = 1
    END DO
    zwork = 0
    !$ACC END KERNELS
    CALL MPI_ALLREDUCE(kice, zwork, jpnij, mpi_integer, mpi_sum, mpi_comm_oce, ierr)
    ndim_rank_ice = SUM(zwork)
    IF (ALLOCATED(nrank_ice)) DEALLOCATE(nrank_ice)
    ALLOCATE(nrank_ice(ndim_rank_ice))
    ii = 0
    nrank_ice = 0
    DO jjproc = 1, jpnij
      IF (zwork(jjproc) == 1) THEN
        ii = ii + 1
        nrank_ice(ii) = jjproc - 1
      END IF
    END DO
    CALL MPI_COMM_GROUP(mpi_comm_oce, ngrp_iworld, ierr)
    CALL MPI_GROUP_INCL(ngrp_iworld, ndim_rank_ice, nrank_ice, ngrp_ice, ierr)
    CALL MPI_COMM_CREATE(mpi_comm_oce, ngrp_ice, ncomm_ice, ierr)
    CALL MPI_GROUP_FREE(ngrp_ice, ierr)
    CALL MPI_GROUP_FREE(ngrp_iworld, ierr)
    DEALLOCATE(kice, zwork)
  END SUBROUTINE mpp_ini_ice
  SUBROUTINE mpp_ini_znl(kumout)
    INTEGER, INTENT(IN) :: kumout
    INTEGER :: jproc
    INTEGER :: ierr, ii
    INTEGER, ALLOCATABLE, DIMENSION(:) :: kwork
    ALLOCATE(kwork(jpnij), STAT = ierr)
    IF (ierr /= 0) THEN
      WRITE(kumout, cform_err)
      WRITE(kumout, FMT = *) 'mpp_ini_znl : failed to allocate 1D array of length jpnij'
      CALL mppstop
    END IF
    IF (jpnj == 1) THEN
      ngrp_znl = ngrp_world
      ncomm_znl = mpi_comm_oce
    ELSE
      CALL MPI_ALLGATHER(njmpp, 1, mpi_integer, kwork, 1, mpi_integer, mpi_comm_oce, ierr)
      ndim_rank_znl = 0
      DO jproc = 1, jpnij
        IF (kwork(jproc) == njmpp) THEN
          ndim_rank_znl = ndim_rank_znl + 1
        END IF
      END DO
      IF (ALLOCATED(nrank_znl)) DEALLOCATE(nrank_znl)
      ALLOCATE(nrank_znl(ndim_rank_znl))
      !$ACC KERNELS
      ii = 0
      nrank_znl(:) = 0
      !$ACC END KERNELS
      DO jproc = 1, jpnij
        IF (kwork(jproc) == njmpp) THEN
          ii = ii + 1
          nrank_znl(ii) = jproc - 1
        END IF
      END DO
      CALL MPI_COMM_GROUP(mpi_comm_oce, ngrp_opa, ierr)
      CALL MPI_GROUP_INCL(ngrp_opa, ndim_rank_znl, nrank_znl, ngrp_znl, ierr)
      CALL MPI_COMM_CREATE(mpi_comm_oce, ngrp_znl, ncomm_znl, ierr)
    END IF
    IF (jpni == 1) THEN
      l_znl_root = .TRUE.
    ELSE
      l_znl_root = .FALSE.
      kwork(1) = nimpp
      CALL mpp_min(kwork(1), kcom = ncomm_znl)
      IF (nimpp == kwork(1)) l_znl_root = .TRUE.
    END IF
    DEALLOCATE(kwork)
  END SUBROUTINE mpp_ini_znl
  SUBROUTINE mpp_ini_north
    INTEGER :: ierr
    INTEGER :: jjproc
    INTEGER :: ii, ji
    njmppmax = MAXVAL(njmppt)
    ndim_rank_north = 0
    DO jjproc = 1, jpnij
      IF (njmppt(jjproc) == njmppmax) ndim_rank_north = ndim_rank_north + 1
    END DO
    IF (ALLOCATED(nrank_north)) DEALLOCATE(nrank_north)
    ALLOCATE(nrank_north(ndim_rank_north))
    ii = 0
    DO ji = 1, jpnij
      IF (njmppt(ji) == njmppmax) THEN
        ii = ii + 1
        nrank_north(ii) = ji - 1
      END IF
    END DO
    CALL MPI_COMM_GROUP(mpi_comm_oce, ngrp_world, ierr)
    CALL MPI_GROUP_INCL(ngrp_world, ndim_rank_north, nrank_north, ngrp_north, ierr)
    CALL MPI_COMM_CREATE(mpi_comm_oce, ngrp_north, ncomm_north, ierr)
  END SUBROUTINE mpp_ini_north
  SUBROUTINE mpi_init_oce(ldtxt, ksft, code)
    CHARACTER(LEN = *), DIMENSION(:), INTENT(OUT) :: ldtxt
    INTEGER, INTENT(INOUT) :: ksft
    INTEGER, INTENT(OUT) :: code
    INTEGER :: ierr, ji
    LOGICAL :: mpi_was_called
    CALL mpi_initialized(mpi_was_called, code)
    IF (code /= MPI_SUCCESS) THEN
      DO ji = 1, SIZE(ldtxt)
        IF (TRIM(ldtxt(ji)) /= '') WRITE(*, FMT = *) ldtxt(ji)
      END DO
      WRITE(*, cform_err)
      WRITE(*, FMT = *) ' lib_mpp: Error in routine mpi_initialized'
      CALL mpi_abort(mpi_comm_world, code, ierr)
    END IF
    IF (.NOT. mpi_was_called) THEN
      CALL mpi_init(code)
      CALL mpi_comm_dup(mpi_comm_world, mpi_comm_oce, code)
      IF (code /= MPI_SUCCESS) THEN
        DO ji = 1, SIZE(ldtxt)
          IF (TRIM(ldtxt(ji)) /= '') WRITE(*, FMT = *) ldtxt(ji)
        END DO
        WRITE(*, cform_err)
        WRITE(*, FMT = *) ' lib_mpp: Error in routine mpi_comm_dup'
        CALL mpi_abort(mpi_comm_world, code, ierr)
      END IF
    END IF
    IF (nn_buffer > 0) THEN
      WRITE(ldtxt(ksft), FMT = *) 'mpi_bsend, buffer allocation of  : ', nn_buffer
      ksft = ksft + 1
      ALLOCATE(tampon(nn_buffer), STAT = ierr)
      IF (ierr /= 0) THEN
        DO ji = 1, SIZE(ldtxt)
          IF (TRIM(ldtxt(ji)) /= '') WRITE(*, FMT = *) ldtxt(ji)
        END DO
        WRITE(*, cform_err)
        WRITE(*, FMT = *) ' lib_mpp: Error in ALLOCATE', ierr
        CALL mpi_abort(mpi_comm_world, code, ierr)
      END IF
      CALL mpi_buffer_attach(tampon, nn_buffer, code)
    END IF
  END SUBROUTINE mpi_init_oce
  SUBROUTINE DDPDD_MPI(ydda, yddb, ilen, itype)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(IN) :: ilen, itype
    COMPLEX(KIND = wp), DIMENSION(ilen), INTENT(IN) :: ydda
    COMPLEX(KIND = wp), DIMENSION(ilen), INTENT(INOUT) :: yddb
    REAL(KIND = wp) :: zerr, zt1, zt2
    INTEGER :: ji, ztmp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ddpdd_mpi', 'r0', 0, 0)
    ztmp = itype
    DO ji = 1, ilen
      zt1 = REAL(ydda(ji)) + REAL(yddb(ji))
      zerr = zt1 - REAL(ydda(ji))
      zt2 = ((REAL(yddb(ji)) - zerr) + (REAL(ydda(ji)) - (zt1 - zerr))) + AIMAG(ydda(ji)) + AIMAG(yddb(ji))
      yddb(ji) = CMPLX(zt1 + zt2, zt2 - ((zt1 + zt2) - zt1), wp)
    END DO
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE DDPDD_MPI
  SUBROUTINE mpp_lbc_north_icb(pt2d, cd_type, psgn, kextj)
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pt2d
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: kextj
    INTEGER :: ji, jj, jr
    INTEGER :: ierr, itaille, ildi, ilei, iilb
    INTEGER :: ipj, ij, iproc
    REAL(KIND = wp), DIMENSION(:, :), ALLOCATABLE :: ztab_e, znorthloc_e
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: znorthgloio_e
    ipj = 4
    ALLOCATE(ztab_e(jpiglo, 1 - kextj : ipj + kextj), znorthloc_e(jpimax, 1 - kextj : ipj + kextj), znorthgloio_e(jpimax, 1 - &
&kextj : ipj + kextj, jpni))
    !$ACC KERNELS
    ztab_e(:, :) = 0._wp
    znorthloc_e(:, :) = 0._wp
    ij = 1 - kextj
    DO jj = jpj - ipj + 1 - kextj, jpj + kextj
      znorthloc_e(1 : jpi, ij) = pt2d(1 : jpi, jj)
      ij = ij + 1
    END DO
    itaille = jpimax * (ipj + 2 * kextj)
    !$ACC END KERNELS
    CALL MPI_ALLGATHER(znorthloc_e(1, 1 - kextj), itaille, MPI_DOUBLE_PRECISION, znorthgloio_e(1, 1 - kextj, 1), itaille, &
&MPI_DOUBLE_PRECISION, ncomm_north, ierr)
    !$ACC KERNELS
    DO jr = 1, ndim_rank_north
      iproc = nrank_north(jr) + 1
      ildi = nldit(iproc)
      ilei = nleit(iproc)
      iilb = nimppt(iproc)
      !$ACC LOOP INDEPENDENT COLLAPSE(2)
      DO jj = 1 - kextj, ipj + kextj
        DO ji = ildi, ilei
          ztab_e(ji + iilb - 1, jj) = znorthgloio_e(ji, jj, jr)
        END DO
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_nfd(ztab_e(:, 1 - kextj : ipj + kextj), cd_type, psgn, kextj)
    !$ACC KERNELS
    ij = 1 - kextj
    DO jj = jpj - ipj + 1 - kextj, jpj + kextj
      DO ji = 1, jpi
        pt2d(ji, jj) = ztab_e(ji + nimpp - 1, ij)
      END DO
      ij = ij + 1
    END DO
    !$ACC END KERNELS
    DEALLOCATE(ztab_e, znorthloc_e, znorthgloio_e)
  END SUBROUTINE mpp_lbc_north_icb
  SUBROUTINE mpp_lnk_2d_icb(pt2d, cd_type, psgn, kexti, kextj)
    REAL(KIND = wp), DIMENSION(1 - kexti : jpi + kexti, 1 - kextj : jpj + kextj), INTENT(INOUT) :: pt2d
    CHARACTER(LEN = 1), INTENT(IN) :: cd_type
    REAL(KIND = wp), INTENT(IN) :: psgn
    INTEGER, INTENT(IN) :: kexti
    INTEGER, INTENT(IN) :: kextj
    INTEGER :: jl
    INTEGER :: imigr, iihom, ijhom
    INTEGER :: ipreci, iprecj
    INTEGER :: ml_req1, ml_req2, ml_err
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ml_stat
    REAL(KIND = wp), DIMENSION(1 - kexti : jpi + kexti, nn_hls + kextj, 2) :: r2dns, r2dsn
    REAL(KIND = wp), DIMENSION(1 - kextj : jpj + kextj, nn_hls + kexti, 2) :: r2dwe, r2dew
    ipreci = nn_hls + kexti
    iprecj = nn_hls + kextj
    IF (l_Iperio) THEN
      !$ACC KERNELS
      pt2d(1 - kexti : 1, :) = pt2d(jpim1 - kexti : jpim1, :)
      pt2d(jpi : jpi + kexti, :) = pt2d(2 : 2 + kexti, :)
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      IF (.NOT. cd_type == 'F') pt2d(1 - kexti : nn_hls, :) = 0._wp
      pt2d(jpi - nn_hls + 1 : jpi + kexti, :) = 0._wp
      !$ACC END KERNELS
    END IF
    IF (l_Jperio) THEN
      !$ACC KERNELS
      pt2d(:, 1 - kextj : 1) = pt2d(:, jpjm1 - kextj : jpjm1)
      pt2d(:, jpj : jpj + kextj) = pt2d(:, 2 : 2 + kextj)
      !$ACC END KERNELS
    ELSE
      !$ACC KERNELS
      IF (.NOT. cd_type == 'F') pt2d(:, 1 - kextj : nn_hls) = 0._wp
      pt2d(:, jpj - nn_hls + 1 : jpj + kextj) = 0._wp
      !$ACC END KERNELS
    END IF
    IF (npolj /= 0) THEN
      SELECT CASE (jpni)
      CASE (1)
        CALL lbc_nfd(pt2d(1 : jpi, 1 : jpj + kextj), cd_type, psgn, kextj)
      CASE DEFAULT
        CALL mpp_lbc_north_icb(pt2d(1 : jpi, 1 : jpj + kextj), cd_type, psgn, kextj)
      END SELECT
    END IF
    !$ACC KERNELS
    SELECT CASE (nbondi)
    CASE (- 1, 0, 1)
      iihom = jpi - nreci - kexti
      DO jl = 1, ipreci
        r2dew(:, jl, 1) = pt2d(nn_hls + jl, :)
        r2dwe(:, jl, 1) = pt2d(iihom + jl, :)
      END DO
    END SELECT
    imigr = ipreci * (jpj + 2 * kextj)
    !$ACC END KERNELS
    SELECT CASE (nbondi)
    CASE (- 1)
      CALL mppsend(2, r2dwe(1 - kextj, 1, 1), imigr, noea, ml_req1)
      CALL mpprecv(1, r2dew(1 - kextj, 1, 2), imigr, noea)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(1, r2dew(1 - kextj, 1, 1), imigr, nowe, ml_req1)
      CALL mppsend(2, r2dwe(1 - kextj, 1, 1), imigr, noea, ml_req2)
      CALL mpprecv(1, r2dew(1 - kextj, 1, 2), imigr, noea)
      CALL mpprecv(2, r2dwe(1 - kextj, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(1, r2dew(1 - kextj, 1, 1), imigr, nowe, ml_req1)
      CALL mpprecv(2, r2dwe(1 - kextj, 1, 2), imigr, nowe)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    !$ACC KERNELS
    iihom = jpi - nn_hls
    SELECT CASE (nbondi)
    CASE (- 1)
      DO jl = 1, ipreci
        pt2d(iihom + jl, :) = r2dew(:, jl, 2)
      END DO
    CASE (0)
      DO jl = 1, ipreci
        pt2d(jl - kexti, :) = r2dwe(:, jl, 2)
        pt2d(iihom + jl, :) = r2dew(:, jl, 2)
      END DO
    CASE (1)
      DO jl = 1, ipreci
        pt2d(jl - kexti, :) = r2dwe(:, jl, 2)
      END DO
    END SELECT
    IF (nbondj /= 2) THEN
      ijhom = jpj - nrecj - kextj
      DO jl = 1, iprecj
        r2dsn(:, jl, 1) = pt2d(:, ijhom + jl)
        r2dns(:, jl, 1) = pt2d(:, nn_hls + jl)
      END DO
    END IF
    imigr = iprecj * (jpi + 2 * kexti)
    !$ACC END KERNELS
    SELECT CASE (nbondj)
    CASE (- 1)
      CALL mppsend(4, r2dsn(1 - kexti, 1, 1), imigr, nono, ml_req1)
      CALL mpprecv(3, r2dns(1 - kexti, 1, 2), imigr, nono)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    CASE (0)
      CALL mppsend(3, r2dns(1 - kexti, 1, 1), imigr, noso, ml_req1)
      CALL mppsend(4, r2dsn(1 - kexti, 1, 1), imigr, nono, ml_req2)
      CALL mpprecv(3, r2dns(1 - kexti, 1, 2), imigr, nono)
      CALL mpprecv(4, r2dsn(1 - kexti, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      IF (l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
    CASE (1)
      CALL mppsend(3, r2dns(1 - kexti, 1, 1), imigr, noso, ml_req1)
      CALL mpprecv(4, r2dsn(1 - kexti, 1, 2), imigr, noso)
      IF (l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
    END SELECT
    !$ACC KERNELS
    ijhom = jpj - nn_hls
    SELECT CASE (nbondj)
    CASE (- 1)
      DO jl = 1, iprecj
        pt2d(:, ijhom + jl) = r2dns(:, jl, 2)
      END DO
    CASE (0)
      DO jl = 1, iprecj
        pt2d(:, jl - kextj) = r2dsn(:, jl, 2)
        pt2d(:, ijhom + jl) = r2dns(:, jl, 2)
      END DO
    CASE (1)
      DO jl = 1, iprecj
        pt2d(:, jl - kextj) = r2dsn(:, jl, 2)
      END DO
    END SELECT
    !$ACC END KERNELS
  END SUBROUTINE mpp_lnk_2d_icb
  SUBROUTINE ctl_stop(cd1, cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd1, cd2, cd3, cd4, cd5
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd6, cd7, cd8, cd9, cd10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_stop', 'r0', 0, 0)
    nstop = nstop + 1
    IF (lwp) THEN
      WRITE(numout, cform_err)
      IF (PRESENT(cd1)) WRITE(numout, FMT = *) cd1
      IF (PRESENT(cd2)) WRITE(numout, FMT = *) cd2
      IF (PRESENT(cd3)) WRITE(numout, FMT = *) cd3
      IF (PRESENT(cd4)) WRITE(numout, FMT = *) cd4
      IF (PRESENT(cd5)) WRITE(numout, FMT = *) cd5
      IF (PRESENT(cd6)) WRITE(numout, FMT = *) cd6
      IF (PRESENT(cd7)) WRITE(numout, FMT = *) cd7
      IF (PRESENT(cd8)) WRITE(numout, FMT = *) cd8
      IF (PRESENT(cd9)) WRITE(numout, FMT = *) cd9
      IF (PRESENT(cd10)) WRITE(numout, FMT = *) cd10
    END IF
    CALL FLUSH(numout)
    IF (numstp /= - 1) CALL FLUSH(numstp)
    IF (numrun /= - 1) CALL FLUSH(numrun)
    IF (numevo_ice /= - 1) CALL FLUSH(numevo_ice)
    IF (cd1 == 'STOP') THEN
      IF (lwp) WRITE(numout, FMT = *) 'huge E-R-R-O-R : immediate stop'
      CALL mppstop
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ctl_stop
  SUBROUTINE ctl_warn(cd1, cd2, cd3, cd4, cd5, cd6, cd7, cd8, cd9, cd10)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd1, cd2, cd3, cd4, cd5
    CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: cd6, cd7, cd8, cd9, cd10
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_warn', 'r0', 0, 0)
    nwarn = nwarn + 1
    IF (lwp) THEN
      WRITE(numout, cform_war)
      IF (PRESENT(cd1)) WRITE(numout, FMT = *) cd1
      IF (PRESENT(cd2)) WRITE(numout, FMT = *) cd2
      IF (PRESENT(cd3)) WRITE(numout, FMT = *) cd3
      IF (PRESENT(cd4)) WRITE(numout, FMT = *) cd4
      IF (PRESENT(cd5)) WRITE(numout, FMT = *) cd5
      IF (PRESENT(cd6)) WRITE(numout, FMT = *) cd6
      IF (PRESENT(cd7)) WRITE(numout, FMT = *) cd7
      IF (PRESENT(cd8)) WRITE(numout, FMT = *) cd8
      IF (PRESENT(cd9)) WRITE(numout, FMT = *) cd9
      IF (PRESENT(cd10)) WRITE(numout, FMT = *) cd10
    END IF
    CALL FLUSH(numout)
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ctl_warn
  SUBROUTINE ctl_opn(knum, cdfile, cdstat, cdform, cdacce, klengh, kout, ldwp, karea)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(OUT) :: knum
    CHARACTER(LEN = *), INTENT(IN) :: cdfile
    CHARACTER(LEN = *), INTENT(IN) :: cdstat
    CHARACTER(LEN = *), INTENT(IN) :: cdform
    CHARACTER(LEN = *), INTENT(IN) :: cdacce
    INTEGER, INTENT(IN) :: klengh
    INTEGER, INTENT(IN) :: kout
    LOGICAL, INTENT(IN) :: ldwp
    INTEGER, OPTIONAL, INTENT(IN) :: karea
    CHARACTER(LEN = 80) :: clfile
    INTEGER :: iost
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_opn', 'r0', 0, 0)
    clfile = TRIM(cdfile)
    IF (PRESENT(karea)) THEN
      IF (karea > 1) WRITE(clfile, FMT = "(a,'_',i4.4)") TRIM(clfile), karea - 1
    END IF
    knum = get_unit()
    iost = 0
    IF (cdacce(1 : 6) == 'DIRECT') THEN
      OPEN(UNIT = knum, FILE = clfile, FORM = cdform, ACCESS = cdacce, STATUS = cdstat, RECL = klengh, ERR = 100, IOSTAT = iost)
    ELSE
      OPEN(UNIT = knum, FILE = clfile, FORM = cdform, ACCESS = cdacce, STATUS = cdstat, ERR = 100, IOSTAT = iost)
    END IF
    IF (iost == 0) THEN
      IF (ldwp) THEN
        WRITE(kout, FMT = *) '     file   : ', clfile, ' open ok'
        WRITE(kout, FMT = *) '     unit   = ', knum
        WRITE(kout, FMT = *) '     status = ', cdstat
        WRITE(kout, FMT = *) '     form   = ', cdform
        WRITE(kout, FMT = *) '     access = ', cdacce
        WRITE(kout, FMT = *)
      END IF
    END IF
100 CONTINUE
    IF (iost /= 0) THEN
      IF (ldwp) THEN
        WRITE(kout, FMT = *)
        WRITE(kout, FMT = *) ' ===>>>> : bad opening file: ', clfile
        WRITE(kout, FMT = *) ' =======   ===  '
        WRITE(kout, FMT = *) '           unit   = ', knum
        WRITE(kout, FMT = *) '           status = ', cdstat
        WRITE(kout, FMT = *) '           form   = ', cdform
        WRITE(kout, FMT = *) '           access = ', cdacce
        WRITE(kout, FMT = *) '           iostat = ', iost
        WRITE(kout, FMT = *) '           we stop. verify the file '
        WRITE(kout, FMT = *)
      ELSE
        WRITE(*, FMT = *)
        WRITE(*, FMT = *) ' ===>>>> : bad opening file: ', clfile
        WRITE(*, FMT = *) ' =======   ===  '
        WRITE(*, FMT = *) '           unit   = ', knum
        WRITE(*, FMT = *) '           status = ', cdstat
        WRITE(*, FMT = *) '           form   = ', cdform
        WRITE(*, FMT = *) '           access = ', cdacce
        WRITE(*, FMT = *) '           iostat = ', iost
        WRITE(*, FMT = *) '           we stop. verify the file '
        WRITE(*, FMT = *)
      END IF
      CALL FLUSH(kout)
      STOP 'ctl_opn bad opening'
    END IF
    CALL profile_psy_data0 % PostEnd
  END SUBROUTINE ctl_opn
  SUBROUTINE ctl_nam(kios, cdnam, ldwp)
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    INTEGER, INTENT(INOUT) :: kios
    CHARACTER(LEN = *), INTENT(IN) :: cdnam
    CHARACTER(LEN = 5) :: clios
    LOGICAL, INTENT(IN) :: ldwp
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('ctl_nam', 'r0', 0, 0)
    WRITE(clios, FMT = '(I5.0)') kios
    IF (kios < 0) THEN
      CALL ctl_warn('end of record or file while reading namelist ' // TRIM(cdnam) // ' iostat = ' // TRIM(clios))
    END IF
    IF (kios > 0) THEN
      CALL ctl_stop('misspelled variable in namelist ' // TRIM(cdnam) // ' iostat = ' // TRIM(clios))
    END IF
    kios = 0
    CALL profile_psy_data0 % PostEnd
    RETURN
  END SUBROUTINE ctl_nam
  INTEGER FUNCTION get_unit()
    USE profile_psy_data_mod, ONLY: profile_PSyDataType
    LOGICAL :: llopn
    TYPE(profile_PSyDataType), TARGET, SAVE :: profile_psy_data0
    CALL profile_psy_data0 % PreStart('get_unit', 'r0', 0, 0)
    get_unit = 15
    llopn = .TRUE.
    DO WHILE ((get_unit < 998) .AND. llopn)
      get_unit = get_unit + 1
      INQUIRE(UNIT = get_unit, OPENED = llopn)
    END DO
    IF ((get_unit == 999) .AND. llopn) THEN
      CALL ctl_stop('get_unit: All logical units until 999 are used...')
      get_unit = - 1
    END IF
    CALL profile_psy_data0 % PostEnd
  END FUNCTION get_unit
END MODULE lib_mpp
