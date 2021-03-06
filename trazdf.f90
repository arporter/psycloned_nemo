MODULE trazdf
  USE oce
  USE dom_oce
  USE domvvl
  USE phycst
  USE zdf_oce
  USE sbc_oce
  USE ldftra
  USE ldfslp
  USE trd_oce
  USE trdtra
  USE in_out_manager
  USE prtctl
  USE lbclnk
  USE lib_mpp
  USE timing
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: tra_zdf
  PUBLIC :: tra_zdf_imp
  CONTAINS
  SUBROUTINE tra_zdf(kt)
    INTEGER, INTENT(IN) :: kt
    INTEGER :: jk
    REAL(KIND = wp), DIMENSION(:, :, :), ALLOCATABLE :: ztrdt, ztrds
    IF (ln_timing) CALL timing_start('tra_zdf')
    IF (kt == nit000) THEN
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'tra_zdf : implicit vertical mixing on T & S'
      IF (lwp) WRITE(numout, FMT = *) '~~~~~~~ '
    END IF
    IF (neuler == 0 .AND. kt == nit000) THEN
      r2dt = rdt
    ELSE IF (kt <= nit000 + 1) THEN
      r2dt = 2. * rdt
    END IF
    IF (l_trdtra) THEN
      ALLOCATE(ztrdt(jpi, jpj, jpk), ztrds(jpi, jpj, jpk))
      !$ACC KERNELS
      ztrdt(:, :, :) = tsa(:, :, :, jp_tem)
      ztrds(:, :, :) = tsa(:, :, :, jp_sal)
      !$ACC END KERNELS
    END IF
    CALL tra_zdf_imp(kt, nit000, 'TRA', r2dt, tsb, tsa, jpts)
    WHERE (tsa(:, :, :, jp_sal) < 0._wp) tsa(:, :, :, jp_sal) = 0.1_wp
    IF (l_trdtra) THEN
      !$ACC KERNELS
      DO jk = 1, jpkm1
        ztrdt(:, :, jk) = ((tsa(:, :, jk, jp_tem) * e3t_a(:, :, jk) - tsb(:, :, jk, jp_tem) * e3t_b(:, :, jk)) / (e3t_n(:, :, jk) * r2dt)) - ztrdt(:, :, jk)
        ztrds(:, :, jk) = ((tsa(:, :, jk, jp_sal) * e3t_a(:, :, jk) - tsb(:, :, jk, jp_sal) * e3t_b(:, :, jk)) / (e3t_n(:, :, jk) * r2dt)) - ztrds(:, :, jk)
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk_multi(ztrdt, 'T', 1., ztrds, 'T', 1.)
      CALL trd_tra(kt, 'TRA', jp_tem, jptra_zdf, ztrdt)
      CALL trd_tra(kt, 'TRA', jp_sal, jptra_zdf, ztrds)
      DEALLOCATE(ztrdt, ztrds)
    END IF
    IF (ln_ctl) CALL prt_ctl(tab3d_1 = tsa(:, :, :, jp_tem), clinfo1 = ' zdf  - Ta: ', mask1 = tmask, tab3d_2 = tsa(:, :, :, jp_sal), clinfo2 = ' Sa: ', mask2 = tmask, clinfo3 = 'tra')
    IF (ln_timing) CALL timing_stop('tra_zdf')
  END SUBROUTINE tra_zdf
  SUBROUTINE tra_zdf_imp(kt, kit000, cdtype, p2dt, ptb, pta, kjpt)
    INTEGER, INTENT(IN   ) :: kt
    INTEGER, INTENT(IN   ) :: kit000
    CHARACTER(LEN = 3), INTENT(IN   ) :: cdtype
    INTEGER, INTENT(IN   ) :: kjpt
    REAL(KIND = wp), INTENT(IN   ) :: p2dt
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(IN   ) :: ptb
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk, kjpt), INTENT(INOUT) :: pta
    INTEGER :: ji, jj, jk, jn
    REAL(KIND = wp) :: zrhs
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zwi, zwt, zwd, zws
    DO jn = 1, kjpt
      IF ((cdtype == 'TRA' .AND. (jn == jp_tem .OR. (jn == jp_sal .AND. ln_zdfddm))) .OR. (cdtype == 'TRC' .AND. jn == 1)) THEN
        IF (cdtype == 'TRA' .AND. jn == jp_tem) THEN
          !$ACC KERNELS
          zwt(:, :, 2 : jpk) = avt(:, :, 2 : jpk)
          !$ACC END KERNELS
        ELSE
          !$ACC KERNELS
          zwt(:, :, 2 : jpk) = avs(:, :, 2 : jpk)
          !$ACC END KERNELS
        END IF
        !$ACC KERNELS
        zwt(:, :, 1) = 0._wp
        !$ACC END KERNELS
        IF (l_ldfslp) THEN
          IF (ln_traldf_msc) THEN
            !$ACC KERNELS
            DO jk = 2, jpkm1
              DO jj = 2, jpjm1
                DO ji = 2, jpim1
                  zwt(ji, jj, jk) = zwt(ji, jj, jk) + akz(ji, jj, jk)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          ELSE
            !$ACC KERNELS
            DO jk = 2, jpkm1
              DO jj = 2, jpjm1
                DO ji = 2, jpim1
                  zwt(ji, jj, jk) = zwt(ji, jj, jk) + ah_wslp2(ji, jj, jk)
                END DO
              END DO
            END DO
            !$ACC END KERNELS
          END IF
        END IF
        !$ACC KERNELS
        DO jk = 1, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zwi(ji, jj, jk) = - p2dt * zwt(ji, jj, jk) / e3w_n(ji, jj, jk)
              zws(ji, jj, jk) = - p2dt * zwt(ji, jj, jk + 1) / e3w_n(ji, jj, jk + 1)
              zwd(ji, jj, jk) = e3t_a(ji, jj, jk) - zwi(ji, jj, jk) - zws(ji, jj, jk)
            END DO
          END DO
        END DO
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zwt(ji, jj, 1) = zwd(ji, jj, 1)
          END DO
        END DO
        DO jk = 2, jpkm1
          DO jj = 2, jpjm1
            DO ji = 2, jpim1
              zwt(ji, jj, jk) = zwd(ji, jj, jk) - zwi(ji, jj, jk) * zws(ji, jj, jk - 1) / zwt(ji, jj, jk - 1)
            END DO
          END DO
        END DO
        !$ACC END KERNELS
      END IF
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          pta(ji, jj, 1, jn) = e3t_b(ji, jj, 1) * ptb(ji, jj, 1, jn) + p2dt * e3t_n(ji, jj, 1) * pta(ji, jj, 1, jn)
        END DO
      END DO
      DO jk = 2, jpkm1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            zrhs = e3t_b(ji, jj, jk) * ptb(ji, jj, jk, jn) + p2dt * e3t_n(ji, jj, jk) * pta(ji, jj, jk, jn)
            pta(ji, jj, jk, jn) = zrhs - zwi(ji, jj, jk) / zwt(ji, jj, jk - 1) * pta(ji, jj, jk - 1, jn)
          END DO
        END DO
      END DO
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          pta(ji, jj, jpkm1, jn) = pta(ji, jj, jpkm1, jn) / zwt(ji, jj, jpkm1) * tmask(ji, jj, jpkm1)
        END DO
      END DO
      DO jk = jpk - 2, 1, - 1
        DO jj = 2, jpjm1
          DO ji = 2, jpim1
            pta(ji, jj, jk, jn) = (pta(ji, jj, jk, jn) - zws(ji, jj, jk) * pta(ji, jj, jk + 1, jn)) / zwt(ji, jj, jk) * tmask(ji, jj, jk)
          END DO
        END DO
      END DO
      !$ACC END KERNELS
    END DO
  END SUBROUTINE tra_zdf_imp
END MODULE trazdf