MODULE icedyn_adv_umx
  USE phycst
  USE dom_oce
  USE sbc_oce, ONLY: nn_fsbc
  USE ice
  USE in_out_manager
  USE lib_mpp
  USE lib_fortran
  USE lbclnk
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ice_dyn_adv_umx
  REAL(KIND = wp) :: z1_6 = 1._wp / 6._wp
  REAL(KIND = wp) :: z1_120 = 1._wp / 120._wp
  CONTAINS
  SUBROUTINE ice_dyn_adv_umx(k_order, kt, pu_ice, pv_ice, pato_i, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pe_s, pe_i)
    INTEGER, INTENT(IN   ) :: k_order
    INTEGER, INTENT(IN   ) :: kt
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN   ) :: pu_ice
    REAL(KIND = wp), DIMENSION(:, :), INTENT(IN   ) :: pv_ice
    REAL(KIND = wp), DIMENSION(:, :), INTENT(INOUT) :: pato_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pv_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pv_s
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: psv_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: poa_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pa_i
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pa_ip
    REAL(KIND = wp), DIMENSION(:, :, :), INTENT(INOUT) :: pv_ip
    REAL(KIND = wp), DIMENSION(:, :, :, :), INTENT(INOUT) :: pe_s
    REAL(KIND = wp), DIMENSION(:, :, :, :), INTENT(INOUT) :: pe_i
    INTEGER :: ji, jj, jk, jl, jt
    INTEGER :: initad
    REAL(KIND = wp) :: zcfl, zusnit, zdt
    REAL(KIND = wp), ALLOCATABLE, DIMENSION(:, :) :: zudy, zvdx, zcu_box, zcv_box
    IF (kt == nit000 .AND. lwp) WRITE(numout, FMT = *) '-- ice_dyn_adv_umx: Ultimate-Macho advection scheme'
    ALLOCATE(zudy(jpi, jpj), zvdx(jpi, jpj), zcu_box(jpi, jpj), zcv_box(jpi, jpj))
    zcfl = MAXVAL(ABS(pu_ice(:, :)) * rdt_ice * r1_e1u(:, :))
    zcfl = MAX(zcfl, MAXVAL(ABS(pv_ice(:, :)) * rdt_ice * r1_e2v(:, :)))
    IF (lk_mpp) CALL mpp_max(zcfl)
    IF (zcfl > 0.5) THEN
      initad = 2
      zusnit = 0.5_wp
    ELSE
      initad = 1
      zusnit = 1.0_wp
    END IF
    !$ACC KERNELS
    zdt = rdt_ice / REAL(initad)
    zudy(:, :) = pu_ice(:, :) * e2u(:, :)
    zvdx(:, :) = pv_ice(:, :) * e1v(:, :)
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        IF (pu_ice(ji, jj) * pu_ice(ji - 1, jj) <= 0._wp) THEN
          zcu_box(ji, jj) = 0._wp
        ELSE IF (pu_ice(ji, jj) > 0._wp) THEN
          zcu_box(ji, jj) = pu_ice(ji - 1, jj)
        ELSE
          zcu_box(ji, jj) = pu_ice(ji, jj)
        END IF
        IF (pv_ice(ji, jj) * pv_ice(ji, jj - 1) <= 0._wp) THEN
          zcv_box(ji, jj) = 0._wp
        ELSE IF (pv_ice(ji, jj) > 0._wp) THEN
          zcv_box(ji, jj) = pv_ice(ji, jj - 1)
        ELSE
          zcv_box(ji, jj) = pv_ice(ji, jj)
        END IF
      END DO
    END DO
    !$ACC END KERNELS
    DO jt = 1, initad
      CALL adv_umx(k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pato_i(:, :))
      DO jl = 1, jpl
        CALL adv_umx(k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pa_i(:, :, jl))
        CALL adv_umx(k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pv_i(:, :, jl))
        CALL adv_umx(k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, psv_i(:, :, jl))
        CALL adv_umx(k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, poa_i(:, :, jl))
        DO jk = 1, nlay_i
          CALL adv_umx(k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pe_i(:, :, jk, jl))
        END DO
        CALL adv_umx(k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pv_s(:, :, jl))
        DO jk = 1, nlay_s
          CALL adv_umx(k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pe_s(:, :, jk, jl))
        END DO
        IF (ln_pnd_H12) THEN
          CALL adv_umx(k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pa_ip(:, :, jl))
          CALL adv_umx(k_order, kt, zdt, zudy, zvdx, zcu_box, zcv_box, pv_ip(:, :, jl))
        END IF
      END DO
    END DO
    DEALLOCATE(zudy, zvdx, zcu_box, zcv_box)
  END SUBROUTINE ice_dyn_adv_umx
  SUBROUTINE adv_umx(k_order, kt, pdt, puc, pvc, pubox, pvbox, ptc)
    INTEGER, INTENT(IN   ) :: k_order
    INTEGER, INTENT(IN   ) :: kt
    REAL(KIND = wp), INTENT(IN   ) :: pdt
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: puc, pvc
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: pubox, pvbox
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: ptc
    INTEGER :: ji, jj
    REAL(KIND = wp) :: ztra
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zfu_ups, zfu_ho, zt_u, zt_ups
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zfv_ups, zfv_ho, zt_v, ztrd
    !$ACC KERNELS
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        zfu_ups(ji, jj) = MAX(puc(ji, jj), 0._wp) * ptc(ji, jj) + MIN(puc(ji, jj), 0._wp) * ptc(ji + 1, jj)
        zfv_ups(ji, jj) = MAX(pvc(ji, jj), 0._wp) * ptc(ji, jj) + MIN(pvc(ji, jj), 0._wp) * ptc(ji, jj + 1)
      END DO
    END DO
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        ztra = - (zfu_ups(ji, jj) - zfu_ups(ji - 1, jj) + zfv_ups(ji, jj) - zfv_ups(ji, jj - 1)) * r1_e1e2t(ji, jj)
        ztrd(ji, jj) = ztra
        zt_ups(ji, jj) = (ptc(ji, jj) + pdt * ztra) * tmask(ji, jj, 1)
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk(zt_ups, 'T', 1.)
    SELECT CASE (k_order)
    CASE (20)
      !$ACC KERNELS
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zfu_ho(ji, jj) = 0.5 * puc(ji, jj) * (ptc(ji, jj) + ptc(ji + 1, jj))
          zfv_ho(ji, jj) = 0.5 * pvc(ji, jj) * (ptc(ji, jj) + ptc(ji, jj + 1))
        END DO
      END DO
      !$ACC END KERNELS
    CASE (1 : 5)
      CALL macho(k_order, kt, pdt, ptc, puc, pvc, pubox, pvbox, zt_u, zt_v)
      !$ACC KERNELS
      DO jj = 1, jpjm1
        DO ji = 1, jpim1
          zfu_ho(ji, jj) = puc(ji, jj) * zt_u(ji, jj)
          zfv_ho(ji, jj) = pvc(ji, jj) * zt_v(ji, jj)
        END DO
      END DO
      !$ACC END KERNELS
    END SELECT
    !$ACC KERNELS
    DO jj = 1, jpjm1
      DO ji = 1, jpim1
        zfu_ho(ji, jj) = zfu_ho(ji, jj) - zfu_ups(ji, jj)
        zfv_ho(ji, jj) = zfv_ho(ji, jj) - zfv_ups(ji, jj)
      END DO
    END DO
    !$ACC END KERNELS
    CALL nonosc_2d(ptc, zfu_ho, zfv_ho, zt_ups, pdt)
    !$ACC KERNELS
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        ztra = ztrd(ji, jj) - (zfu_ho(ji, jj) - zfu_ho(ji - 1, jj) + zfv_ho(ji, jj) - zfv_ho(ji, jj - 1)) * r1_e1e2t(ji, jj)
        ptc(ji, jj) = (ptc(ji, jj) + pdt * ztra) * tmask(ji, jj, 1)
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk(ptc, 'T', 1.)
  END SUBROUTINE adv_umx
  SUBROUTINE macho(k_order, kt, pdt, ptc, puc, pvc, pubox, pvbox, pt_u, pt_v)
    INTEGER, INTENT(IN   ) :: k_order
    INTEGER, INTENT(IN   ) :: kt
    REAL(KIND = wp), INTENT(IN   ) :: pdt
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: ptc
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: puc, pvc
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: pubox, pvbox
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(  OUT) :: pt_u, pt_v
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zc_box
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zzt
    IF (MOD((kt - 1) / nn_fsbc, 2) == 0) THEN
      CALL ultimate_x(k_order, pdt, ptc, puc, pt_u)
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zzt(ji, jj) = ptc(ji, jj) - pubox(ji, jj) * pdt * (pt_u(ji, jj) - pt_u(ji - 1, jj)) * r1_e1t(ji, jj) - ptc(ji, jj) * pdt * (puc(ji, jj) - puc(ji - 1, jj)) * r1_e1e2t(ji, jj)
          zzt(ji, jj) = zzt(ji, jj) * tmask(ji, jj, 1)
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk(zzt, 'T', 1.)
      CALL ultimate_y(k_order, pdt, zzt, pvc, pt_v)
    ELSE
      CALL ultimate_y(k_order, pdt, ptc, pvc, pt_v)
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = 2, jpim1
          zzt(ji, jj) = ptc(ji, jj) - pvbox(ji, jj) * pdt * (pt_v(ji, jj) - pt_v(ji, jj - 1)) * r1_e2t(ji, jj) - ptc(ji, jj) * pdt * (pvc(ji, jj) - pvc(ji, jj - 1)) * r1_e1e2t(ji, jj)
          zzt(ji, jj) = zzt(ji, jj) * tmask(ji, jj, 1)
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk(zzt, 'T', 1.)
      CALL ultimate_x(k_order, pdt, zzt, puc, pt_u)
    END IF
  END SUBROUTINE macho
  SUBROUTINE ultimate_x(k_order, pdt, pt, puc, pt_u)
    INTEGER, INTENT(IN   ) :: k_order
    REAL(KIND = wp), INTENT(IN   ) :: pdt
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: puc
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: pt
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(  OUT) :: pt_u
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zcu, zdx2, zdx4
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztu1, ztu2, ztu3, ztu4
    !$ACC KERNELS
    DO jj = 2, jpjm1
      DO ji = 1, jpim1
        ztu1(ji, jj) = (pt(ji + 1, jj) - pt(ji, jj)) * r1_e1u(ji, jj) * umask(ji, jj, 1)
      END DO
      DO ji = 2, jpim1
        ztu2(ji, jj) = (ztu1(ji, jj) - ztu1(ji - 1, jj)) * r1_e1t(ji, jj)
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk(ztu2, 'T', 1.)
    !$ACC KERNELS
    DO jj = 2, jpjm1
      DO ji = 1, jpim1
        ztu3(ji, jj) = (ztu2(ji + 1, jj) - ztu2(ji, jj)) * r1_e1u(ji, jj) * umask(ji, jj, 1)
      END DO
      DO ji = 2, jpim1
        ztu4(ji, jj) = (ztu3(ji, jj) - ztu3(ji - 1, jj)) * r1_e1t(ji, jj)
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk(ztu4, 'T', 1.)
    SELECT CASE (k_order)
    CASE (1)
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = 1, jpim1
          pt_u(ji, jj) = 0.5_wp * umask(ji, jj, 1) * (pt(ji + 1, jj) + pt(ji, jj) - SIGN(1._wp, puc(ji, jj)) * (pt(ji + 1, jj) - pt(ji, jj)))
        END DO
      END DO
      !$ACC END KERNELS
    CASE (2)
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = 1, jpim1
          zcu = puc(ji, jj) * r1_e2u(ji, jj) * pdt * r1_e1u(ji, jj)
          pt_u(ji, jj) = 0.5_wp * umask(ji, jj, 1) * (pt(ji + 1, jj) + pt(ji, jj) - zcu * (pt(ji + 1, jj) - pt(ji, jj)))
        END DO
      END DO
      !$ACC END KERNELS
    CASE (3)
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = 1, jpim1
          zcu = puc(ji, jj) * r1_e2u(ji, jj) * pdt * r1_e1u(ji, jj)
          zdx2 = e1u(ji, jj) * e1u(ji, jj)
          pt_u(ji, jj) = 0.5_wp * umask(ji, jj, 1) * ((pt(ji + 1, jj) + pt(ji, jj) - zcu * (pt(ji + 1, jj) - pt(ji, jj))) + z1_6 * zdx2 * (zcu * zcu - 1._wp) * (ztu2(ji + 1, jj) + ztu2(ji, jj) - SIGN(1._wp, zcu) * (ztu2(ji + 1, jj) - ztu2(ji, jj))))
        END DO
      END DO
      !$ACC END KERNELS
    CASE (4)
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = 1, jpim1
          zcu = puc(ji, jj) * r1_e2u(ji, jj) * pdt * r1_e1u(ji, jj)
          zdx2 = e1u(ji, jj) * e1u(ji, jj)
          pt_u(ji, jj) = 0.5_wp * umask(ji, jj, 1) * ((pt(ji + 1, jj) + pt(ji, jj) - zcu * (pt(ji + 1, jj) - pt(ji, jj))) + z1_6 * zdx2 * (zcu * zcu - 1._wp) * (ztu2(ji + 1, jj) + ztu2(ji, jj) - 0.5_wp * zcu * (ztu2(ji + 1, jj) - ztu2(ji, jj))))
        END DO
      END DO
      !$ACC END KERNELS
    CASE (5)
      !$ACC KERNELS
      DO jj = 2, jpjm1
        DO ji = 1, jpim1
          zcu = puc(ji, jj) * r1_e2u(ji, jj) * pdt * r1_e1u(ji, jj)
          zdx2 = e1u(ji, jj) * e1u(ji, jj)
          zdx4 = zdx2 * zdx2
          pt_u(ji, jj) = 0.5_wp * umask(ji, jj, 1) * ((pt(ji + 1, jj) + pt(ji, jj) - zcu * (pt(ji + 1, jj) - pt(ji, jj))) + z1_6 * zdx2 * (zcu * zcu - 1._wp) * (ztu2(ji + 1, jj) + ztu2(ji, jj) - 0.5_wp * zcu * (ztu2(ji + 1, jj) - ztu2(ji, jj))) + z1_120 * zdx4 * (zcu * zcu - 1._wp) * (zcu * zcu - 4._wp) * (ztu4(ji + 1, jj) + ztu4(ji, jj) - SIGN(1._wp, zcu) * (ztu4(ji + 1, jj) - ztu4(ji, jj))))
        END DO
      END DO
      !$ACC END KERNELS
    END SELECT
  END SUBROUTINE ultimate_x
  SUBROUTINE ultimate_y(k_order, pdt, pt, pvc, pt_v)
    INTEGER, INTENT(IN   ) :: k_order
    REAL(KIND = wp), INTENT(IN   ) :: pdt
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: pvc
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: pt
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(  OUT) :: pt_v
    INTEGER :: ji, jj
    REAL(KIND = wp) :: zcv, zdy2, zdy4
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: ztv1, ztv2, ztv3, ztv4
    !$ACC KERNELS
    DO jj = 1, jpjm1
      DO ji = 2, jpim1
        ztv1(ji, jj) = (pt(ji, jj + 1) - pt(ji, jj)) * r1_e2v(ji, jj) * vmask(ji, jj, 1)
      END DO
    END DO
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        ztv2(ji, jj) = (ztv1(ji, jj) - ztv1(ji, jj - 1)) * r1_e2t(ji, jj)
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk(ztv2, 'T', 1.)
    !$ACC KERNELS
    DO jj = 1, jpjm1
      DO ji = 2, jpim1
        ztv3(ji, jj) = (ztv2(ji, jj + 1) - ztv2(ji, jj)) * r1_e2v(ji, jj) * vmask(ji, jj, 1)
      END DO
    END DO
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        ztv4(ji, jj) = (ztv3(ji, jj) - ztv3(ji, jj - 1)) * r1_e2t(ji, jj)
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk(ztv4, 'T', 1.)
    SELECT CASE (k_order)
    CASE (1)
      !$ACC KERNELS
      DO jj = 1, jpjm1
        DO ji = 2, jpim1
          pt_v(ji, jj) = 0.5_wp * vmask(ji, jj, 1) * ((pt(ji, jj + 1) + pt(ji, jj)) - SIGN(1._wp, pvc(ji, jj)) * (pt(ji, jj + 1) - pt(ji, jj)))
        END DO
      END DO
      !$ACC END KERNELS
    CASE (2)
      !$ACC KERNELS
      DO jj = 1, jpjm1
        DO ji = 2, jpim1
          zcv = pvc(ji, jj) * r1_e1v(ji, jj) * pdt * r1_e2v(ji, jj)
          pt_v(ji, jj) = 0.5_wp * vmask(ji, jj, 1) * ((pt(ji, jj + 1) + pt(ji, jj)) - zcv * (pt(ji, jj + 1) - pt(ji, jj)))
        END DO
      END DO
      !$ACC END KERNELS
      CALL lbc_lnk(pt_v, 'V', 1.)
    CASE (3)
      !$ACC KERNELS
      DO jj = 1, jpjm1
        DO ji = 2, jpim1
          zcv = pvc(ji, jj) * r1_e1v(ji, jj) * pdt * r1_e2v(ji, jj)
          zdy2 = e2v(ji, jj) * e2v(ji, jj)
          pt_v(ji, jj) = 0.5_wp * vmask(ji, jj, 1) * ((pt(ji, jj + 1) + pt(ji, jj) - zcv * (pt(ji, jj + 1) - pt(ji, jj))) + z1_6 * zdy2 * (zcv * zcv - 1._wp) * (ztv2(ji, jj + 1) + ztv2(ji, jj) - SIGN(1._wp, zcv) * (ztv2(ji, jj + 1) - ztv2(ji, jj))))
        END DO
      END DO
      !$ACC END KERNELS
    CASE (4)
      !$ACC KERNELS
      DO jj = 1, jpjm1
        DO ji = 2, jpim1
          zcv = pvc(ji, jj) * r1_e1v(ji, jj) * pdt * r1_e2v(ji, jj)
          zdy2 = e2v(ji, jj) * e2v(ji, jj)
          pt_v(ji, jj) = 0.5_wp * vmask(ji, jj, 1) * ((pt(ji, jj + 1) + pt(ji, jj) - zcv * (pt(ji, jj + 1) - pt(ji, jj))) + z1_6 * zdy2 * (zcv * zcv - 1._wp) * (ztv2(ji, jj + 1) + ztv2(ji, jj) - 0.5_wp * zcv * (ztv2(ji, jj + 1) - ztv2(ji, jj))))
        END DO
      END DO
      !$ACC END KERNELS
    CASE (5)
      !$ACC KERNELS
      DO jj = 1, jpjm1
        DO ji = 2, jpim1
          zcv = pvc(ji, jj) * r1_e1v(ji, jj) * pdt * r1_e2v(ji, jj)
          zdy2 = e2v(ji, jj) * e2v(ji, jj)
          zdy4 = zdy2 * zdy2
          pt_v(ji, jj) = 0.5_wp * vmask(ji, jj, 1) * ((pt(ji, jj + 1) + pt(ji, jj) - zcv * (pt(ji, jj + 1) - pt(ji, jj))) + z1_6 * zdy2 * (zcv * zcv - 1._wp) * (ztv2(ji, jj + 1) + ztv2(ji, jj) - 0.5_wp * zcv * (ztv2(ji, jj + 1) - ztv2(ji, jj))) + z1_120 * zdy4 * (zcv * zcv - 1._wp) * (zcv * zcv - 4._wp) * (ztv4(ji, jj + 1) + ztv4(ji, jj) - SIGN(1._wp, zcv) * (ztv4(ji, jj + 1) - ztv4(ji, jj))))
        END DO
      END DO
      !$ACC END KERNELS
    END SELECT
  END SUBROUTINE ultimate_y
  SUBROUTINE nonosc_2d(pbef, paa, pbb, paft, pdt)
    REAL(KIND = wp), INTENT(IN   ) :: pdt
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: pbef, paft
    REAL(KIND = wp), DIMENSION(jpi, jpj), INTENT(INOUT) :: paa, pbb
    INTEGER :: ji, jj
    INTEGER :: ikm1
    REAL(KIND = wp) :: zpos, zneg, zbt, za, zb, zc, zbig, zsml, z1_dt
    REAL(KIND = wp) :: zau, zbu, zcu, zav, zbv, zcv, zup, zdo
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zbetup, zbetdo, zbup, zbdo, zmsk, zdiv
    !$ACC KERNELS
    zbig = 1.E+40_wp
    zsml = 1.E-15_wp
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zdiv(ji, jj) = - (paa(ji, jj) - paa(ji - 1, jj) + pbb(ji, jj) - pbb(ji, jj - 1))
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk(zdiv, 'T', 1.)
    WHERE (pbef(:, :) == 0._wp .AND. paft(:, :) == 0._wp .AND. zdiv(:, :) == 0._wp)
      zmsk(:, :) = 0._wp
    ELSEWHERE
      zmsk(:, :) = 1._wp * tmask(:, :, 1)
    END WHERE
    !$ACC KERNELS
    zbup(:, :) = MAX(pbef(:, :) * zmsk(:, :) - zbig * (1.E0 - zmsk(:, :)), paft(:, :) * zmsk(:, :) - zbig * (1.E0 - zmsk(:, :)))
    zbdo(:, :) = MIN(pbef(:, :) * zmsk(:, :) + zbig * (1.E0 - zmsk(:, :)), paft(:, :) * zmsk(:, :) + zbig * (1.E0 - zmsk(:, :)))
    z1_dt = 1._wp / pdt
    DO jj = 2, jpjm1
      DO ji = 2, jpim1
        zup = MAX(zbup(ji, jj), zbup(ji - 1, jj), zbup(ji + 1, jj), zbup(ji, jj - 1), zbup(ji, jj + 1))
        zdo = MIN(zbdo(ji, jj), zbdo(ji - 1, jj), zbdo(ji + 1, jj), zbdo(ji, jj - 1), zbdo(ji, jj + 1))
        zpos = MAX(0., paa(ji - 1, jj)) - MIN(0., paa(ji, jj)) + MAX(0., pbb(ji, jj - 1)) - MIN(0., pbb(ji, jj))
        zneg = MAX(0., paa(ji, jj)) - MIN(0., paa(ji - 1, jj)) + MAX(0., pbb(ji, jj)) - MIN(0., pbb(ji, jj - 1))
        zbt = e1e2t(ji, jj) * z1_dt
        zbetup(ji, jj) = (zup - paft(ji, jj)) / (zpos + zsml) * zbt
        zbetdo(ji, jj) = (paft(ji, jj) - zdo) / (zneg + zsml) * zbt
      END DO
    END DO
    !$ACC END KERNELS
    CALL lbc_lnk_multi(zbetup, 'T', 1., zbetdo, 'T', 1.)
    !$ACC KERNELS
    DO jj = 2, jpjm1
      DO ji = 1, jpim1
        zau = MIN(1._wp, zbetdo(ji, jj), zbetup(ji + 1, jj))
        zbu = MIN(1._wp, zbetup(ji, jj), zbetdo(ji + 1, jj))
        zcu = 0.5 + SIGN(0.5, paa(ji, jj))
        paa(ji, jj) = paa(ji, jj) * (zcu * zau + (1._wp - zcu) * zbu)
      END DO
    END DO
    DO jj = 1, jpjm1
      DO ji = 2, jpim1
        zav = MIN(1._wp, zbetdo(ji, jj), zbetup(ji, jj + 1))
        zbv = MIN(1._wp, zbetup(ji, jj), zbetdo(ji, jj + 1))
        zcv = 0.5 + SIGN(0.5, pbb(ji, jj))
        pbb(ji, jj) = pbb(ji, jj) * (zcv * zav + (1._wp - zcv) * zbv)
      END DO
    END DO
    !$ACC END KERNELS
  END SUBROUTINE nonosc_2d
END MODULE icedyn_adv_umx