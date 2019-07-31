MODULE diawri
  USE oce
  USE dom_oce
  USE phycst
  USE dianam
  USE diahth
  USE dynadv, ONLY: ln_dynadv_vec
  USE icb_oce
  USE icbdia
  USE ldftra
  USE ldfdyn
  USE sbc_oce
  USE sbc_ice
  USE sbcssr
  USE sbcwave
  USE wet_dry
  USE zdf_oce
  USE zdfdrg
  USE zdfmxl
  USE lbclnk
  USE in_out_manager
  USE diatmb
  USE dia25h
  USE iom
  USE ioipsl
  USE icewri
  USE lib_mpp
  USE timing
  USE diurnal_bulk
  USE cool_skin
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: dia_wri
  PUBLIC :: dia_wri_state
  PUBLIC :: dia_wri_alloc
  INTEGER :: nid_T, nz_T, nh_T, ndim_T, ndim_hT
  INTEGER :: nb_T, ndim_bT
  INTEGER :: nid_U, nz_U, nh_U, ndim_U, ndim_hU
  INTEGER :: nid_V, nz_V, nh_V, ndim_V, ndim_hV
  INTEGER :: nid_W, nz_W, nh_W
  INTEGER :: ndex(1)
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_hT, ndex_hU, ndex_hV
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_T, ndex_U, ndex_V
  INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_bT
  CONTAINS
  INTEGER FUNCTION dia_wri_alloc()
    INTEGER, DIMENSION(2) :: ierr
    ierr = 0
    ALLOCATE(ndex_hT(jpi * jpj), ndex_T(jpi * jpj * jpk), ndex_hU(jpi * jpj), ndex_U(jpi * jpj * jpk), ndex_hV(jpi * jpj), ndex_V(jpi * jpj * jpk), STAT = ierr(1))
    dia_wri_alloc = MAXVAL(ierr)
    IF (lk_mpp) CALL mpp_sum(dia_wri_alloc)
  END FUNCTION dia_wri_alloc
  SUBROUTINE dia_wri(kt)
    INTEGER, INTENT( IN ) :: kt
    LOGICAL :: ll_print = .FALSE.
    CHARACTER(LEN = 40) :: clhstnam, clop, clmx
    INTEGER :: inum = 11
    INTEGER :: ji, jj, jk
    INTEGER :: ierr
    INTEGER :: iimi, iima, ipk, it, itmod, ijmi, ijma
    INTEGER :: jn, ierror
    REAL(KIND = wp) :: zsto, zout, zmax, zjulian
    REAL(KIND = wp), DIMENSION(jpi, jpj) :: zw2d
    REAL(KIND = wp), DIMENSION(jpi, jpj, jpk) :: zw3d
    RETURN
    IF (ln_timing) CALL timing_start('dia_wri')
    IF (ninist == 1) THEN
      CALL dia_wri_state('output.init', kt)
      ninist = 0
    END IF
    ll_print = .FALSE.
    ll_print = ll_print .AND. lwp
    clop = "x"
    zsto = rdt
    clop = "ave(" // TRIM(clop) // ")"
    zout = nwrite * rdt
    zmax = (nitend - nit000 + 1) * rdt
    iimi = 1
    iima = jpi
    ijmi = 1
    ijma = jpj
    ipk = jpk
    it = kt
    itmod = kt - nit000 + 1
    IF (kt == nit000) THEN
      CALL ymds2ju(nyear, nmonth, nday, rdt, zjulian)
      zjulian = zjulian - adatrj
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'Date 0 used :', nit000, ' YEAR ', nyear, ' MONTH ', nmonth, ' DAY ', nday, 'Julian day : ', zjulian
      IF (lwp) WRITE(numout, FMT = *) ' indexes of zoom = ', iimi, iima, ijmi, ijma, ' limit storage in depth = ', ipk
      IF (lwp) THEN
        CALL dia_nam(clhstnam, nwrite, ' ')
        CALL ctl_opn(inum, 'date.file', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', - 1, numout, lwp, narea)
        WRITE(inum, FMT = *) clhstnam
        CLOSE(UNIT = inum)
      END IF
      CALL dia_nam(clhstnam, nwrite, 'grid_T')
      IF (lwp) WRITE(numout, FMT = *) " Name of NETCDF file ", clhstnam
      CALL histbeg(clhstnam, jpi, glamt, jpj, gphit, iimi, iima - iimi + 1, ijmi, ijma - ijmi + 1, nit000 - 1, zjulian, rdt, nh_T, nid_T, domain_id = nidom, snc4chunks = snc4set)
      CALL histvert(nid_T, "deptht", "Vertical T levels", "m", ipk, gdept_1d, nz_T, "down")
      CALL wheneq(jpi * jpj * ipk, tmask, 1, 1., ndex_T, ndim_T)
      CALL wheneq(jpi * jpj, tmask, 1, 1., ndex_hT, ndim_hT)
      IF (ln_icebergs) THEN
        ALLOCATE(ndex_bT(jpi * jpj * nclasses), STAT = ierror)
        IF (lk_mpp) CALL mpp_sum(ierror)
        IF (ierror /= 0) THEN
          CALL ctl_stop('dia_wri: failed to allocate iceberg diagnostic array')
          RETURN
        END IF
        CALL histvert(nid_T, "class", "Iceberg class", "number", nclasses, class_num, nb_T)
        ndim_bT = 3
        DO jn = 1, nclasses
          ndex_bT((jn - 1) * jpi * jpj + 1 : jn * jpi * jpj) = ndex_hT(1 : jpi * jpj)
        END DO
      END IF
      CALL dia_nam(clhstnam, nwrite, 'grid_U')
      IF (lwp) WRITE(numout, FMT = *) " Name of NETCDF file ", clhstnam
      CALL histbeg(clhstnam, jpi, glamu, jpj, gphiu, iimi, iima - iimi + 1, ijmi, ijma - ijmi + 1, nit000 - 1, zjulian, rdt, nh_U, nid_U, domain_id = nidom, snc4chunks = snc4set)
      CALL histvert(nid_U, "depthu", "Vertical U levels", "m", ipk, gdept_1d, nz_U, "down")
      CALL wheneq(jpi * jpj * ipk, umask, 1, 1., ndex_U, ndim_U)
      CALL wheneq(jpi * jpj, umask, 1, 1., ndex_hU, ndim_hU)
      CALL dia_nam(clhstnam, nwrite, 'grid_V')
      IF (lwp) WRITE(numout, FMT = *) " Name of NETCDF file ", clhstnam
      CALL histbeg(clhstnam, jpi, glamv, jpj, gphiv, iimi, iima - iimi + 1, ijmi, ijma - ijmi + 1, nit000 - 1, zjulian, rdt, nh_V, nid_V, domain_id = nidom, snc4chunks = snc4set)
      CALL histvert(nid_V, "depthv", "Vertical V levels", "m", ipk, gdept_1d, nz_V, "down")
      CALL wheneq(jpi * jpj * ipk, vmask, 1, 1., ndex_V, ndim_V)
      CALL wheneq(jpi * jpj, vmask, 1, 1., ndex_hV, ndim_hV)
      CALL dia_nam(clhstnam, nwrite, 'grid_W')
      IF (lwp) WRITE(numout, FMT = *) " Name of NETCDF file ", clhstnam
      CALL histbeg(clhstnam, jpi, glamt, jpj, gphit, iimi, iima - iimi + 1, ijmi, ijma - ijmi + 1, nit000 - 1, zjulian, rdt, nh_W, nid_W, domain_id = nidom, snc4chunks = snc4set)
      CALL histvert(nid_W, "depthw", "Vertical W levels", "m", ipk, gdepw_1d, nz_W, "down")
      CALL histdef(nid_T, "votemper", "Temperature", "C", jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout)
      CALL histdef(nid_T, "vosaline", "Salinity", "PSU", jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout)
      IF (.NOT. ln_linssh) THEN
        CALL histdef(nid_T, "vovvle3t", "Level thickness", "m", jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout)
        CALL histdef(nid_T, "vovvldep", "T point depth", "m", jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout)
        CALL histdef(nid_T, "vovvldef", "Squared level deformation", "%^2", jpi, jpj, nh_T, ipk, 1, ipk, nz_T, 32, clop, zsto, zout)
      END IF
      CALL histdef(nid_T, "sosstsst", "Sea Surface temperature", "C", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "sosaline", "Sea Surface Salinity", "PSU", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "sossheig", "Sea Surface Height", "m", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "sowaflup", "Net Upward Water Flux", "Kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "sorunoff", "River runoffs", "Kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "sosfldow", "downward salt flux", "PSU/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      IF (ln_linssh) THEN
        CALL histdef(nid_T, "sosst_cd", "Concentration/Dilution term on temperature", "KgC/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "sosss_cd", "Concentration/Dilution term on salinity", "KgPSU/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      END IF
      CALL histdef(nid_T, "sohefldo", "Net Downward Heat Flux", "W/m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "soshfldo", "Shortwave Radiation", "W/m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "somixhgt", "Turbocline Depth", "m", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "somxl010", "Mixed Layer Depth 0.01", "m", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "soicecov", "Ice fraction", "[0,1]", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histdef(nid_T, "sowindsp", "wind speed at 10m", "m/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      IF (ln_icebergs) THEN
        CALL histdef(nid_T, "calving", "calving mass input", "kg/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "calving_heat", "calving heat flux", "XXXX", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "berg_floating_melt", "Melt rate of icebergs + bits", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "berg_stored_ice", "Accumulated ice mass by class", "kg", jpi, jpj, nh_T, nclasses, 1, nclasses, nb_T, 32, clop, zsto, zout)
        IF (ln_bergdia) THEN
          CALL histdef(nid_T, "berg_melt", "Melt rate of icebergs", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "berg_buoy_melt", "Buoyancy component of iceberg melt rate", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "berg_eros_melt", "Erosion component of iceberg melt rate", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "berg_conv_melt", "Convective component of iceberg melt rate", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "berg_virtual_area", "Virtual coverage by icebergs", "m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "bits_src", "Mass source of bergy bits", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "bits_melt", "Melt rate of bergy bits", "kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "bits_mass", "Bergy bit density field", "kg/m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "berg_mass", "Iceberg density field", "kg/m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
          CALL histdef(nid_T, "berg_real_calving", "Calving into iceberg class", "kg/s", jpi, jpj, nh_T, nclasses, 1, nclasses, nb_T, 32, clop, zsto, zout)
        END IF
      END IF
      IF (.NOT. ln_cpl) THEN
        CALL histdef(nid_T, "sohefldp", "Surface Heat Flux: Damping", "W/m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "sowafldp", "Surface Water Flux: Damping", "Kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "sosafldp", "Surface salt flux: damping", "Kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      END IF
      IF (ln_cpl .AND. nn_ice <= 1) THEN
        CALL histdef(nid_T, "sohefldp", "Surface Heat Flux: Damping", "W/m2", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "sowafldp", "Surface Water Flux: Damping", "Kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
        CALL histdef(nid_T, "sosafldp", "Surface salt flux: Damping", "Kg/m2/s", jpi, jpj, nh_T, 1, 1, 1, - 99, 32, clop, zsto, zout)
      END IF
      clmx = "l_max(only(x))"
      CALL histend(nid_T, snc4chunks = snc4set)
      CALL histdef(nid_U, "vozocrtx", "Zonal Current", "m/s", jpi, jpj, nh_U, ipk, 1, ipk, nz_U, 32, clop, zsto, zout)
      IF (ln_wave .AND. ln_sdw) THEN
        CALL histdef(nid_U, "sdzocrtx", "Stokes Drift Zonal Current", "m/s", jpi, jpj, nh_U, ipk, 1, ipk, nz_U, 32, clop, zsto, zout)
      END IF
      CALL histdef(nid_U, "sozotaux", "Wind Stress along i-axis", "N/m2", jpi, jpj, nh_U, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histend(nid_U, snc4chunks = snc4set)
      CALL histdef(nid_V, "vomecrty", "Meridional Current", "m/s", jpi, jpj, nh_V, ipk, 1, ipk, nz_V, 32, clop, zsto, zout)
      IF (ln_wave .AND. ln_sdw) THEN
        CALL histdef(nid_V, "sdmecrty", "Stokes Drift Meridional Current", "m/s", jpi, jpj, nh_V, ipk, 1, ipk, nz_V, 32, clop, zsto, zout)
      END IF
      CALL histdef(nid_V, "sometauy", "Wind Stress along j-axis", "N/m2", jpi, jpj, nh_V, 1, 1, 1, - 99, 32, clop, zsto, zout)
      CALL histend(nid_V, snc4chunks = snc4set)
      CALL histdef(nid_W, "vovecrtz", "Vertical Velocity", "m/s", jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout)
      CALL histdef(nid_W, "votkeavt", "Vertical Eddy Diffusivity", "m2/s", jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout)
      CALL histdef(nid_W, "votkeavm", "Vertical Eddy Viscosity", "m2/s", jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout)
      IF (ln_zdfddm) THEN
        CALL histdef(nid_W, "voddmavs", "Salt Vertical Eddy Diffusivity", "m2/s", jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout)
      END IF
      IF (ln_wave .AND. ln_sdw) THEN
        CALL histdef(nid_W, "sdvecrtz", "Stokes Drift Vertical Current", "m/s", jpi, jpj, nh_W, ipk, 1, ipk, nz_W, 32, clop, zsto, zout)
      END IF
      CALL histend(nid_W, snc4chunks = snc4set)
      IF (lwp) WRITE(numout, FMT = *)
      IF (lwp) WRITE(numout, FMT = *) 'End of NetCDF Initialization'
      IF (ll_print) CALL FLUSH(numout)
    END IF
    IF (lwp .AND. MOD(itmod, nwrite) == 0) THEN
      WRITE(numout, FMT = *) 'dia_wri : write model outputs in NetCDF files at ', kt, 'time-step'
      WRITE(numout, FMT = *) '~~~~~~ '
    END IF
    IF (.NOT. ln_linssh) THEN
      CALL histwrite(nid_T, "votemper", it, tsn(:, :, :, jp_tem) * e3t_n(:, :, :), ndim_T, ndex_T)
      CALL histwrite(nid_T, "vosaline", it, tsn(:, :, :, jp_sal) * e3t_n(:, :, :), ndim_T, ndex_T)
      CALL histwrite(nid_T, "sosstsst", it, tsn(:, :, 1, jp_tem) * e3t_n(:, :, 1), ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "sosaline", it, tsn(:, :, 1, jp_sal) * e3t_n(:, :, 1), ndim_hT, ndex_hT)
    ELSE
      CALL histwrite(nid_T, "votemper", it, tsn(:, :, :, jp_tem), ndim_T, ndex_T)
      CALL histwrite(nid_T, "vosaline", it, tsn(:, :, :, jp_sal), ndim_T, ndex_T)
      CALL histwrite(nid_T, "sosstsst", it, tsn(:, :, 1, jp_tem), ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "sosaline", it, tsn(:, :, 1, jp_sal), ndim_hT, ndex_hT)
    END IF
    IF (.NOT. ln_linssh) THEN
      !$ACC KERNELS
      zw3d(:, :, :) = ((e3t_n(:, :, :) - e3t_0(:, :, :)) / e3t_0(:, :, :) * 100 * tmask(:, :, :)) ** 2
      !$ACC END KERNELS
      CALL histwrite(nid_T, "vovvle3t", it, e3t_n(:, :, :), ndim_T, ndex_T)
      CALL histwrite(nid_T, "vovvldep", it, gdept_n(:, :, :), ndim_T, ndex_T)
      CALL histwrite(nid_T, "vovvldef", it, zw3d, ndim_T, ndex_T)
    END IF
    CALL histwrite(nid_T, "sossheig", it, sshn, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "sowaflup", it, (emp - rnf), ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "sorunoff", it, rnf, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "sosfldow", it, sfx, ndim_hT, ndex_hT)
    IF (ln_linssh) THEN
      !$ACC KERNELS
      zw2d(:, :) = emp(:, :) * tsn(:, :, 1, jp_tem)
      !$ACC END KERNELS
      CALL histwrite(nid_T, "sosst_cd", it, zw2d, ndim_hT, ndex_hT)
      !$ACC KERNELS
      zw2d(:, :) = emp(:, :) * tsn(:, :, 1, jp_sal)
      !$ACC END KERNELS
      CALL histwrite(nid_T, "sosss_cd", it, zw2d, ndim_hT, ndex_hT)
    END IF
    CALL histwrite(nid_T, "sohefldo", it, qns + qsr, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "soshfldo", it, qsr, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "somixhgt", it, hmld, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "somxl010", it, hmlp, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "soicecov", it, fr_i, ndim_hT, ndex_hT)
    CALL histwrite(nid_T, "sowindsp", it, wndm, ndim_hT, ndex_hT)
    IF (ln_icebergs) THEN
      CALL histwrite(nid_T, "calving", it, berg_grid % calving, ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "calving_heat", it, berg_grid % calving_hflx, ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "berg_floating_melt", it, berg_grid % floating_melt, ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "berg_stored_ice", it, berg_grid % stored_ice, ndim_bT, ndex_bT)
      IF (ln_bergdia) THEN
        CALL histwrite(nid_T, "berg_melt", it, berg_melt, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "berg_buoy_melt", it, buoy_melt, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "berg_eros_melt", it, eros_melt, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "berg_conv_melt", it, conv_melt, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "berg_virtual_area", it, virtual_area, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "bits_src", it, bits_src, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "bits_melt", it, bits_melt, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "bits_mass", it, bits_mass, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "berg_mass", it, berg_mass, ndim_hT, ndex_hT)
        CALL histwrite(nid_T, "berg_real_calving", it, real_calving, ndim_bT, ndex_bT)
      END IF
    END IF
    IF (.NOT. ln_cpl) THEN
      CALL histwrite(nid_T, "sohefldp", it, qrp, ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "sowafldp", it, erp, ndim_hT, ndex_hT)
      IF (ln_ssr) zw2d(:, :) = erp(:, :) * tsn(:, :, 1, jp_sal) * tmask(:, :, 1)
      CALL histwrite(nid_T, "sosafldp", it, zw2d, ndim_hT, ndex_hT)
    END IF
    IF (ln_cpl .AND. nn_ice <= 1) THEN
      CALL histwrite(nid_T, "sohefldp", it, qrp, ndim_hT, ndex_hT)
      CALL histwrite(nid_T, "sowafldp", it, erp, ndim_hT, ndex_hT)
      IF (ln_ssr) zw2d(:, :) = erp(:, :) * tsn(:, :, 1, jp_sal) * tmask(:, :, 1)
      CALL histwrite(nid_T, "sosafldp", it, zw2d, ndim_hT, ndex_hT)
    END IF
    CALL histwrite(nid_U, "vozocrtx", it, un, ndim_U, ndex_U)
    CALL histwrite(nid_U, "sozotaux", it, utau, ndim_hU, ndex_hU)
    CALL histwrite(nid_V, "vomecrty", it, vn, ndim_V, ndex_V)
    CALL histwrite(nid_V, "sometauy", it, vtau, ndim_hV, ndex_hV)
    CALL histwrite(nid_W, "vovecrtz", it, wn, ndim_T, ndex_T)
    CALL histwrite(nid_W, "votkeavt", it, avt, ndim_T, ndex_T)
    CALL histwrite(nid_W, "votkeavm", it, avm, ndim_T, ndex_T)
    IF (ln_zdfddm) THEN
      CALL histwrite(nid_W, "voddmavs", it, avs, ndim_T, ndex_T)
    END IF
    IF (ln_wave .AND. ln_sdw) THEN
      CALL histwrite(nid_U, "sdzocrtx", it, usd, ndim_U, ndex_U)
      CALL histwrite(nid_V, "sdmecrty", it, vsd, ndim_V, ndex_V)
      CALL histwrite(nid_W, "sdvecrtz", it, wsd, ndim_T, ndex_T)
    END IF
    IF (kt == nitend) THEN
      CALL histclo(nid_T)
      CALL histclo(nid_U)
      CALL histclo(nid_V)
      CALL histclo(nid_W)
    END IF
    IF (ln_timing) CALL timing_stop('dia_wri')
  END SUBROUTINE dia_wri
  SUBROUTINE dia_wri_state(cdfile_name, kt)
    CHARACTER(LEN = *), INTENT( IN ) :: cdfile_name
    INTEGER, INTENT( IN ) :: kt
    CHARACTER(LEN = 32) :: clname
    CHARACTER(LEN = 40) :: clop
    INTEGER :: id_i, nz_i, nh_i
    INTEGER, DIMENSION(1) :: idex
    REAL(KIND = wp) :: zsto, zout, zmax, zjulian
    clname = cdfile_name
    IF (.NOT. Agrif_Root()) clname = TRIM(Agrif_CFixed()) // '_' // TRIM(clname)
    zsto = rdt
    clop = "inst(x)"
    zout = rdt
    zmax = (nitend - nit000 + 1) * rdt
    IF (lwp) WRITE(numout, FMT = *)
    IF (lwp) WRITE(numout, FMT = *) 'dia_wri_state : single instantaneous ocean state'
    IF (lwp) WRITE(numout, FMT = *) '~~~~~~~~~~~~~   and forcing fields file created '
    IF (lwp) WRITE(numout, FMT = *) '                and named :', clname, '.nc'
    CALL ymds2ju(nyear, nmonth, nday, rdt, zjulian)
    zjulian = zjulian - adatrj
    CALL histbeg(clname, jpi, glamt, jpj, gphit, 1, jpi, 1, jpj, nit000 - 1, zjulian, rdt, nh_i, id_i, domain_id = nidom, snc4chunks = snc4set)
    CALL histvert(id_i, "deptht", "Vertical T levels", "m", jpk, gdept_1d, nz_i, "down")
    CALL histdef(id_i, "vosaline", "Salinity", "PSU", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
    CALL histdef(id_i, "votemper", "Temperature", "C", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
    CALL histdef(id_i, "sossheig", "Sea Surface Height", "m", jpi, jpj, nh_i, 1, 1, 1, nz_i, 32, clop, zsto, zout)
    CALL histdef(id_i, "vozocrtx", "Zonal Current", "m/s", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
    CALL histdef(id_i, "vomecrty", "Meridional Current", "m/s", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
    CALL histdef(id_i, "vovecrtz", "Vertical Velocity", "m/s", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
    IF (ALLOCATED(ahtu)) THEN
      CALL histdef(id_i, "ahtu", "u-eddy diffusivity", "m2/s", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
      CALL histdef(id_i, "ahtv", "v-eddy diffusivity", "m2/s", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
    END IF
    IF (ALLOCATED(ahmt)) THEN
      CALL histdef(id_i, "ahmt", "t-eddy viscosity", "m2/s", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
      CALL histdef(id_i, "ahmf", "f-eddy viscosity", "m2/s", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
    END IF
    CALL histdef(id_i, "sowaflup", "Net Upward Water Flux", "Kg/m2/S", jpi, jpj, nh_i, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(id_i, "sohefldo", "Net Downward Heat Flux", "W/m2", jpi, jpj, nh_i, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(id_i, "soshfldo", "Shortwave Radiation", "W/m2", jpi, jpj, nh_i, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(id_i, "soicecov", "Ice fraction", "[0,1]", jpi, jpj, nh_i, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(id_i, "sozotaux", "Zonal Wind Stress", "N/m2", jpi, jpj, nh_i, 1, 1, 1, - 99, 32, clop, zsto, zout)
    CALL histdef(id_i, "sometauy", "Meridional Wind Stress", "N/m2", jpi, jpj, nh_i, 1, 1, 1, - 99, 32, clop, zsto, zout)
    IF (.NOT. ln_linssh) THEN
      CALL histdef(id_i, "vovvldep", "T point depth", "m", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
      CALL histdef(id_i, "vovvle3t", "T point thickness", "m", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
    END IF
    IF (ln_wave .AND. ln_sdw) THEN
      CALL histdef(id_i, "sdzocrtx", "Stokes Drift Zonal", "m/s", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
      CALL histdef(id_i, "sdmecrty", "Stokes Drift Merid", "m/s", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
      CALL histdef(id_i, "sdvecrtz", "Stokes Drift Vert", "m/s", jpi, jpj, nh_i, jpk, 1, jpk, nz_i, 32, clop, zsto, zout)
    END IF
    IF (nn_ice == 2) THEN
      CALL ice_wri_state(kt, id_i, nh_i)
    END IF
    idex(1) = 1
    CALL histwrite(id_i, "votemper", kt, tsn(:, :, :, jp_tem), jpi * jpj * jpk, idex)
    CALL histwrite(id_i, "vosaline", kt, tsn(:, :, :, jp_sal), jpi * jpj * jpk, idex)
    CALL histwrite(id_i, "sossheig", kt, sshn, jpi * jpj, idex)
    CALL histwrite(id_i, "vozocrtx", kt, un, jpi * jpj * jpk, idex)
    CALL histwrite(id_i, "vomecrty", kt, vn, jpi * jpj * jpk, idex)
    CALL histwrite(id_i, "vovecrtz", kt, wn, jpi * jpj * jpk, idex)
    IF (ALLOCATED(ahtu)) THEN
      CALL histwrite(id_i, "ahtu", kt, ahtu, jpi * jpj * jpk, idex)
      CALL histwrite(id_i, "ahtv", kt, ahtv, jpi * jpj * jpk, idex)
    END IF
    IF (ALLOCATED(ahmt)) THEN
      CALL histwrite(id_i, "ahmt", kt, ahmt, jpi * jpj * jpk, idex)
      CALL histwrite(id_i, "ahmf", kt, ahmf, jpi * jpj * jpk, idex)
    END IF
    CALL histwrite(id_i, "sowaflup", kt, emp - rnf, jpi * jpj, idex)
    CALL histwrite(id_i, "sohefldo", kt, qsr + qns, jpi * jpj, idex)
    CALL histwrite(id_i, "soshfldo", kt, qsr, jpi * jpj, idex)
    CALL histwrite(id_i, "soicecov", kt, fr_i, jpi * jpj, idex)
    CALL histwrite(id_i, "sozotaux", kt, utau, jpi * jpj, idex)
    CALL histwrite(id_i, "sometauy", kt, vtau, jpi * jpj, idex)
    IF (.NOT. ln_linssh) THEN
      CALL histwrite(id_i, "vovvldep", kt, gdept_n(:, :, :), jpi * jpj * jpk, idex)
      CALL histwrite(id_i, "vovvle3t", kt, e3t_n(:, :, :), jpi * jpj * jpk, idex)
    END IF
    IF (ln_wave .AND. ln_sdw) THEN
      CALL histwrite(id_i, "sdzocrtx", kt, usd, jpi * jpj * jpk, idex)
      CALL histwrite(id_i, "sdmecrty", kt, vsd, jpi * jpj * jpk, idex)
      CALL histwrite(id_i, "sdvecrtz", kt, wsd, jpi * jpj * jpk, idex)
    END IF
    CALL histclo(id_i)
    IF (ninist /= 1) THEN
      CALL histclo(nid_T)
      CALL histclo(nid_U)
      CALL histclo(nid_V)
      CALL histclo(nid_W)
    END IF
  END SUBROUTINE dia_wri_state
END MODULE diawri