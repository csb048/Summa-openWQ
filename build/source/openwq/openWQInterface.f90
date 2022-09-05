interface
    function create_openwq_c() bind(C, name="create_openwq")

        use iso_c_binding
        implicit none
        type(c_ptr) :: create_openwq_c

    end function

    function openwq_decl_c( &
        openWQ, &
        num_hru, &
        nCanopy_2openwq, &
        num_layers_matricHead, &
        nAquifer_2openwq, &
        num_layers_volFracWat, &
        y_direction) bind(C, name="openwq_decl")

        use iso_c_binding
        implicit none
        integer(c_int) :: openwq_decl_c ! returns a return value of 0 (success) or -1 (failure)
        type(c_ptr), intent(in), value :: openWQ
        integer(c_int), intent(in), value  :: num_hru
        integer(c_int), intent(in), value  :: nCanopy_2openwq
        integer(c_int), intent(in), value  :: num_layers_matricHead
        integer(c_int), intent(in), value  :: nAquifer_2openwq
        integer(c_int), intent(in), value  :: num_layers_volFracWat
        integer(c_int), intent(in), value  :: y_direction

    end function

    function openwq_run_time_start_c(&
        openWQ, &
        numHRU, &
        nSnow_2openwq, &
        nSoil_2openwq, &
        simtime, &
        soilMoist_depVar, &
        soilTemp_depVar, &
        airTemp_depVar, &
        sweWatVol_stateVar, &
        canopyWatVol_stateVar, &
        soilWatVol_stateVar, &
        aquiferWatVol_stateVar) bind(C, name="openwq_run_time_start")

        use iso_c_binding
        implicit none
        integer(c_int)                       :: openwq_run_time_start_c ! returns 0 (success) or -1 (failure)
        type(c_ptr),    intent(in), value    :: openWQ
        integer(c_int), intent(in), value    :: numHRU
        integer(c_int), intent(in), value    :: nSnow_2openwq
        integer(c_int), intent(in), value    :: nSoil_2openwq
        integer(c_int), intent(in)           :: simtime(5)
        real(c_double), intent(in)           :: soilMoist_depVar(numHRU, nSoil_2openwq)
        real(c_double), intent(in)           :: soilTemp_depVar(numHRU, nSoil_2openwq)
        real(c_double), intent(in)           :: airTemp_depVar(numHRU)
        real(c_double), intent(in)           :: sweWatVol_stateVar(numHRU, nSnow_2openwq)
        real(c_double), intent(in)           :: canopyWatVol_stateVar(numHRU)
        real(c_double), intent(in)           :: soilWatVol_stateVar(numHRU, nSoil_2openwq)
        real(c_double), intent(in)           :: aquiferWatVol_stateVar(numHRU)

    end function

    function openwq_run_space_c(&
        openWQ, &
        simtime, &
        source,ix_s,iy_s,iz_s, &
        recipient,ix_r,iy_r,iz_r, &
        wflux_s2r, &
        wmass_source) bind(C, name="openwq_run_space")

        use iso_c_binding
        implicit none
        integer(c_int) :: openwq_run_space_c ! returns 0 (success) or -1 (failure)
        type(c_ptr),    intent(in), value      :: openWQ
        integer(c_int), intent(in)             :: simtime(5)
        integer(c_int), intent(in), value      :: source
        integer(c_int), intent(in), value      :: ix_s
        integer(c_int), intent(in), value      :: iy_s 
        integer(c_int), intent(in), value      :: iz_s
        integer(c_int), intent(in), value      :: recipient
        integer(c_int), intent(in), value      :: ix_r
        integer(c_int), intent(in), value      :: iy_r
        integer(c_int), intent(in), value      :: iz_r
        real(c_double), intent(in), value      :: wflux_s2r
        real(c_double), intent(in), value      :: wmass_source

    end function

    function openwq_run_space_in_c( &
        openWQ, &
        simtime, &
        recipient,ix_r,iy_r,iz_r, &
        wflux_s2r) bind(C, name="openwq_run_space_in")

        USE iso_c_binding
        implicit none
        integer(c_int) :: openwq_run_space_in_c
        type(c_ptr), intent(in), value         :: openWQ
        integer(c_int), intent(in)             :: simtime(5)
        integer(c_int), intent(in), value      :: recipient
        integer(c_int), intent(in), value      :: ix_r
        integer(c_int), intent(in), value      :: iy_r
        integer(c_int), intent(in), value      :: iz_r
        real(c_double), intent(in), value      :: wflux_s2r

    end function

    function openwq_run_time_end_c( &
        openWQ, &
        simtime) bind(C, name="openwq_run_time_end")

        USE iso_c_binding
        implicit none
        integer(c_int) :: openwq_run_time_end_c ! returns 0 (success) or -1 (failure)
        type(c_ptr),    intent(in), value   :: openWQ
        integer(c_int), intent(in)          :: simtime(5)

    end function

end interface