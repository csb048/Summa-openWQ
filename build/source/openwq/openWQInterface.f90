interface
    function create_openwq_c() bind(C, name="create_openwq")
        use iso_c_binding
        implicit none
        type(c_ptr) :: create_openwq_c
    end function

    function openwq_decl_c(openWQ, num_hru, num_layers_canopy, num_layers_matricHead, &
        num_layers_aquifer, num_layers_volFracWat) bind(C, name="openwq_decl")
        use iso_c_binding
        implicit none
        integer(c_int) :: openwq_decl_c ! returns a return value of 0 (success) or -1 (failure)
        type(c_ptr), intent(in), value :: openWQ
        integer(c_int), intent(in), value  :: num_hru
        integer(c_int), intent(in), value  :: num_layers_canopy
        integer(c_int), intent(in), value  :: num_layers_matricHead
        integer(c_int), intent(in), value  :: num_layers_aquifer
        integer(c_int), intent(in), value  :: num_layers_volFracWat
    end function

    function openwq_run_time_start_c(openWQ, numHRU, simtime, &
        soilMoisture, soilTemp, airTemp, swe_vol, canopyWat_vol, &
        matricHead_vol, aquiferStorage_vol) bind(C, name="openwq_run_time_start")
        use iso_c_binding
        implicit none
        integer(c_int)                       :: openwq_run_time_start_c ! returns 0 (success) or -1 (failure)
        type(c_ptr),    intent(in), value    :: openWQ
        integer(c_int), intent(in), value    :: numHRU
        integer(c_int), intent(in)           :: simtime(5)
        real(c_double), intent(in)           :: soilMoisture(numHRU)
        real(c_double), intent(in)           :: soilTemp(numHRU)
        real(c_double), intent(in)           :: airTemp(numHRU)
        real(c_double), intent(in)           :: swe_vol(numHRU)
        real(c_double), intent(in)           :: canopyWat_vol(numHRU)
        real(c_double), intent(in)           :: matricHead_vol(numHRU)
        real(c_double), intent(in)           :: aquiferStorage_vol(numHRU)
    end function

    function openwq_run_space_c(openWQ,simtime,source,ix_s,iy_s,iz_s,recipient,ix_r,iy_r,iz_r,wflux_s2r,wmass_source) bind(C, name="openwq_run_space")
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

    function openwq_run_time_end_c(openWQ,simtime) bind(C, name="openwq_run_time_end")
        USE iso_c_binding
        implicit none
        integer(c_int) :: openwq_run_time_end_c ! returns 0 (success) or -1 (failure)
        type(c_ptr),    intent(in), value   :: openWQ
        integer(c_int), intent(in)          :: simtime(5)
    end function
end interface