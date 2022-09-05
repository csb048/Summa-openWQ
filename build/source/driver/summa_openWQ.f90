module summa_openWQ
  USE nrtype
  USE openWQ, only:ClassWQ_OpenWQ
  USE data_types, only:gru_hru_doubleVec
  implicit none
  private
  ! Subroutines
  public :: init_openwq
  public :: run_time_start
  public :: run_time_start_go
  public :: run_space_step
  public :: run_time_end

  ! Global Data for prognostic Variables of HRUs
  type(gru_hru_doubleVec),save,public   :: progStruct_timestep_start ! copy of progStruct at the start of timestep for passing fluxes


  contains

  ! Subroutine to initalize the openWQ object
  ! putting it here to keep the SUMMA_Driver clean
subroutine init_openwq(err, message)

  USE globalData,only:openWQ_obj
  USE globalData,only:gru_struc                               ! gru-hru mapping structures
  USE globalData,only:prog_meta
  USE allocspace_module,only:allocGlobal                      ! module to allocate space for global data structures

  implicit none

  integer(i4b),intent(inout)                      :: err
  character(*),intent(inout)                      :: message         ! error messgage
  integer(i4b)                                    :: hruCount
  integer(i4b)                                    :: nCanopy_2openwq
  integer(i4b)                                    :: nSnow_2openwq
  integer(i4b)                                    :: nSoil_2openwq
  integer(i4b)                                    :: nAquifer_2openwq
  integer(i4b)                                    :: nYdirec_2openwq     ! number of layers in the y-dir (not used in summa)
  integer(i4b)                                    :: iGRU, iHRU          ! indices of GRUs and HRUs

  openwq_obj = ClassWQ_OpenWQ() ! initalize openWQ object

  ! nx -> num of HRUs)
  hruCount = sum( gru_struc(:)%hruCount )

  ! ny -> this seems to be fixes because SUMMA is based on the HRU concept, so grids are serialized)
  nYdirec_2openwq = 1

  ! Openwq nz (number of layers)
  nCanopy_2openwq = 1       ! Cannopy has only 1 layer
  nAquifer_2openwq = 1      ! GW has only 1 layer
  nSoil_2openwq = 0   ! Soil may have multiple layers, and gru-hrus may have different values
  nSnow_2openwq = 0   ! Snow has multiple layers, and gru-hrus may have different values (up to 5 layers)
  do iGRU = 1, size(gru_struc(:))
    do iHRU = 1, gru_struc(iGRU)%hruCount
      nSoil_2openwq = max( gru_struc(iGRU)%hruInfo(iHRU)%nSoil, nSoil_2openwq )
      !num_layers_volFracWat = max( gru_struc(iGRU)%hruInfo(iHRU)%nSoil, num_layers_volFracWat )
    enddo
  enddo
  nSnow_2openwq = 5 ! maximum number of snow layers

  ! intialize openWQ
  err=openwq_obj%decl(    &
    hruCount,             & ! num HRU
    nCanopy_2openwq,      & ! num layers of canopy (fixed to 1)
    nSnow_2openwq,        & ! num layers of snow (fixed to max of 5 because it varies)
    nSoil_2openwq,        & ! num layers of snoil (variable)
    nAquifer_2openwq,     & ! num layers of aquifer (fixed to 1)
    nYdirec_2openwq)             ! num of layers in y-dir (set to 1 because not used in summa)
  
  ! Create copy of state information, needed for passing to openWQ with fluxes that require
  ! the previous time_steps volume
  call allocGlobal(prog_meta, progStruct_timestep_start, err, message) 

end subroutine init_openwq
  
! Subroutine that SUMMA calls to pass varialbes that need to go to
! openWQ - the copy of progStruct is done in here
subroutine run_time_start(  &
    openWQ_obj,             & ! passing openwq object
    summa1_struc)

  USE summa_type, only: summa1_type_dec            ! master summa data type
  USE globalData, only: gru_struc

  implicit none

  ! Dummy Varialbes
  class(ClassWQ_OpenWQ), intent(in)  :: openWQ_obj
  type(summa1_type_dec), intent(in)  :: summa1_struc
  ! local variables
  integer(i4b)                       :: iGRU
  integer(i4b)                       :: iHRU
  integer(i4b)                       :: nSoil_2openwq ! maximum number of layers for soil
  integer(i4b)                       :: nSnow_2openwq ! maximum number of layers for snow)
  integer(i4b)                       :: err

  ! Get number of soil and snow layers
  ! Needs to be isolated because explicit-shaped arrays can only be defined with parameters 
  ! or int passed as an argument
  nSoil_2openwq = 0
  nSnow_2openwq = 0
  do iGRU = 1, size(gru_struc(:))
    do iHRU = 1, gru_struc(iGRU)%hruCount
      ! nSnow_2openwq = max( gru_struc(iGRU)%hruInfo(iHRU)%nSnow, nSnow_2openwq )
      nSoil_2openwq = max( gru_struc(iGRU)%hruInfo(iHRU)%nSoil, nSoil_2openwq )
    enddo
  enddo

  call run_time_start_go( &
    openwq_obj,           &
    summa1_struc,         &
    nSnow_2openwq,        &
    nSoil_2openwq)

end subroutine

subroutine run_time_start_go( &
    openWQ_obj,               &
    summa1_struc,             &
    nSnow_2openwq,            &
    nSoil_2openwq)

  USE summa_type, only: summa1_type_dec            ! master summa data type
  USE globalData, only: gru_struc
  USE var_lookup, only: iLookPROG  ! named variables for state variables
  USE var_lookup, only: iLookTIME  ! named variables for time data structure
  USE var_lookup, only: iLookATTR  ! named variables for real valued attribute data structure

  implicit none

  ! Dummy Varialbes
  class(ClassWQ_OpenWQ), intent(in)   :: openWQ_obj
  type(summa1_type_dec), intent(in)   :: summa1_struc
  ! local variables
  integer(i4b), intent(in)            :: nSnow_2openwq
  integer(i4b), intent(in)            :: nSoil_2openwq
  integer(i4b)                        :: iGRU
  integer(i4b)                        :: iHRU
  integer(i4b)                        :: ilay
  integer(i4b)                        :: iVar
  integer(i4b)                        :: iDat
  integer(i4b)                        :: openWQArrayIndex
  integer(i4b)                        :: simtime(5) ! 5 time values yy-mm-dd-hh-min
  real(rkind)                         :: airTemp_K_depVar(sum(gru_struc(:)%hruCount))
  real(rkind)                         :: canopyWatVol_stateVar(sum(gru_struc(:)%hruCount))
  real(rkind)                         :: aquiferWatVol_stateVar(sum(gru_struc(:)%hruCount))
  real(rkind)                         :: sweWatVol_stateVar(sum(gru_struc(:)%hruCount), nSnow_2openwq)
  real(rkind)                         :: soilWatVol_stateVar(sum(gru_struc(:)%hruCount), nSoil_2openwq)
  real(rkind)                         :: soilTemp_K_depVar(sum(gru_struc(:)%hruCount), nSoil_2openwq)
  real(rkind)                         :: soilMoist_depVar(sum(gru_struc(:)%hruCount), nSoil_2openwq)
  integer(i4b)                        :: err

  summaVars: associate(&
      progStruct     => summa1_struc%progStruct             , &
      timeStruct     => summa1_struc%timeStruct             , &
      attrStruct     => summa1_struc%attrStruct               &
  )

  ! Update dependencies and storage volumes
  ! Assemble the data to send to openWQ

  openWQArrayIndex = 0 ! index into the arrays that are being passed to openWQ

  do iGRU = 1, size(gru_struc(:))
      do iHRU = 1, gru_struc(iGRU)%hruCount

        openWQArrayIndex = openWQArrayIndex + 1 

        ! ############################
        ! Update unlayered variables and dependencies 
        ! (1 layer only)
        ! ############################

        ! Tair 
        ! (Summa in K)
        airTemp_K_depVar(openWQArrayIndex) = &
          progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanairTemp)%dat(1)

        ! Vegetation
        ! unit for volume = m3 (summa-to-openwq unit conversions needed)
        ! scalarCanopyWat [kg m-2], so needs to  to multiply by hru area [m2] and divide by water density
        canopyWatVol_stateVar(openWQArrayIndex) = &
          progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyWat)%dat(1) &
          * attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea) / 1000

        ! Aquifer
        ! unit for volume = m3 (summa-to-openwq unit conversions needed)
        ! scalarAquiferStorage [m], so needs to  to multiply by hru area [m2] only
        aquiferWatVol_stateVar(openWQArrayIndex) = &
          progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarAquiferStorage)%dat(1) &
          * attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)
        
        ! ############################
        ! Update layered variables and dependenecies
        ! ############################

        ! Soil
        do ilay = 1, nSoil_2openwq
          
          ! Tsoil
          ! (Summa in K)
          soilTemp_K_depVar(openWQArrayIndex, ilay) = &
            progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerTemp)%dat(ilay)

          soilMoist_depVar(openWQArrayIndex, ilay) = 0     ! TODO: Find the value for this varaibles

          ! Soil
          ! unit for volume = m3 (summa-to-openwq unit conversions needed)
          ! mLayerMatricHead [m], so needs to  to multiply by hru area [m2]
          soilWatVol_stateVar(openWQArrayIndex, ilay) = &
            progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerMatricHead)%dat(ilay) &
            * attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)
    
        enddo

        ! Snow
        do ilay = 1, nSnow_2openwq
          
          ! Snow
          ! unit for volume = m3 (summa-to-openwq unit conversions needed)
          ! scalarSWE [kg m-2], so needs to  to multiply by hru area [m2] and divide by water density
          sweWatVol_stateVar(openWQArrayIndex, ilay) = &
            progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSWE)%dat(ilay) &
            * attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea) / 1000

        enddo

        ! **************  
        ! Fluxes
        !************** 

        ! TODO

        ! Copy the prog structure
        do iVar = 1, size(progStruct%gru(iGRU)%hru(iHRU)%var)
          do iDat = 1, size(progStruct%gru(iGRU)%hru(iHRU)%var(iVar)%dat)
            progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iVar)%dat(iDat) = progStruct%gru(iGRU)%hru(iHRU)%var(iVar)%dat(iDat)
          end do
        end do

      end do
  end do

  ! add the time values to the array
  simtime(1) = timeStruct%var(iLookTIME%iyyy)  ! Year
  simtime(2) = timeStruct%var(iLookTIME%im)    ! month
  simtime(3) = timeStruct%var(iLookTIME%id)    ! hour
  simtime(4) = timeStruct%var(iLookTIME%ih)    ! day
  simtime(5) = timeStruct%var(iLookTIME%imin)  ! minute
  
  err=openWQ_obj%run_time_start(&
        sum(gru_struc(:)%hruCount),             & ! total HRUs
        nSnow_2openwq,                          &
        nSoil_2openwq,                          &
        simtime,                                &
        soilMoist_depVar,                       &                    
        soilTemp_K_depVar,                        &
        airTemp_K_depVar,                         &
        sweWatVol_stateVar,                     &
        canopyWatVol_stateVar,                  &
        soilWatVol_stateVar,                    &
        aquiferWatVol_stateVar)

  ! copy progStruct values to progStruct_timestep_start


  end associate summaVars

end subroutine


subroutine run_space_step(  &
    timeStruct,             &
    fluxStruct,             &
    nGRU)

  USE var_lookup,   only: iLookPROG  ! named variables for state variables
  USE var_lookup,   only: iLookTIME  ! named variables for time data structure
  USE var_lookup,   only: iLookFLUX  ! named varaibles for flux data
  USE globalData,   only: openWQ_obj
  USE data_types,   only: var_dlength,var_i
  USE globalData,   only: gru_struc
  implicit none

  type(var_i),             intent(in)    :: timeStruct 
  type(gru_hru_doubleVec), intent(in)    :: fluxStruct
  integer(i4b),            intent(in)    :: nGRU

  integer(i4b)                           :: hru_index ! needed because openWQ saves hrus as a single array
  integer(i4b)                           :: iHRU      ! variable needed for looping
  integer(i4b)                           :: iGRU      ! variable needed for looping

  integer(i4b)                           :: simtime(5) ! 5 time values yy-mm-dd-hh-min
  integer(i4b)                           :: err
  ! compartment indexes
  integer(i4b)                           :: scalarCanopyWat=0 ! SUMMA Side units: kg m-2
  integer(i4b)                           :: mLayerMatricHead=1 ! SUMMA Side units: m
  integer(i4b)                           :: scalarAquifer=2 ! SUMMA Side units: m
  integer(i4b)                           :: mLayerVolFracWat=3 ! SUMMA Side units: ????
  integer(i4b)                           :: iy_r
  integer(i4b)                           :: iz_r
  integer(i4b)                           :: iy_s
  integer(i4b)                           :: iz_s

  ! Fluxes leaving the canopy
  real(rkind)                            :: scalarCanopySnowUnloading ! kg m-2 s-1
  real(rkind)                            :: scalarCanopyLiqDrainage   ! kg m_2 s-1



  simtime(1) = timeStruct%var(iLookTIME%iyyy)  ! Year
  simtime(2) = timeStruct%var(iLookTIME%im)    ! month
  simtime(3) = timeStruct%var(iLookTIME%id)    ! hour
  simtime(4) = timeStruct%var(iLookTIME%ih)    ! day
  simtime(5) = timeStruct%var(iLookTIME%imin)  ! minute
  
  hru_index=1
  do iGRU=1,nGRU
    do iHRU=1,gru_struc(iGRU)%hruCount
      ! Canopy Fluxes
      scalarCanopySnowUnloading = fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1)
      scalarCanopyLiqDrainage = fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)
      
      iy_s = 1
      iz_s = 1
      iy_r = 1
      iz_r = 1

      err=openwq_obj%run_space(simtime,                                                                         &
                           scalarCanopyWat, hru_index, iy_s, iz_s,                                              &
                           mLayerVolFracWat, hru_index, iy_r, iz_r,                                             &
                           scalarCanopySnowUnloading,                                                           &
                           progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyWat)%dat(1))
      
      err=openwq_obj%run_space(simtime,                                                                         &
                           scalarCanopyWat, hru_index, iy_s, iz_s,                                              &
                           mLayerVolFracWat, hru_index, iy_r, iz_r,                                             &
                           scalarCanopyLiqDrainage,                                                             &
                           progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyWat)%dat(1))



      hru_index = hru_index + hru_index
    end do
  end do



end subroutine run_space_step


subroutine run_time_end( &
  openWQ_obj, &
  summa1_struc)

  USE summa_type, only:summa1_type_dec            ! master summa data type
  
  USE var_lookup, only: iLookTIME  ! named variables for time data structure

  implicit none

  ! Dummy Varialbes
  class(ClassWQ_OpenWQ), intent(in)  :: openWQ_obj
  type(summa1_type_dec), intent(in)  :: summa1_struc

  ! Local Variables
  integer(i4b)                       :: simtime(5) ! 5 time values yy-mm-dd-hh-min
  integer(i4b)                       :: err ! error control

  summaVars: associate(&
      timeStruct     => summa1_struc%timeStruct       &       
  )

  simtime(1) = timeStruct%var(iLookTIME%iyyy)  ! Year
  simtime(2) = timeStruct%var(iLookTIME%im)    ! month
  simtime(3) = timeStruct%var(iLookTIME%id)    ! hour
  simtime(4) = timeStruct%var(iLookTIME%ih)    ! day
  simtime(5) = timeStruct%var(iLookTIME%imin)  ! minute

  err=openwq_obj%run_time_end(simtime)           ! minute

  end associate summaVars
end subroutine



end module summa_openWQ