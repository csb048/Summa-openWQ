module openwq
   
 USE, intrinsic :: iso_c_binding
 USE nrtype
 private
 public :: ClassWQ_OpenWQ

 include "openWQInterface.f90"

 type ClassWQ_OpenWQ
    private
    type(c_ptr) :: ptr ! pointer to openWQ class

 contains
   !  procedure :: get_num => openWQ_get_num
    procedure :: decl => openWQ_init
    procedure :: run_time_start => openWQ_run_time_start
    procedure :: run_space => openWQ_run_space
    procedure :: run_space_in => openWQ_run_space_in
    procedure :: run_time_end => openWQ_run_time_end

 end type

 interface ClassWQ_OpenWQ
    procedure create_openwq
 end interface
 contains
    function create_openwq()
        implicit none
        type(ClassWQ_OpenWQ) :: create_openwq
        create_openwq%ptr = create_openwq_c()
    end function

    ! supposed to be decl but needed to openWQ_decl in the interface file
    ! returns integer of either a failure(-1) or success(0)
   integer function openWQ_init( &
      this,                      & ! openwq object
      num_hru,                   & ! num HRU
      nCanopy_2openwq,           & ! num layers of canopy (fixed to 1)
      nSnow_2openwq,             & ! num layers of snow (fixed to max of 5 because it varies)
      nSoil_2openwq,             & ! num layers of snoil (variable)
      nAquifer_2openwq,          & ! num layers of aquifer (fixed to 1)
      nYdirec_2openwq)                 ! num of layers in y-dir (set to 1 because not used in summa)
      
      implicit none
      class(ClassWQ_OpenWQ) :: this
      integer(i4b), intent(in) :: num_hru
      integer(i4b), intent(in) :: nCanopy_2openwq
      integer(i4b), intent(in) :: nSnow_2openwq
      integer(i4b), intent(in) :: nSoil_2openwq
      integer(i4b), intent(in) :: nAquifer_2openwq
      
      integer(i4b), intent(in) :: nYdirec_2openwq

      openWQ_init = openwq_decl_c(  &
         this%ptr,                  & ! openwq object
         num_hru,                   & ! num HRU
         nCanopy_2openwq,         & ! num layers of canopy (fixed to 1)
         nSnow_2openwq,           & ! num layers of snow (fixed to max of 5 because it varies)
         nSoil_2openwq,           & ! num layers of snoil (variable)
         nAquifer_2openwq,        & ! num layers of aquifer (fixed to 1)
         nYdirec_2openwq)                 ! num of layers in y-dir (set to 1 because not used in summa)

    end function
!  ! Globaly accessible variable

   integer function openWQ_run_time_start(   &
      this,                                  &
      numHRU,                                &
      nSnow_2openwq,                         &
      nSoil_2openwq,                         &
      simtime,                               &
      soilMoist_depVar,                      &
      soilTemp_depVar,                       &
      airTemp_K_depVar,                        &
      sweWatVol_stateVar,                    &
      canopyWatVol_stateVar,                 &
      soilWatVol_stateVar,                   &
      aquiferWatVol_stateVar)
      
      implicit none
      class(ClassWQ_OpenWQ)      :: this
      integer(i4b), intent(in)   :: numHRU
      integer(i4b), intent(in)   :: nSnow_2openwq
      integer(i4b), intent(in)   :: nSoil_2openwq
      integer(i4b), intent(in)   :: simtime(5) ! 5 is the number of timevars
      real(rkind),  intent(in)   :: airTemp_K_depVar(numHRU)
      real(rkind),  intent(in)   :: soilTemp_depVar(numHRU, nSoil_2openwq)
      real(rkind),  intent(in)   :: soilMoist_depVar(numHRU, nSoil_2openwq)
      real(rkind),  intent(in)   :: canopyWatVol_stateVar(numHRU)
      real(rkind),  intent(in)   :: sweWatVol_stateVar(numHRU, nSnow_2openwq)
      real(rkind),  intent(in)   :: soilWatVol_stateVar(numHRU,nSoil_2openwq)
      real(rkind),  intent(in)   :: aquiferWatVol_stateVar(numHRU)

      openWQ_run_time_start = openwq_run_time_start_c( &
         this%ptr,               & 
         numHRU,                 &
         nSnow_2openwq,          &
         nSoil_2openwq,          &
         simtime,                &
         soilMoist_depVar,       &
         soilTemp_depVar,        &
         airTemp_K_depVar,         &
         sweWatVol_stateVar,     &
         canopyWatVol_stateVar,  &
         soilWatVol_stateVar,    &
         aquiferWatVol_stateVar)
   
      end function

   integer function openWQ_run_space(  &
      this,                            &
      simtime,                         &
      source,ix_s,iy_s,iz_s,           &
      recipient,ix_r,iy_r,iz_r,        &
      wflux_s2r,wmass_source)

      implicit none
      class(ClassWQ_OpenWQ)      :: this
      integer(i4b), intent(in)   :: simtime(5) ! 5 is the number of timevars
      integer(i4b), intent(in)   :: source
      integer(i4b), intent(in)   :: ix_s
      integer(i4b), intent(in)   :: iy_s
      integer(i4b), intent(in)   :: iz_s
      integer(i4b), intent(in)   :: recipient
      integer(i4b), intent(in)   :: ix_r
      integer(i4b), intent(in)   :: iy_r
      integer(i4b), intent(in)   :: iz_r
      real(rkind),  intent(in)   :: wflux_s2r
      real(rkind),  intent(in)   :: wmass_source

      openWQ_run_space = openwq_run_space_c( &
         this%ptr,                           &
         simtime,                            &
         source,ix_s,iy_s,iz_s,              &
         recipient,ix_r,iy_r,iz_r,           &
         wflux_s2r,wmass_source)
   
   end function

   integer function openWQ_run_space_in(  &
      this,                               &
      simtime,                            &
      recipient,ix_r,iy_r,iz_r,           &
      wflux_s2r)

      implicit none
      class(ClassWQ_OpenWQ)      :: this
      integer(i4b), intent(in)   :: simtime(5) ! 5 is the number of timevars
      integer(i4b), intent(in)   :: recipient
      integer(i4b), intent(in)   :: ix_r
      integer(i4b), intent(in)   :: iy_r
      integer(i4b), intent(in)   :: iz_r
      real(rkind),  intent(in)   :: wflux_s2r

      openWQ_run_space_in = openwq_run_space_in_c( &
         this%ptr,                                 &
         simtime,                                  &
         recipient,ix_r,iy_r,ix_r,                 &
         wflux_s2r)

   end function


   integer function openWQ_run_time_end(  &
      this,                               &
      simtime)

      implicit none
      class(ClassWQ_OpenWQ)      :: this
      integer(i4b), intent(in)   :: simtime(5) ! 5 is the number of timevars

      openWQ_run_time_end = openWQ_run_time_end_c( &
         this%ptr,                                 &
         simtime)

   end function

end module openwq