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
      nCanopy_2openwq,         & ! num layers of canopy (fixed to 1)
      nSnow_2openwq,           & ! num layers of snow (fixed to max of 5 because it varies)
      nSoil_2openwq,           & ! num layers of snoil (variable)
      nAquifer_2openwq,        & ! num layers of aquifer (fixed to 1)
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

   integer function openWQ_run_time_start(this, numHRU, maxNumLayers_snow, maxNumLayers_soil, simtime, &
      soilMoisture, soilTemp, airTemp, swe_vol, canopyWat_vol, matricHead_vol, aquiferStorage_vol)
      
      implicit none
      class(ClassWQ_OpenWQ)      :: this
      integer(i4b), intent(in)   :: numHRU
      integer(i4b), intent(in)   :: maxNumLayers_snow
      integer(i4b), intent(in)   :: maxNumLayers_soil
      integer(i4b), intent(in)   :: simtime(5) ! 5 is the number of timevars
      real(rkind),  intent(in)   :: soilMoisture(numHRU, maxNumLayers_soil)
      real(rkind),  intent(in)   :: soilTemp(numHRU, maxNumLayers_soil)
      real(rkind),  intent(in)   :: airTemp(numHRU)
      real(rkind),  intent(in)   :: swe_vol(numHRU, maxNumLayers_snow)
      real(rkind),  intent(in)   :: canopyWat_vol(numHRU)
      real(rkind),  intent(in)   :: matricHead_vol(numHRU,maxNumLayers_soil)
      real(rkind),  intent(in)   :: aquiferStorage_vol(numHRU)

      openWQ_run_time_start = openwq_run_time_start_c(this%ptr, numHRU, maxNumLayers_snow, maxNumLayers_soil, simtime, &
         soilMoisture, soilTemp, airTemp, swe_vol, canopyWat_vol, matricHead_vol, aquiferStorage_vol)
   
      end function

   integer function openWQ_run_space(this,simtime,source,ix_s,iy_s,iz_s, &
         recipient,ix_r,iy_r,iz_r,wflux_s2r,wmass_source)
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

      openWQ_run_space = openwq_run_space_c(this%ptr,simtime,source,ix_s,iy_s,iz_s,recipient,ix_r,iy_r,iz_r,wflux_s2r,wmass_source)
   
   end function

   integer function openWQ_run_space_in(this,simtime,recipient,ix_r,iy_r,iz_r,wflux_s2r)
      implicit none
      class(ClassWQ_OpenWQ)      :: this
      integer(i4b), intent(in)   :: simtime(5) ! 5 is the number of timevars
      integer(i4b), intent(in)   :: recipient
      integer(i4b), intent(in)   :: ix_r
      integer(i4b), intent(in)   :: iy_r
      integer(i4b), intent(in)   :: iz_r
      real(rkind),  intent(in)   :: wflux_s2r

      openWQ_run_space_in = openwq_run_space_in_c(this%ptr,simtime,recipient,ix_r,iy_r,ix_r,wflux_s2r)

   end function


   integer function openWQ_run_time_end(this, simtime)
      implicit none
      class(ClassWQ_OpenWQ)      :: this
      integer(i4b), intent(in)   :: simtime(5) ! 5 is the number of timevars

      openWQ_run_time_end = openWQ_run_time_end_c(this%ptr,simtime)

   end function

end module openwq