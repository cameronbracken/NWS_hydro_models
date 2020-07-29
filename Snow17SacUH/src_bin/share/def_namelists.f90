! A.Wood-2016:  large edits for distributed SAC
module def_namelists
  use nrtype
  implicit none

  ! initialization variables
  ! these are in the namelist bloc &INIT_CONTROL
  character(len = 20)   :: main_id              ! ID string used for main/combined output
  integer(I4B)		:: n_hrus		! number of HRU areas in parameter files
  character(len = 1024) :: forcing_root		! base name of forcing data file root
  character(len = 1024) :: output_root		! base name for output files
  character(len = 1024) :: sac_param_file	! name for sac-sma parameters
  character(len = 1024) :: snow17_param_file	! name for snow17 parameters
  character(len = 1024) :: uh_param_file	! name for uh and pet parameters
  character(len = 1024) :: forcing_adjust_file  ! file with forcing adjustments (only used if use_forcing_adjust=1)
  character(len = 1024)	:: uh_state_out_root	! name for uh state output root
  character(len = 1024)	:: uh_state_in_root	! name for uh state input root
  character(len = 1024)	:: snow_state_out_root  ! name for snow17 state output root
  character(len = 1024)	:: sac_state_out_root	! name for sac state output root
  character(len = 1024)	:: snow_state_in_root	! name for snow17 state input root
  character(len = 1024)	:: sac_state_in_root	! name for sac state input root
  integer(I4B)		:: output_hrus 		! output HRU results? (1=yes; 0=no)

  integer(I4B)    :: dt    ! model timestep in seconds
  integer(I4B)    :: calibration_mode    ! calibration mode turns off most output except routed flow or runoff
  integer(I4B)    :: use_forcing_adjust  ! turn on or off input of a file "forcing_adjust_file"

  integer(I4B)		:: start_month		! starting month 
  integer(I4B)    :: start_hour    ! starting hour
  integer(I4B)		:: start_day		! starting day
  integer(I4B)		:: start_year		! starting year
  integer(I4B)		:: end_month		! ending month 
  integer(I4B)		:: end_day		! ending day
  integer(I4B)    :: end_hour    ! ending day
  integer(I4B)		:: end_year		! ending year
  integer(I4B)		:: warm_start_run	! warm restart run flag
  integer(I4B)		:: write_states	        ! flag to write states for a warm start run

  real(sp)		:: init_swe             ! initial states if for a cold start run
  real(sp)		:: init_uztwc           ! used in all model HRUs
  real(sp)		:: init_uzfwc           ! model state variables not listed start at 0
  real(sp)		:: init_lztwc
  real(sp)		:: init_lzfsc
  real(sp)		:: init_lzfpc
  real(sp)		:: init_adimc

  ! SAC_model params & other key inputs in the sace param file
  character(len = 20), dimension(:), allocatable :: hru_id   ! local hru id
  real(dp), dimension(:), allocatable :: hru_area   ! sq-km, needed for combination & routing conv.
  real(sp), dimension(:), allocatable :: uztwm,uzfwm,uzk,pctim,adimp,zperc,rexp
  real(sp), dimension(:), allocatable :: lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree
  real(sp), dimension(:), allocatable :: riva,side,rserv,peadj,pxadj

  ! UH params in the UH param file
  real(sp), dimension(:), allocatable :: unit_shape,unit_scale  !unit hydrograph parameters

  ! Snow17_model params in the snow param file
  real(dp), dimension(:), allocatable :: latitude   ! decimal degrees
  real(dp), dimension(:), allocatable :: elev       ! m
  real(sp), dimension(:), allocatable :: scf,mfmax,mfmin,uadj,si,pxtemp
  real(sp), dimension(:), allocatable :: nmf,tipm,mbase,plwhc,daygm
  real(sp), dimension(11)  :: adc  ! different for each hru
  real(dp), dimension(:), allocatable :: adc_a, adc_b, adc_c ! areal depletion curve parameters ax^b+(1-a)x^c

  ! forcing adjustment parameters for map, mat, pet, ptps
  double precision, dimension(:), allocatable:: map_adj_jan, map_adj_feb, map_adj_mar, &
    map_adj_apr, map_adj_may, map_adj_jun, map_adj_jul, map_adj_aug, &
    map_adj_sep, map_adj_oct, map_adj_nov, map_adj_dec, &
    mat_adj_jan, mat_adj_feb, mat_adj_mar, &
    mat_adj_apr, mat_adj_may, mat_adj_jun, mat_adj_jul, mat_adj_aug, &
    mat_adj_sep, mat_adj_oct, mat_adj_nov, mat_adj_dec, &
    pet_adj_jan, pet_adj_feb, pet_adj_mar, &
    pet_adj_apr, pet_adj_may, pet_adj_jun, pet_adj_jul, pet_adj_aug, &
    pet_adj_sep, pet_adj_oct, pet_adj_nov, pet_adj_dec, &
    ptps_adj_jan, ptps_adj_feb, ptps_adj_mar, &
    ptps_adj_apr, ptps_adj_may, ptps_adj_jun, ptps_adj_jul, ptps_adj_aug, &
    ptps_adj_sep, ptps_adj_oct, ptps_adj_nov, ptps_adj_dec

  ! namelist elements to be shared
  namelist / INIT_CONTROL / forcing_root, output_root, main_id, n_hrus, output_hrus, &
                          start_hour,start_day,start_month,start_year,end_year,end_month, &
                          end_hour,end_day,dt,init_swe,init_uztwc,init_uzfwc,init_lztwc,init_lzfsc, &
                          init_lzfpc,init_adimc,sac_param_file,snow17_param_file,uh_param_file, &
                          forcing_adjust_file, &
                          uh_state_in_root, warm_start_run, write_states, &
			  snow_state_out_root,sac_state_out_root,snow_state_in_root, &
			  sac_state_in_root, uh_state_out_root, calibration_mode, use_forcing_adjust
  save
end module
