module CropHeatStress

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !MODULE: CropHeatStress
  !
  ! !DESCRIPTION:
  ! Module for implementing the effect of heat stress on crop production
  ! Added by SdR. 
  !
  ! !USES:
  !------------------------------------------------------------------------
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use shr_sys_mod       , only : shr_sys_flush
  use shr_infnan_mod    , only : isnan => shr_infnan_isnan, isinf => shr_infnan_isinf, nan => shr_infnan_nan, assignment(=)
  use clm_varcon, only : spval
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: crop_heatstress_ndays       ! SdR: checks for number of days above tcrit for crop heat stress
  public  :: calc_HS_factor              ! SdR: calculates heat stress magnitude affecting leaf area decline (grainfill phase)
  public  :: crop_heatstress_reset       ! Unsets variables related to crop heat stress
  public  :: calc_TVDAY_peak
  public  :: check_min_TVpeak_years

  !S
  ! !PUBLIC FOR UNIT TESTING
  real(r8), public, parameter :: tcrit_min = 296.15_r8
  real(r8), public, parameter :: HS_ndays_min = 3._r8

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains

  subroutine crop_heatstress_reset(HS_ndays, heatwave_crop)

    ! !DESCRIPTION
    ! Unsets variables related to crop heat stress

    ! !ARGUMENTS:
    real(r8), intent(inout)    :: HS_ndays ! number of crop heat stress days (ndays) should be integer at final implementation
    real(r8), intent(inout)    :: heatwave_crop ! keep track if heatwave is activated

    HS_ndays = 0._r8
    heatwave_crop = 0._r8

  end subroutine crop_heatstress_reset

  !------------------------------------------------------------------------
  subroutine crop_heatstress_ndays(HS_ndays, heatwave_crop, t_veg_day, croplive, peakTVDAY_years)

    ! !DESCRIPTION:
    ! added by SdR for heat stress implementation
    ! function to keep track of critical temperature for crops that is exceeded at the end of each day
    ! needs a minimum of 3 consecutive days before heat stress is activated

    ! !ARGUMENTS:
    real(r8),        intent(inout)    :: HS_ndays              ! number of crop heat stress days (ndays) should be integer at final implementation
    real(r8),        intent(inout)    :: heatwave_crop         ! keep track if heatwave is activated
    real(r8),        intent(in)       :: t_veg_day
    real(r8),        intent(in)       :: peakTVDAY_years       ! peak vegetation temperature, minimum over simulated year (Kelvin)
    logical,         intent(in)       :: croplive              ! crop between sowing and harvest

    ! !LOCAL VARIABLES:
    real(r8)   :: tcrit


    !----------------------------------------------------------------------
    tcrit  = peakTVDAY_years

    if (tcrit < tcrit_min) then
      tcrit = tcrit_min
    end if

    ! No heat stress if crop isn't alive
    if (.not. croplive) then
       call crop_heatstress_reset(HS_ndays, heatwave_crop)
       return
    end if

    ! Don't do anything if it's not a real temperature
    if (t_veg_day > spval / 1000._r8) then
       return
    end if

    ! check if tcrit is exceeded and count days
    if (t_veg_day >= tcrit .and. t_veg_day >= tcrit_min) then
         HS_ndays = HS_ndays + 1.0_r8
    else
         HS_ndays = 0.0_r8
    end if

    ! check if a heatwave is occurring
    if (HS_ndays >= (HS_ndays_min - 0.2_r8)) then
         heatwave_crop = 1.0_r8 
    else
         heatwave_crop = 0.0_r8
    end if

  end subroutine crop_heatstress_ndays


  subroutine calc_HS_factor(HS_factor, HS_ndays, t_veg_day, croplive,peakTVDAY_years)

    ! !DESCRIPTION:
    ! function to calculate heat stress instensity by applying a factor to bglfr (increasing LAI decline). function based on Apsim-Nwheat model Asseng et al. 2011:https://doi.org/10.1111/j.1365-2486.2010.02262.x
    ! WIP: different tcrit, tmax values for different crops or climatologies

    ! !ARGUMENTS:
    real(r8),        intent(inout)     :: HS_factor         ! keep track if heatwave is activated
    real(r8),        intent(in)        :: HS_ndays          ! number of crop heat stress days (ndays) should be integer at final implementation
    real(r8),        intent(in)        :: t_veg_day         ! daily vegetation temperature (Kelvin)
    real(r8),        intent(in)        :: peakTVDAY_years   ! peak vegetation temperature, minimum over simulated year (Kelvin)
    logical,         intent(in)        :: croplive          ! crop between sowing and harvest

    ! !LOCAL VARIABLES:
    integer  :: day_min
    real(r8) :: tcrit, tmax, Fheat_max, onset_jump

    !-----------------------------------------------------------------------

    day_min = 3
    tcrit  = peakTVDAY_years

    ! define Tmax based on Tcrit value. Larger values for Tmax results in steeper slopes to Tmax. Tmax range is between 35 and 49 degrees Celsius
    if (tcrit < tcrit_min) then
      tcrit = tcrit_min
      tmax  = 273.15_r8 + 35._r8
    else if (tcrit > tcrit_min .and. tcrit <= 318.15_r8) then
      !Tcrit smaller or eq to 45degreesC
      tmax = (273.15_r8 + 35._r8) + (7._r8 / 11._r8) * (tcrit - tcrit_min)
    else if (tcrit > 318.15_r8) then
      tcrit = 318.15_r8
      tmax = 273.15_r8 + 49._r8
    end if

    ! function parameters to be tested
    Fheat_max  = 0.6_r8   ! lai(5, 15, 25) rep(0.3,0.6, 0.9)
    onset_jump = 0.2_r8 * Fheat_max


    !check  if stress occurs
    if (HS_ndays >= day_min .and. HS_ndays < (day_min + 1) .and. croplive .and. t_veg_day > tcrit .and. t_veg_day > tcrit_min) then
      ! onset heatwave first day
      HS_factor = onset_jump
    else if (HS_ndays > day_min .and. croplive .and. t_veg_day > tcrit .and. t_veg_day > tcrit_min) then
      if (t_veg_day <= tmax ) then
        HS_factor = onset_jump + ((Fheat_max - onset_jump) * ((t_veg_day - tcrit) / (tmax - tcrit)))
      else
        HS_factor = Fheat_max 
      end if
    else
      HS_factor = 0._r8  ! no stress values lai:1._r8 ; rep:0._r8
    end if


  end subroutine calc_HS_factor

  subroutine calc_TVDAY_peak(peakTVDAY, t_veg_day, croplive)
    ! !DESCRIPTION:
    ! Keeps track of maximum vegetation temperature when crop is alive
    ! !ARGUMENTS:
    real(r8),        intent(inout)     :: peakTVDAY   ! peak vegetation temperature during crop growing season (Kelvin)
    real(r8),        intent(in)        :: t_veg_day   ! daily vegetation temperature (Kelvin)
    logical,         intent(in)        :: croplive    ! crop between sowing and harvest

    if (croplive .and. t_veg_day > peakTVDAY) then ! add .and. ((crop_phase == cphase_leafemerge) .or. (crop_phase == cphase_grainfill))
      peakTVDAY = t_veg_day
    end if

  end subroutine calc_TVDAY_peak

  subroutine check_min_TVpeak_years(peakTVDAY_years,peakTVDAY, npeakyears)
    ! !DESCRIPTION:
    ! called at the end of the growing season, which is assumed to be once a year
    ! Keeps track of mminimum TVDAY_peak during simulation years: maximum daytime vegetation
    ! temperature during growing season that occurs at least once every year
    ! also resets peakTVDAY to initial value, to start calc_TVDAY_peak in the new growing season from scratch
    ! !ARGUMENTS:
    real(r8),        intent(inout)  :: peakTVDAY_years  ! peak daily vegetation temperature over simulation years (Kelvin)
    real(r8),        intent(inout)  :: peakTVDAY        ! peak daily vegetation temperature during crop growing season (Kelvin)
    integer,         intent(inout)  :: npeakyears       ! keep track of number of updates only update for first 5 years
        
    ! !LOCAL VARIABLES:
    integer  :: max_nyears

    !-----------------------------------------------------------------------
    ! Initialized at 0, first year is 1 
    max_nyears = 5
    npeakyears = npeakyears + 1

    if (peakTVDAY > 2._r8 .and. npeakyears == 1) then
      peakTVDAY_years = peakTVDAY
    else if (peakTVDAY > 2._r8 .and. npeakyears <= max_nyears) then
      if (peakTVDAY_years > peakTVDAY .or. peakTVDAY_years < 2._r8) then
        peakTVDAY_years = peakTVDAY
      else
        peakTVDAY_years = peakTVDAY_years
      end if
    end if

    ! re-initialize peakTVDAY for new season
    peakTVDAY = 1._r8

  end subroutine check_min_TVpeak_years


end module CropHeatStress
