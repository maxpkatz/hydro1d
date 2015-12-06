module datatypes_module

  implicit none

  integer, parameter :: dp_t = selected_real_kind(15,307)

  real (kind=dp_t), parameter :: ZERO  = 0.0_dp_t
  real (kind=dp_t), parameter :: HALF  = 0.5_dp_t
  real (kind=dp_t), parameter :: ONE   = 1.0_dp_t
  real (kind=dp_t), parameter :: TWO   = 2.0_dp_t
  real (kind=dp_t), parameter :: THREE = 3.0_dp_t
  real (kind=dp_t), parameter :: FOUR  = 4.0_dp_t
  real (kind=dp_t), parameter :: SIX   = 6.0_dp_t

  real (kind=dp_t), parameter :: FOURTH = ONE / FOUR
  real (kind=dp_t), parameter :: TWO3RD = TWO / THREE
  real (kind=dp_t), parameter :: SIXTH  = ONE / SIX
  
end module datatypes_module
