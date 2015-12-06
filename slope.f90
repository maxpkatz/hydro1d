module slope_module

  use datatypes_module

  implicit none

contains

  ! 4th order MC limiting.  See Colella (1985) Eq. 2.5 and 2.6,
  ! Colella (1990) page 191 (with the delta a terms all equal) or
  ! Saltzman 1994, page 156    

  function slope(U, dx) result(dUdx)

    implicit none

    real(kind=dp_t) :: U(-2:2), dx(-2:2)
    
    real(kind=dp_t) :: dUvl_l, dUvl_r

    real(kind=dp_t) :: dU(-1:1), dUdx

    dU(-1) = (U(0) - U(-1)) / (HALF * (dx(0) + dx(-1)))
    du( 0) = (U(1) - U(-1)) / (dx(0) + HALF * (dx(1) + dx(-1)))
    du( 1) = (U(1) - U( 0)) / (HALF * (dx(1) + dx( 0)))

    dUvl_l = van_leer_slope(U(-2:0), dx(-2:0))
    dUvl_r = van_leer_slope(U( 0:2), dx( 0:2))

    if (dU(-1) * dU(1) > ZERO) then

       dUdx = min( TWO3RD * abs(TWO * dU(0) - FOURTH * (dUvl_l + dUvl_r)), &
                   TWO * abs(dU(1)), TWO * abs(dU(-1)) ) * sign(ONE, dU(0))

    else

       dUdx = ZERO

    endif

  end function slope


  ! 2nd order MC limiting.

  function van_leer_slope(U, dx) result(dUdx)

    implicit none

    real(kind=dp_t) :: U(-1:1), dx(-1:1)

    real(kind=dp_t) :: dU(-1:1), dUdx

    dU(-1) = (U(0) - U(-1)) / (HALF * (dx(0) + dx(-1)))
    du( 0) = (U(1) - U(-1)) / (dx(0) + HALF * (dx(1) + dx(-1)))
    du( 1) = (U(1) - U( 0)) / (HALF * (dx(1) + dx( 0)))

          
    if (dU(-1) * dU(1) > ZERO) then

       dUdx = min(abs(dU(0)), TWO * abs(dU(1)), TWO * abs(dU(-1))) * sign(ONE,dU(0))

    else

       dUdx = ZERO

    endif

  end function van_leer_slope

end module slope_module
