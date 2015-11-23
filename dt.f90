! dt module
!
! compute the timestep as the minimum dx / (|U| + cs)

module dt_module

  use datatypes_module
  use grid_module
  use params_module
  use eos_module
  use variables_module

  implicit none

  private

  public :: compute_dt

contains

  subroutine compute_dt(U, vf, n, dt)

    type(gridvar_t),     intent(in   ) :: U
    type(gridedgevar_t), intent(in   ) :: vf
    integer,             intent(inout) :: n
    real (kind=dp_t),    intent(inout) :: dt

    integer :: i
    real (kind=dp_t) :: cs, p, e, rho

    
    dt = huge(0.0_dp_t)

    do i = U%grid%lo, U%grid%hi

       ! compute cs (soundspeed)
       e = (U%data(i,iuener) - 0.5_dp_t*U%data(i,iumomx)**2/U%data(i,iudens)) / &
            U%data(i,iudens)
       
       rho = U%data(i,iudens)
       call eos(eos_input_e, p, e, rho)

       cs = sqrt(gamma*p/U%data(i,iudens))

       ! If the zone faces move in a Lagrangian manner, then the timestep restriction
       ! is the smaller of dx / cs and dx / vf.

       if (lagrange_remap) then
          dt = min(dt, U%grid%dx / cs, U%grid%dx / maxval(abs(vf%data(i:i+1,1))))
       else
          dt = min(dt, U%grid%dx/(abs(U%data(i,iumomx)/U%data(i,iudens)) + cs))
       endif
       
       
    enddo    

    dt = cfl*dt

    if (n == 0) dt = init_shrink*dt

  end subroutine compute_dt

end module dt_module
