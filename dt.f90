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
    real (kind=dp_t) :: cs, p, e, rho, mom, vel
    
    dt = huge(0.0_dp_t)

    do i = U%grid%lo, U%grid%hi

       rho = U%data(i,iudens)
       mom = U%data(i,iumomx)
       vel = mom / rho
       
       ! compute cs (soundspeed)
       e = (U%data(i,iuener) - 0.5_dp_t*rho*vel**2) / rho
       
       call eos(eos_input_e, p, e, rho)

       cs = sqrt(gamma*p/rho)

       dt = min(dt, U%grid%dx/(abs(vel) + cs))
       
    enddo    

    dt = cfl*dt

    if (n == 0) dt = init_shrink*dt

  end subroutine compute_dt

end module dt_module
