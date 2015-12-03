module interface_states_godunov_module

  use datatypes_module
  use grid_module
  use variables_module

  implicit none

  private

  public :: make_interface_states_godunov

contains
  
  subroutine make_interface_states_godunov(U, U_l, U_r, vf, dt)

    type(gridvar_t),     intent(in   ) :: U
    type(gridedgevar_t), intent(inout) :: U_l, U_r, vf
    real (kind=dp_t),    intent(in   ) :: dt

    integer :: i

    real (kind=dp_t) :: rho_e_l, rho_e_r

    ! piecewise constant slopes

    ! loop over interfaces and fill the states 
    do i = U%grid%lo, U%grid%hi+1
       U_l%data(i,:) = U%data(i-1,:)
       U_r%data(i,:) = U%data(i,:)

       ! subtract off the face velocities
       U_l%data(i,2) = U_l%data(i,2) - U%data(i-1,1) * vf%data(i,1)
       U_r%data(i,2) = U_r%data(i,2) - U%data(i  ,1) * vf%data(i,1)


       ! Subtract off corresponding kinetic energies. We do this by saving the
       ! old (rho * e) and adding to it the kinetic energy in the co-moving frame.

       rho_e_l = U%data(i-1,3) - 0.5_dp_t * U%data(i-1,2)**2 / U%data(i-1,1)
       rho_e_r = U%data(i  ,3) - 0.5_dp_t * U%data(i  ,2)**2 / U%data(i  ,1)

       U_l%data(i,3) = rho_e_l + 0.5_dp_t * U_l%data(i,2)**2 / U_l%data(i,1)
       U_r%data(i,3) = rho_e_r + 0.5_dp_t * U_r%data(i,2)**2 / U_r%data(i,1)
    enddo

  end subroutine make_interface_states_godunov

end module interface_states_godunov_module
