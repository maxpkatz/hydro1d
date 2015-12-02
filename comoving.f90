module interface_velocities_module

  use datatypes_module
  use grid_module
  use variables_module
  use params_module
  
  implicit none

  private

  public :: make_interface_velocities

contains
  
  subroutine make_interface_velocities(U, vf)

    type(gridvar_t),     intent(in   ) :: U
    type(gridedgevar_t), intent(inout) :: vf
    
    integer :: i

    do i = U%grid%lo-1, U%grid%hi+2
       if (invariant_hydro) then
          vf%data(i,1) = 0.5 * (U%data(i,2) / U%data(i,1) + U%data(i-1,2) / U%data(i-1,1))
       else
          vf%data(i,1) = 0.0
       endif
    enddo

  end subroutine make_interface_velocities

end module interface_velocities_module
