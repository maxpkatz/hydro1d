module init_module

  use datatypes_module
  use probparams_module
  use grid_module
  use eos_module
  use variables_module

  implicit none

  private

  public :: init_data

contains
  
  subroutine init_data(U)

    type(gridvar_t), intent(inout) :: U

    integer :: i
    real (kind=dp_t) :: rho, ux, p, e

    real (kind=dp_t) :: xcenter

    xcenter = 0.5_dp_t*(U%grid%xmin + U%grid%xmax)

    do i = U%grid%lo, U%grid%hi

       if (abs(U%grid%x(i) - xcenter) < width) then

          ! ambient state
          rho = rho_t
          ux  = v_adv
          p   = p_t
          
       else

          ! tophat state
          rho = rho_l
          ux  = v_adv
          p   = p_l

       endif

       ! get the internal energy from the EOS
       call eos(eos_input_p, p, e, rho)
       
       ! store the conserved state
       U%data(i,iudens) = rho
       U%data(i,iumomx) = rho*ux
       U%data(i,iuener) = 0.5_dp_t*rho*ux*ux + rho*e
          
    enddo

  end subroutine init_data

end module init_module

