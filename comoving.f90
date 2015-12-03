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

    real(kind=dp_t) :: v_l, v_r

    vf % data(:,:) = ZERO
    
    if (invariant_hydro) then
    
       do i = U%grid%lo-U%grid%ng+1, U%grid%hi+U%grid%ng

          ! Edge based indexing
          v_l = U % data(i-1,2) / U % data(i-1,1)
          v_r = U % data(i  ,2) / U % data(i  ,1)
          
          vf % data(i,1) = HALF * (v_l + v_r)
          
       enddo

    endif
       
  end subroutine make_interface_velocities

end module interface_velocities_module
