module update_module

  use datatypes_module
  use grid_module
  use variables_module
  use params_module
  use bcs_module
  
  implicit none

  private

  public :: update

contains
  
  subroutine update(U, vf, fluxes, dt)

    type(gridvar_t),     intent(inout) :: U
    type(gridedgevar_t), intent(in   ) :: fluxes, vf
    real (kind=dp_t),    intent(in   ) :: dt

    type(gridvar_t) :: Uold, U_temp

    real (kind=dp_t) :: dtdx, f(-1:1)
    
    integer :: i, n

    ! store the old state
    call build(Uold, U%grid, U%nvar)
    do i = U%grid%lo-1, U%grid%hi+1
       Uold%data(i,:) = U%data(i,:)
    enddo
    
    ! update
    do i = U%grid%lo-1, U%grid%hi+1
       U%data(i,:) = U%data(i,:) + &
            (dt/U%grid%dx)*(fluxes%data(i,:) - fluxes%data(i+1,:))
    enddo    
    
    ! If we're using the Lagrange + remap method, do the remap now.
    ! At present this is only intended to work with periodic BCs.

    if (invariant_hydro .eq. 2) then

       ! Build a temporary copy of the data.

       call build(U_temp, U%grid, U%nvar)
       do i = U%grid%lo, U%grid%hi
          U_temp%data(i,:) = U%data(i,:)
       enddo

       call fillBCs(U_temp)

       ! Remap assumes uniform grid spacing, but it's straightforward to adjust
       ! for non-uniform spacing.

       dtdx = dt / U%grid%dx

       do i = U%grid%lo, U%grid%hi

          ! If the left interface is moving to the left,
          ! we get no contribution from the left zone;
          ! otherwise, we get v * dt / dx.

          f(-1) =  vf%data(i  ,1) * dtdx

          ! If the right interface is moving to the right,
          ! we get no contribution from the right zone;
          ! otherwise, we get v * dt / dx.

          f( 1) = -vf%data(i+1,1) * dtdx

          ! Since we're conservative, the contribution from the
          ! center zone is just the original zone minus
          ! the fractions of the left and right zones.

          f(0) = ONE - max(ZERO,f(-1)) - max(ZERO,f(1))

          f = max(ZERO,f)
                    
       
          if (remap_order .eq. 0) then

             do n = 1, 3
                U % data(i,n) = sum(f * U_temp % data(i-1:i+1,n))
             enddo

          else

             print *, "Cannot handle remap_order > 1 at present."
             stop

          endif

       enddo

       call destroy(U_temp)

    endif
       
    ! time-centered source terms 
    do i = U%grid%lo, U%grid%hi
       U%data(i,iumomx) = U%data(i,iumomx) + &
            0.5_dp_t*dt*(U%data(i,iudens) + Uold%data(i,iudens))*grav

       U%data(i,iuener) = U%data(i,iuener) + &
            0.5_dp_t*dt*(U%data(i,iumomx) + Uold%data(i,iumomx))*grav
    enddo

    ! clean-up
    call destroy(Uold)

  end subroutine update

end module update_module
