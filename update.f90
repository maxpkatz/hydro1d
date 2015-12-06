module update_module

  use datatypes_module
  use grid_module
  use variables_module
  use params_module
  use bcs_module
  use slope_module
  
  implicit none

  private

  public :: update

contains
  
  subroutine update(U, vf, fluxes, dt)

    type(gridvar_t),     intent(inout) :: U
    type(gridedgevar_t), intent(in   ) :: fluxes, vf
    real (kind=dp_t),    intent(in   ) :: dt

    type(gridvar_t) :: Uold, U_temp

    real (kind=dp_t) :: dtdx
    real (kind=dp_t) :: U_av(-1:1, U%nvar), f(-1:1)
    real (kind=dp_t) :: ul(-2:2), uc(-2:2), ur(-2:2)
    real (kind=dp_t) :: dxl(-2:2), dxc(-2:2), dxr(-2:2)
    real (kind=dp_t) :: dul, duc, dur
    real (kind=dp_t) :: dx
    
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

       U_temp = U

       call fillBCs(U_temp)

       ! Remap assumes uniform grid spacing, but it's straightforward to adjust
       ! for non-uniform spacing.

       dx = U % grid % dx

       dtdx = dt / dx

       do i = U%grid%lo, U%grid%hi

          ! If the left interface is moving to the left,
          ! we get no contribution from the left zone;
          ! otherwise, we get v * dt / dx.

          f(-1) = max(ZERO,  vf%data(i  ,1) * dtdx)

          ! If the right interface is moving to the right,
          ! we get no contribution from the right zone;
          ! otherwise, we get -v * dt / dx.

          f( 1) = max(ZERO, -vf%data(i+1,1) * dtdx)

          ! Since we're conservative, the contribution from the
          ! center zone is just the original zone minus
          ! the fractions of the left and right zones.

          f( 0) = ONE - f(-1) - f(1)                       

          if (remap_order .eq. 0) then

             do n = 1, U % nvar
                U_av(-1:1,n) = U_temp % data(i-1:i+1,n)
             enddo

          else if (remap_order .eq. 1) then

             do n = 1, U % nvar
             
                ! Our strategy is to construct a linear reconstruction of the
                ! data in each zone, and then integrate over each portion of the
                ! Lagrangian zones that overlap with the original Eulerian zone.

                ul = U_temp % data(i-3:i+1,n)
                uc = U_temp % data(i-2:i+2,n)
                ur = U_temp % data(i-1:i+3,n)

                dxl = (ONE + vf % data(i-2:i+2,1) * dtdx - vf % data(i-3:i+1,1) * dtdx) * dx
                dxc = (ONE + vf % data(i-1:i+3,1) * dtdx - vf % data(i-2:i+2,1) * dtdx) * dx
                dxr = (ONE + vf % data(i  :i+4,1) * dtdx - vf % data(i-1:i+3,1) * dtdx) * dx
                
                ! Construct the limited slopes

                dul = slope(ul, dxl)
                duc = slope(uc, dxc)
                dur = slope(ur, dxr)

                ! Now that we have the slopes, integrate the portion of each zone that overlaps.

                ! First, the contribution from the zone on the left. If it has moved to the right, then
                ! the integral of the state in this region is equal to (v_face * dt)
                ! multiplied by the average value of the reconstructed profile between
                ! x_{i-1/2} and x_{i-1/2} + (v_face * dt). Since the slope is linear,
                ! this is simply equal to the value at x_{i-1/2} + (v_face * dt) / 2.

                U_av(-1, n) = ul(0) + HALF * (dxl(0) - f(-1) * dx) * dul

                ! Contribution from the zone in the center.

                U_av( 0, n) = uc(0) + HALF * (f(-1) * dx - f( 1) * dx) * duc
                
                ! Now, the contribution from the zone on the right.

                U_av( 1, n) = ur(0) - HALF * (dxr(0) - f( 1) * dx) * dur

             enddo

          else

             print *, "Cannot handle remap_order > 1 at present."
             stop

          endif

          ! Now sum up all the contributions.

          do n = 1, U % nvar
             U % data(i,n) = dot_product(f(:), U_av(:,n))
          enddo

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
