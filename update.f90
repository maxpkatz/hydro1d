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

    real (kind=dp_t) :: dtdx
    real (kind=dp_t) :: U_av(-1:1, U%nvar), f(-1:1)
    real (kind=dp_t) :: ull, ul, uc, ur, urr
    real (kind=dp_t) :: dul, duc, dur
    real (kind=dp_t) :: dul_l, dul_c, dul_r
    real (kind=dp_t) :: duc_l, duc_c, duc_r
    real (kind=dp_t) :: dur_l, dur_c, dur_r
    real (kind=dp_t) :: dx, dxll, dxl, dxc, dxr, dxrr
    
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

             dxll = (ONE + vf % data(i-1,1) * dtdx - vf % data(i-2,1) * dtdx) * dx
             dxl  = (ONE + vf % data(i  ,1) * dtdx - vf % data(i-1,1) * dtdx) * dx
             dxc  = (ONE + vf % data(i+1,1) * dtdx - vf % data(i  ,1) * dtdx) * dx
             dxr  = (ONE + vf % data(i+2,1) * dtdx - vf % data(i+1,1) * dtdx) * dx
             dxrr = (ONE + vf % data(i+3,1) * dtdx - vf % data(i+2,1) * dtdx) * dx             
             
             do n = 1, U % nvar
             
                ! Our strategy is to construct a linear reconstruction of the
                ! data in each zone, and then integrate over each portion of the
                ! Lagrangian zones that overlap with the original Eulerian zone.

                ull = U_temp % data(i-2,n)
                ul  = U_temp % data(i-1,n)
                uc  = U_temp % data(i  ,n)
                ur  = U_temp % data(i+1,n)
                urr = U_temp % data(i+2,n)

                ! For example, in zone i, the slope is u'(x) = (u_{i+1} - u_{i-1}) / (x_{i+1} - x_{i-1}),
                ! where the zone indices are evaluated in the Lagrangian sense.

                dul_l = (ul  - ull) / (HALF * dxl + HALF * dxll)
                dul_c = (uc  - ull) / (dxl + HALF * (dxll + dxc))
                dul_r = (uc  - ul ) / (HALF * dxl + HALF * dxc )

                dul = TWO * min(abs(dul_r), abs(dul_l))
                
                if (dul_r * dul_l > ZERO) then
                   dul = min(abs(dul_c), abs(dul)) * sign(ONE, dul_c)
                else
                   dul = ZERO
                endif
                
                duc_l = (uc  - ul ) / (HALF * dxc + HALF * dxl )
                duc_c = (ur  - ul ) / (dxc + HALF * (dxl  + dxr))
                duc_r = (ur  - uc ) / (HALF * dxr + HALF * dxc )

                duc = TWO * min(abs(duc_r), abs(duc_l))
                
                if (duc_r * duc_l > ZERO) then
                   duc = min(abs(duc_c), abs(duc)) * sign(ONE, duc_c)
                else
                   duc = ZERO
                endif
                
                dur_l = (ur  - uc ) / (HALF * dxr + HALF * dxc )
                dur_c = (urr - uc ) / (dxr + HALF * (dxc  + dxrr))
                dur_r = (urr - ur ) / (HALF * dxrr + HALF * dxr)

                dur = TWO * min(abs(dur_r), abs(dur_l))
                
                if (dur_r * dur_l > ZERO) then
                   dur = min(abs(dur_c), abs(dur)) * sign(ONE, dur_c)
                else
                   dur = ZERO
                endif

                ! Now that we have the slopes, integrate the portion of each zone that overlaps.

                ! First, the contribution from the zone on the left. If it has moved to the right, then
                ! the integral of the state in this region is equal to (v_face * dt)
                ! multiplied by the average value of the reconstructed profile between
                ! x_{i-1/2} and x_{i-1/2} + (v_face * dt). Since the slope is linear,
                ! this is simply equal to the value at x_{i-1/2} + (v_face * dt) / 2.

                U_av(-1, n) = ul + HALF * (dxl - f(-1) * dx) * dul

                ! Contribution from the zone in the center.

                U_av( 0, n) = uc + HALF * (f(-1) * dx - f( 1) * dx) * duc
                
                ! Now, the contribution from the zone on the right.

                U_av( 1, n) = ur - HALF * (dxr - f( 1) * dx) * dur

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
