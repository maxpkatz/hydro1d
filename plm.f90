module interface_states_plm_module

  use datatypes_module
  use grid_module
  use variables_module
  use params_module
  use eos_module
  use eigen_module
  use flatten_module

  implicit none

  private

  public :: make_interface_states_plm

contains
  
  subroutine make_interface_states_plm(U, U_l, U_r, vf, dt)

    type(gridvar_t),     intent(in   ) :: U
    type(gridedgevar_t), intent(inout) :: U_l, U_r, vf
    real (kind=dp_t),    intent(in   ) :: dt

    type(gridvar_t) :: Q
    type(gridedgevar_t) :: Q_l, Q_r
    type(gridvar_t) :: ldelta_m, ldelta_p, dQ_vl

    type(gridvar_t) :: Q_temp
    
    real (kind=dp_t) :: rvec_m(nwaves,nprim), lvec_m(nwaves,nprim), eval_m(nwaves)
    real (kind=dp_t) :: rvec_p(nwaves,nprim), lvec_p(nwaves,nprim), eval_p(nwaves)    
    real (kind=dp_t) :: dQ_m(nprim), dQ_p(nprim)

    real (kind=dp_t) :: r, ux, ux_m, ux_p, p, cs
    real (kind=dp_t) :: ldr_m, ldu_m, ldp_m
    real (kind=dp_t) :: ldr_p, ldu_p, ldp_p    
    real (kind=dp_t) :: r_xm, r_xp, u_xm, u_xp, p_xm, p_xp
    real (kind=dp_t) :: beta_xm(nwaves), beta_xp(nwaves)
    real (kind=dp_t) :: sum_m, sum_p
    real (kind=dp_t) :: Q_m, Q_c, Q_p
    real (kind=dp_t) :: dQ_l, dQ_c, dQ_r, dQ_lim
    
    type(gridvar_t) :: xi_m, xi_p

    real (kind=dp_t) :: dtdx

    real (kind=dp_t) :: e

    real (kind=dp_t) :: reference_fac_p, reference_fac_m

    integer :: i, m, n

    ! piecewise linear slopes 
    !
    ! This is a 1-d version of the piecewise linear Godunov method
    ! detailed in Colella (1990).  See also Colella & Glaz and
    ! Saltzman (1994).
    !
    ! We wish to solve
    !
    !   U_t + [F(U)]_x = H
    !
    ! we want U_{i+1/2}^{n+1/2} -- the interface values that are input
    ! to the Riemann problem through the faces for each zone.
    !
    ! First we convert from the conserved variables, U = (rho, rho u, rho E)
    ! to the primitive variables, Q = (rho, u, p).  
    !
    ! The system of equations in primitive form appear as:
    !
    !   Q_t + A(Q) Q_x = H'
    !
    ! Then we taylor expand the primitive variable from the
    ! cell-center to the interface at the half-time:
    !
    !  n+1/2        n           dq           dq  
    ! q          = q  + 0.5 dx  --  + 0.5 dt --  
    !  i+1/2,L      i           dx           dt   
    ! 
    !   
    !               n           dq               dq 
    !            = q  + 0.5 dx  --  - 0.5 dt ( A -- - H' ) 
    !               i           dx               dx  
    !      
    !       
    !               n                dt     dq 
    !            = q  + 0.5 dx ( 1 - -- A ) -- + 0.5 dt H'
    !               i                dx     dx      
    !   
    !   
    !               n              dt     _ 
    !            = q  + 0.5  ( 1 - -- A ) Dq + 0.5 dt H'   
    !               i              dx     
    !     
    !                   +---------+---------+ +---+---+    
    !                             |               |   
    !               
    !                 this is the monotonized   source term    
    !                 central difference term     
  

    
    ! sanity check
    if (U%grid%ng < 4) then
       print *, "ERROR: ng < 4 in plm states"
       stop
    endif


    !-------------------------------------------------------------------------
    ! convert to primitve variables
    !-------------------------------------------------------------------------
    call build(Q, U%grid, nprim)

    do i = U%grid%lo-U%grid%ng, U%grid%hi+U%grid%ng
       ! density
       Q%data(i,iqdens) = U%data(i,iudens)

       ! velocity
       Q%data(i,iqxvel) = U%data(i,iumomx)/U%data(i,iudens)

       ! pressure
       e = (U%data(i,iuener) - &
            HALF*U%data(i,iumomx)**2/U%data(i,iudens))/U%data(i,iudens)
       call eos(eos_input_e, Q%data(i,iqpres), e, Q%data(i,iqdens))
    enddo


    !-------------------------------------------------------------------------
    ! compute the flattening coefficients
    !-------------------------------------------------------------------------

    if (use_flattening) then
       
       call build(xi_m, U%grid, 1)
       call build(xi_p, U%grid, 1)

       Q_temp = Q
    
       do i = U%grid%lo - U%grid%ng+1, U%grid%hi+U%grid%ng
          Q_temp%data(i,2) = Q_temp%data(i,2) - vf % data(i,1)
       enddo
    
       call flatten(Q_temp, xi_m)

       Q_temp = Q
    
       do i = U%grid%lo - U%grid%ng+1, U%grid%hi+U%grid%ng
          Q_temp%data(i,2) = Q_temp%data(i,2) - vf % data(i+1,1)
       enddo
       
       call flatten(Q_temp, xi_p)


    endif
       
    !-------------------------------------------------------------------------
    ! compute the monotonized central differences
    !-------------------------------------------------------------------------

    ! 4th order MC limiting.  See Colella (1985) Eq. 2.5 and 2.6,
    ! Colella (1990) page 191 (with the delta a terms all equal) or
    ! Saltzman 1994, page 156
    
    call build(dQ_vl, U%grid, 1)
    call build(ldelta_m, U%grid, nprim)
    call build(ldelta_p, U%grid, nprim)

    ! First do the left-edge slopes

    do n = 1, nprim

       ! first do the normal MC limiting 
       do i = Q%grid%lo-3, Q%grid%hi+3

          Q_m = Q % data(i-1,n)
          Q_c = Q % data(i  ,n)
          Q_p = Q % data(i+1,n)

          if (n == iqxvel) then
             Q_m = Q_m - vf % data(i,1)
             Q_c = Q_c - vf % data(i,1)
             Q_p = Q_p - vf % data(i,1)
          endif

          dQ_l = Q_c - Q_m
          dQ_c = HALF * (Q_p - Q_m)
          dQ_r = Q_p - Q_c

          dQ_lim = 2.0_dp_t * min(abs(dQ_r), abs(dQ_l))
          
          if (dQ_r * dQ_l > ZERO) then
             dQ_vl % data(i,1) = min(abs(dQ_c), abs(dQ_lim)) * sign(ONE,dQ_c)
          else
             dQ_vl % data(i,1) = ZERO
          endif

          ldelta_m % data(i,n) = dQ_vl % data(i,1)
          
       enddo


       ! now do the fourth order part

       if (use_higher_order_limiter) then

          do i = Q%grid%lo-2, Q%grid%hi+2

             Q_m = Q%data(i-1,n)
             Q_c = Q%data(i  ,n)
             Q_p = Q%data(i+1,n)

             if (n == iqxvel) then
                Q_m = Q_m - vf % data(i,1)
                Q_c = Q_c - vf % data(i,1)
                Q_p = Q_p - vf % data(i,1)
             endif

             dQ_l = Q_c - Q_m
             dQ_c = HALF * (Q_p - Q_m)
             dQ_r = Q_p - Q_c          

             if (dQ_l * dQ_r > ZERO) then
                ldelta_m % data(i,n) = &
                     min((2.0_dp_t/3.0_dp_t)*abs(2.0_dp_t * dQ_c - &
                     0.25_dp_t*(dQ_vl % data(i+1,1) + dQ_vl % data(i-1,1))), &
                     min(2.0*abs(dQ_r), 2.0*abs(dQ_l))) * &
                     sign(ONE, Q_p - Q_m)
             else
                ldelta_m % data(i,n) = ZERO
             endif

          enddo

       endif
          
    enddo


    ! apply flattening to the slopes
    if (use_flattening) then
       do n = 1, nprim
          do i = Q%grid%lo-2, Q%grid%hi+2
             ldelta_m%data(i,n) = xi_m%data(i,1)*ldelta_m%data(i,n)
          enddo
       enddo
    endif



    ! Now do the right-edge slopes
    
    do n = 1, nprim

       ! first do the normal MC limiting 
       do i = Q%grid%lo-3, Q%grid%hi+3

          Q_m = Q % data(i-1,n)
          Q_c = Q % data(i  ,n)
          Q_p = Q % data(i+1,n)

          if (n == iqxvel) then
             Q_m = Q_m - vf % data(i+1,1)
             Q_c = Q_c - vf % data(i+1,1)
             Q_p = Q_p - vf % data(i+1,1)
          endif

          dQ_l = Q_c - Q_m
          dQ_c = HALF * (Q_p - Q_m)
          dQ_r = Q_p - Q_c                    

          dQ_lim = 2.0_dp_t * min(abs(dQ_r), abs(dQ_l))          
          
          if (dQ_l * dQ_r > ZERO) then
             dQ_vl % data(i,1) = min(abs(dQ_c), abs(dQ_lim)) * sign(ONE, dQ_c)
          else
             dQ_vl % data(i,1) = ZERO
          endif

          ldelta_p % data(i,n) = dQ_vl % data(i,1)
          
       enddo
    

       ! now do the fourth order part

       if (use_higher_order_limiter) then
          
          do i = Q%grid%lo-2, Q%grid%hi+2

             Q_m = Q%data(i-1,n)
             Q_c = Q%data(i  ,n)
             Q_p = Q%data(i+1,n)

             if (n == iqxvel) then
                Q_m = Q_m - vf % data(i+1,1)
                Q_c = Q_c - vf % data(i+1,1)
                Q_p = Q_p - vf % data(i+1,1)
             endif

             dQ_l = Q_c - Q_m
             dQ_c = HALF * (Q_p - Q_m)
             dQ_r = Q_p - Q_c                    

             if (dQ_l * dQ_r > ZERO) then
                ldelta_p % data(i,n) = &
                     min((2.0_dp_t/3.0_dp_t)*abs(2.0_dp_t * dQ_c - &
                     0.25_dp_t*(dQ_vl % data(i+1,1) + dQ_vl % data(i-1,1))), &
                     min(2.0*abs(dQ_r), 2.0*abs(dQ_l))) * &
                     sign(ONE, dQ_c)
             else
                ldelta_p % data(i,n) = ZERO
             endif

          enddo

       endif
          
    enddo


    ! apply flattening to the slopes
    if (use_flattening) then
       do n = 1, nprim
          do i = Q%grid%lo-2, Q%grid%hi+2
             ldelta_p%data(i,n) = xi_p%data(i,1)*ldelta_p%data(i,n)
          enddo
       enddo
    endif
       

    
    call destroy(dQ_vl)
    if (use_flattening) then
       call destroy(xi_m)
       call destroy(xi_p)
    endif

    !-------------------------------------------------------------------------
    ! compute left and right primitive variable states
    !-------------------------------------------------------------------------
    call build(Q_l, U%grid, nprim)
    call build(Q_r, U%grid, nprim)

    dtdx = dt/U%grid%dx

    ! The basic idea here is that we do a characteristic
    ! decomposition.  The jump in primitive variables (Q) can be
    ! transformed to a jump in characteristic variables using the left
    ! and right eigenvectors.  Then each wave tells us how much of
    ! each characteristic quantity reaches the interface over dt/2.
    ! We only add the quantity if it moves toward the interface.
    !
    ! Following Colella & Glaz, and Colella (1990), we pick a
    ! reference state to minimize the size of the jump that the
    ! projection operates on---this is because our equations are not
    ! linear.
    !
    ! So 
    !
    !  n+1/2                n                     dt     -
    ! q         - q     =  q   - q    + 0.5 ( 1 - -- A ) Dq
    !  i+1/2,L     ref      i     ref             dx
    !
    !
    ! The reference state is chosen as (Colella Eq. 2.11; Colella &
    ! Glaz, p. 278):
    !
    !         ~                   dt           +       -
    ! q     = q  = q   + 0.5 [1 - -- max(lambda , 0) ] Dq
    !  ref     L    i             dx
    !
    ! We project the RHS using the left and right eigenvectors, and only
    ! consider those waves moving toward the interface.  This gives:
    !
    !  n+1/2      ~        dt                        +           -
    ! q         = q  + 0.5 --  sum { l  . [max(lambda , 0) - A ] Dq r  }
    !  i+1/2,L     L       dx   i     i                              i
    !
    ! since l A = lambda l, we have:
    !
    !  n+1/2      ~        dt                   +                     -
    ! q         = q  + 0.5 --  sum { [max(lambda , 0) - lambda ] (l . Dq) r  }
    !  i+1/2,L     L       dx   i                             i    i       i
    !
    ! See Miller & Colella (2002) for more details.  This expression is 
    ! found in Colella (1990) at the bottom of p. 191.

    do i = U%grid%lo-1, U%grid%hi+1

       r    = Q%data(i,iqdens)
       ux   = Q%data(i,iqxvel)
       ux_m = Q%data(i,iqxvel) - vf % data(i  ,1)
       ux_p = Q%data(i,iqxvel) - vf % data(i+1,1)
       p    = Q%data(i,iqpres)

       ldr_m = ldelta_m%data(i,iqdens)
       ldu_m = ldelta_m%data(i,iqxvel)
       ldp_m = ldelta_m%data(i,iqpres)

       dQ_m(:) = [ ldr_m, ldu_m, ldp_m ]

       ldr_p = ldelta_p%data(i,iqdens)
       ldu_p = ldelta_p%data(i,iqxvel)
       ldp_p = ldelta_p%data(i,iqpres)

       dQ_p(:) = [ ldr_p, ldu_p, ldp_p ]

       ! Subtract off the 
       
       ! compute the sound speed
       cs = sqrt(gamma*p/r)


       ! get the eigenvalues and eigenvectors
       call eigen(r, ux, p, cs, lvec_p, rvec_p, eval_p)
       call eigen(r, ux, p, cs, lvec_m, rvec_m, eval_m)

       ! Define the reference states (here xp is the right interface
       ! for the current zone and xm is the left interface for the
       ! current zone)
       !                           ~
       ! These expressions are the V_{L,R} in Colella (1990) at the 
       ! bottom of page 191.  They are also in Saltzman (1994) as
       ! V_ref on page 161.

       if (reference_state == 1) then
          reference_fac_p = HALF * (ONE - dtdx * max(eval_p(3), ZERO))
          reference_fac_m = HALF * (ONE + dtdx * min(eval_m(1), ZERO))
       else
          reference_fac_p = ZERO
          reference_fac_m = ZERO
       endif
       
       r_xp = r    + reference_fac_p * ldr_p
       r_xm = r    - reference_fac_m * ldr_m

       u_xp = ux_p + reference_fac_p * ldu_p
       u_xm = ux_m - reference_fac_m * ldu_m

       p_xp = p    + reference_fac_p * ldp_p
       p_xm = p    - reference_fac_m * ldp_m

       !                                                   ^
       ! Now compute the interface states.   These are the V expressions
       ! in Colella (1990) page 191, and the interface state expressions
       ! -V and +V in Saltzman (1994) on pages 161.
       
       ! first compute beta_xm and beta_xp -- these are the
       ! coefficients to the right eigenvectors in the eigenvector
       ! expansion (see Colella 1990, page 191)
       do m = 1, nwaves

          ! dot product of the current left eigenvector with the
          ! primitive variable jump
          sum_m = dot_product(lvec_m(m,:),dQ_m(:))
          sum_p = dot_product(lvec_p(m,:),dQ_p(:))
          
          if (use_tracing) then

             if (reference_state == 1) then
                
                ! here the sign() function makes sure we only add the right-moving
                ! waves            
                beta_xp(m) = 0.25_dp_t*dtdx*(eval_p(3) - eval_p(m))* &
                     (sign(ONE, eval_p(m)) + ONE)*sum_p

                ! here the sign() function makes sure we only add the left-moving
                ! waves
                beta_xm(m) = 0.25_dp_t*dtdx*(eval_m(1) - eval_m(m))* &
                     (ONE - sign(ONE, eval_m(m)))*sum_m

             else

                beta_xp(m) =  0.25_dp_t*(ONE - dtdx*eval_p(m))* &
                     (sign(ONE, eval_p(m)) + ONE)*sum_p

                beta_xm(m) = -0.25_dp_t*(ONE - dtdx*eval_m(m))* &
                     (ONE - sign(ONE, eval_m(m)))*sum_m                

             endif
                
          else

             if (reference_state == 1) then
                beta_xp(m) = HALF * dtdx * (eval_p(3) - eval_p(m)) * sum_p
                beta_xm(m) = HALF * dtdx * (eval_m(1) - eval_m(m)) * sum_m
             else
                beta_xp(m) =  HALF * (ONE - dtdx * eval_p(m)) * sum_p
                beta_xm(m) = -HALF * (ONE + dtdx * eval_m(m)) * sum_m
             endif

          endif
             
       enddo


       ! finally, sum up all the jumps 
       
       ! density
       Q_l%data(i+1,iqdens) = r_xp + sum(beta_xp(:) * rvec_p(:,iqdens))
       Q_r%data(i  ,iqdens) = r_xm + sum(beta_xm(:) * rvec_m(:,iqdens))

       ! velocity
       Q_l%data(i+1,iqxvel) = u_xp + sum(beta_xp(:) * rvec_p(:,iqxvel))
       Q_r%data(i  ,iqxvel) = u_xm + sum(beta_xm(:) * rvec_m(:,iqxvel))

       ! pressure
       Q_l%data(i+1,iqpres) = p_xp + sum(beta_xp(:) * rvec_p(:,iqpres))
       Q_r%data(i  ,iqpres) = p_xm + sum(beta_xm(:) * rvec_m(:,iqpres))

       
    enddo

    ! clean-up
    call destroy(ldelta_p)
    call destroy(ldelta_m)


    !-------------------------------------------------------------------------
    ! apply the source terms
    !-------------------------------------------------------------------------
    do i = U%grid%lo, U%grid%hi+1
       Q_l%data(i,iqxvel) = Q_l%data(i,iqxvel) + HALF*dt*grav
       Q_r%data(i,iqxvel) = Q_r%data(i,iqxvel) + HALF*dt*grav
    enddo

    ! special fixes at the boundary -- gravity must be reflected
    ! (correct the above too)
    if (U%grid%xlboundary == "reflect") then
       Q_l%data(U%grid%lo,iqxvel) = &
            Q_l%data(U%grid%lo,iqxvel) - dt*grav
    endif

    if (U%grid%xrboundary == "reflect") then
       Q_r%data(U%grid%hi+1,iqxvel) = &
            Q_r%data(U%grid%hi+1,iqxvel) - dt*grav
    endif


    !-------------------------------------------------------------------------
    ! transform the states into conserved variables
    !-------------------------------------------------------------------------
    do i = U%grid%lo, U%grid%hi+1

       ! density
       U_l%data(i,iudens) = Q_l%data(i,iqdens) 
       U_r%data(i,iudens) = Q_r%data(i,iqdens) 

       ! momentum
       U_l%data(i,iumomx) = Q_l%data(i,iqdens)*Q_l%data(i,iqxvel)
       U_r%data(i,iumomx) = Q_r%data(i,iqdens)*Q_r%data(i,iqxvel)

       ! total energy
       call eos(eos_input_p, Q_l%data(i,iqpres), e, Q_l%data(i,iqdens))
       U_l%data(i,iuener) = Q_l%data(i,iqdens)*e + &
            HALF*Q_l%data(i,iqdens)*Q_l%data(i,iqxvel)**2

       call eos(eos_input_p, Q_r%data(i,iqpres), e, Q_r%data(i,iqdens))
       U_r%data(i,iuener) = Q_r%data(i,iqdens)*e + &
            HALF*Q_r%data(i,iqdens)*Q_r%data(i,iqxvel)**2

    enddo
    

    ! clean-up
    call destroy(Q)
    call destroy(Q_l)
    call destroy(Q_r)


  end subroutine make_interface_states_plm

end module interface_states_plm_module
