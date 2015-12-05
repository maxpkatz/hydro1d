module interface_states_ppm_module

  use datatypes_module
  use grid_module
  use variables_module
  use params_module
  use eos_module
  use eigen_module
  use flatten_module

  implicit none

  private

  public :: make_interface_states_ppm

contains
  
  subroutine make_interface_states_ppm(U, U_l, U_r, vf, dt)

    type(gridvar_t),     intent(in   ) :: U
    type(gridedgevar_t), intent(inout) :: U_l, U_r, vf
    real (kind=dp_t),    intent(in   ) :: dt

    type(gridvar_t) :: Q
    type(gridedgevar_t) :: Q_l, Q_r

    real (kind=dp_t) :: rvec_m(nwaves,nprim), lvec_m(nwaves,nprim), eval_m(nwaves)
    real (kind=dp_t) :: rvec_p(nwaves,nprim), lvec_p(nwaves,nprim), eval_p(nwaves)    

    real (kind=dp_t) :: r, ux_m, ux_p, p, cs
    real (kind=dp_t) :: beta_xm(nwaves), beta_xp(nwaves)

    real (kind=dp_t) :: dq0, dqp
    real (kind=dp_t) :: Iminus(nwaves,nprim), Iplus(nwaves,nprim)
    real (kind=dp_t) :: Qref_xm(nprim), Qref_xp(nprim)

    type(gridedgevar_t) :: Qminus_m, Qplus_m, Qminus_p, Qplus_p
    type(gridvar_t) :: Q6_m, Q6_p

    real (kind=dp_t) :: sigma_m, sigma_p
    real (kind=dp_t) :: Q_mm, Q_m, Q_c, Q_p, Q_pp

    type(gridvar_t) :: Q_temp
    type(gridvar_t) :: xi_m, xi_p

    real (kind=dp_t) :: dtdx

    real (kind=dp_t) :: e

    integer :: i, m, n

    ! piecewise parabolic slopes 
    !
    ! This is a 1-d version of the piecewise parabolic method detailed
    ! Colella & Woodward (1984).  We follow the description of Almgren
    ! et al. 2010 (the CASTRO paper) and Miller and Colella (2002).  
    !
    ! Also note that we do not implement the contact steepening here.
    
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
    ! interpolate the cell-centered data to the edges
    !-------------------------------------------------------------------------

    ! For each cell, we will find the Qminus and Qplus states -- these
    ! are the - and + edges of the cell.
    call build(Qminus_m, U%grid, nprim)
    call build(Qplus_m,  U%grid, nprim)
    
    do n = 1, nprim
       do i = U%grid%lo-2, U%grid%hi+2

          Q_mm = Q % data(i-2,n)
          Q_m  = Q % data(i-1,n)
          Q_c  = Q % data(i  ,n)
          Q_p  = Q % data(i+1,n)
          Q_pp = Q % data(i+2,n)

          if (n == iqxvel) then
             Q_mm = Q_mm - vf % data(i,1)
             Q_m  = Q_m  - vf % data(i,1)
             Q_c  = Q_c  - vf % data(i,1)
             Q_p  = Q_p  - vf % data(i,1)
             Q_pp = Q_pp - vf % data(i,1)
          endif          
          
          ! dq (C&W Eq. 1.7)
          dq0 = HALF*(Q_p - Q_m)
          dqp = HALF*(Q_pp - Q_c)

          ! limiting (C&W Eq. 1.8)
          if ( (Q_p - Q_c) * (Q_c - Q_m) > ZERO) then
             dq0 = sign(ONE,dq0) * min( abs(dq0), TWO*abs(Q_c - Q_m), TWO*abs(Q_p - Q_c) )
          else
             dq0 = ZERO
          endif

          if ( (Q_pp - Q_p) * (Q_p - Q_c) > ZERO) then
             dqp = sign(ONE,dqp) * min( abs(dqp), TWO*abs(Q_p - Q_c), TWO*abs(Q_pp - Q_p) )
          else
             dqp = ZERO
          endif

          ! cubic (C&W Eq. 1.6)
          Qplus_m%data(i,n) = HALF*(Q_c + Q_p) - (ONE/SIX)*(dqp - dq0)

          Qminus_m%data(i+1,n) = Qplus_m%data(i,n)

          ! make sure that we didn't over or undersoot -- this may not
          ! be needed, but is discussed in Colella & Sekora (2008)
          Qplus_m%data(i,n) = max(Qplus_m%data(i,n), min(Q_c,Q_p))
          Qplus_m%data(i,n) = min(Qplus_m%data(i,n), max(Q_c,Q_p))

          Qminus_m%data(i+1,n) = max(Qminus_m%data(i+1,n), min(Q_c,Q_p))
          Qminus_m%data(i+1,n) = min(Qminus_m%data(i+1,n), max(Q_c,Q_p))

       enddo
    enddo


    !-------------------------------------------------------------------------
    ! construct the parameters for the parabolic reconstruction polynomials
    !-------------------------------------------------------------------------

    ! Our parabolic profile has the form:
    !
    !  q(xi) = qminus + xi*(qplus - qminus + q6 * (1-xi) )
    !
    ! with xi = (x - xl)/dx, where xl is the interface of the left
    ! edge of the cell.  qminus and qplus are the values of the 
    ! parabola on the left and right edges of the current cell.

    ! Limit the left and right values of the parabolic interpolant
    ! (C&W Eq. 1.10).  Here the loop is over cells, and considers the
    ! values on either side of the center of the cell (Qminus and
    ! Qplus).
    do n = 1, nprim
       do i = U%grid%lo-1, U%grid%hi+1

          Q_c = Q % data(i,n)

          if (n == iqxvel) then
             Q_c = Q_c - vf % data(i,1)
          endif                    
          
          if ( (Qplus_m%data(i,n) - Q_c) * (Q_c - Qminus_m%data(i,n)) <= ZERO) then
             Qminus_m%data(i,n) = Q_c
             Qplus_m%data(i,n) = Q_c

          else if ( (Qplus_m%data(i,n) - Qminus_m%data(i,n)) * &
                    (Q_c - &
                      HALF*(Qminus_m%data(i,n) + Qplus_m%data(i,n))) > &
                   (Qplus_m%data(i,n) - Qminus_m%data(i,n))**2/SIX ) then

          ! alternate test from Colella & Sekora (2008)
          !else if (abs(Qminus%data(i,n) - Q%data(i,n)) >= &
          !     2.0*abs(Qplus%data(i,n) - Q%data(i,n))) then
             Qminus_m%data(i,n) = THREE*Q_c - TWO*Qplus_m%data(i,n)

          else if (-(Qplus_m%data(i,n) - Qminus_m%data(i,n))**2/SIX > &
                    (Qplus_m%data(i,n) - Qminus_m%data(i,n)) * &
                    (Q_c - &
                           HALF*(Qminus_m%data(i,n) + Qplus_m%data(i,n))) ) then

          !else if (abs(Qplus%data(i,n) - Q%data(i,n)) >= &
          !     2.0*abs(Qminus%data(i,n) - Q%data(i,n))) then
             Qplus_m%data(i,n) = THREE*Q_c - TWO*Qminus_m%data(i,n)

          endif

       enddo
    enddo

    ! define Q6
    call build(Q6_m,  U%grid, nprim)

    do n = 1, nprim
       do i = U%grid%lo-1, U%grid%hi+1

          Q_c = Q % data(i,n)

          if (n == iqxvel) then
             Q_c = Q_c - vf % data(i,1)
          endif
          
          Q6_m%data(i,n) = SIX*Q_c - THREE*(Qminus_m%data(i,n) + Qplus_m%data(i,n))

       enddo
    enddo


    !-------------------------------------------------------------------------
    ! interpolate the cell-centered data to the edges
    !-------------------------------------------------------------------------

    ! For each cell, we will find the Qminus and Qplus states -- these
    ! are the - and + edges of the cell.
    call build(Qminus_p, U%grid, nprim)
    call build(Qplus_p,  U%grid, nprim)
    
    do n = 1, nprim
       do i = U%grid%lo-2, U%grid%hi+2

          Q_mm = Q % data(i-2,n)
          Q_m  = Q % data(i-1,n)
          Q_c  = Q % data(i  ,n)
          Q_p  = Q % data(i+1,n)
          Q_pp = Q % data(i+2,n)

          if (n == iqxvel) then
             Q_mm = Q_mm - vf % data(i+1,1)
             Q_m  = Q_m  - vf % data(i+1,1)
             Q_c  = Q_c  - vf % data(i+1,1)
             Q_p  = Q_p  - vf % data(i+1,1)
             Q_pp = Q_pp - vf % data(i+1,1)
          endif
          
          ! dq (C&W Eq. 1.7)
          dq0 = HALF*(Q_p - Q_m)
          dqp = HALF*(Q_pp - Q_c)

          ! limiting (C&W Eq. 1.8)
          if ( (Q_p - Q_c) * (Q_c - Q_m) > ZERO) then
             dq0 = sign(ONE,dq0) * min( abs(dq0), TWO*abs(Q_c - Q_m), TWO*abs(Q_p - Q_c) )
          else
             dq0 = ZERO
          endif

          if ( (Q_pp - Q_p) * (Q_p - Q_c) > ZERO) then
             dqp = sign(ONE,dqp) * min( abs(dqp), TWO*abs(Q_p - Q_c), TWO*abs(Q_pp - Q_p) )
          else
             dqp = ZERO
          endif

          ! cubic (C&W Eq. 1.6)
          Qplus_p%data(i,n) = HALF*(Q_c + Q_p) - (ONE/SIX)*(dqp - dq0)

          Qminus_p%data(i+1,n) = Qplus_p%data(i,n)

          ! make sure that we didn't over or undersoot -- this may not
          ! be needed, but is discussed in Colella & Sekora (2008)
          Qplus_p%data(i,n) = max(Qplus_p%data(i,n), min(Q_c,Q_p))
          Qplus_p%data(i,n) = min(Qplus_p%data(i,n), max(Q_c,Q_p))

          Qminus_p%data(i+1,n) = max(Qminus_p%data(i+1,n), min(Q_c,Q_p))
          Qminus_p%data(i+1,n) = min(Qminus_p%data(i+1,n), max(Q_c,Q_p))

       enddo
    enddo


    !-------------------------------------------------------------------------
    ! construct the parameters for the parabolic reconstruction polynomials
    !-------------------------------------------------------------------------

    ! Our parabolic profile has the form:
    !
    !  q(xi) = qminus + xi*(qplus - qminus + q6 * (1-xi) )
    !
    ! with xi = (x - xl)/dx, where xl is the interface of the left
    ! edge of the cell.  qminus and qplus are the values of the 
    ! parabola on the left and right edges of the current cell.

    ! Limit the left and right values of the parabolic interpolant
    ! (C&W Eq. 1.10).  Here the loop is over cells, and considers the
    ! values on either side of the center of the cell (Qminus_p and
    ! Qplus).
    do n = 1, nprim
       do i = U%grid%lo-1, U%grid%hi+1

          Q_c = Q % data(i,n)

          if (n == iqxvel) then
             Q_c = Q_c - vf % data(i+1,1)
          endif
          
          if ( (Qplus_p%data(i,n) - Q_c) * &
               (Q_c - Qminus_p%data(i,n)) <= ZERO) then
             Qminus_p%data(i,n) = Q_c
             Qplus_p%data(i,n) = Q_c

          else if ( (Qplus_p%data(i,n) - Qminus_p%data(i,n)) * &
                    (Q_c - &
                      HALF*(Qminus_p%data(i,n) + Qplus_p%data(i,n))) > &
                   (Qplus_p%data(i,n) - Qminus_p%data(i,n))**2/SIX ) then

          ! alternate test from Colella & Sekora (2008)
          !else if (abs(Qminus%data(i,n) - Q_c) >= &
          !     2.0*abs(Qplus%data(i,n) - Q_c)) then
             Qminus_p%data(i,n) = THREE*Q_c - TWO*Qplus_p%data(i,n)

          else if (-(Qplus_p%data(i,n) - Qminus_p%data(i,n))**2/SIX > &
                    (Qplus_p%data(i,n) - Qminus_p%data(i,n)) * &
                    (Q_c - &
                           HALF*(Qminus_p%data(i,n) + Qplus_p%data(i,n))) ) then

          !else if (abs(Qplus%data(i,n) - Q_c) >= &
          !     2.0*abs(Qminus%data(i,n) - Q_c)) then
             Qplus_p%data(i,n) = THREE*Q_c - TWO*Qminus_p%data(i,n)

          endif

       enddo
    enddo

    ! define Q6
    call build(Q6_p,  U%grid, nprim)

    do n = 1, nprim
       do i = U%grid%lo-1, U%grid%hi+1

          Q_c = Q % data(i,n)

          if (n == iqxvel) then
             Q_c = Q_c - vf % data(i+1,1)
          endif
          
          Q6_p%data(i,n) = SIX*Q_c - THREE*(Qminus_p%data(i,n) + Qplus_p%data(i,n))

       enddo
    enddo

    
    

    !-------------------------------------------------------------------------
    ! compute left and right primitive variable states
    !-------------------------------------------------------------------------

    call build(Q_l, U%grid, nprim)
    call build(Q_r, U%grid, nprim)

    dtdx = dt/U%grid%dx

    ! The fluxes are going to be defined on the left edge of the
    ! computational zone.
    !
    !          |             |             |             | 
    !          |             |             |             |  
    !         -+------+------+------+------+------+------+--   
    !          |     i-1     |      i      |     i+1     |   
    !                       * *           * 
    !                   q_l,i  q_r,i   q_l,i+1 
    !           
    ! q_l,i+1 are computed using the information in zone i.
    !

    do i = U%grid%lo-1, U%grid%hi+1

       r    = Q%data(i,iqdens)
       ux_m = Q%data(i,iqxvel) - vf % data(i  ,1)
       ux_p = Q%data(i,iqxvel) - vf % data(i+1,1)
       p    = Q%data(i,iqpres)

       ! compute the sound speed
       cs = sqrt(gamma*p/r)


       ! get the eigenvalues and eigenvectors
       call eigen(r, ux_m, p, cs, lvec_m, rvec_m, eval_m)
       call eigen(r, ux_p, p, cs, lvec_p, rvec_p, eval_p)       


       ! integrate the parabola in the cell from the left interface
       ! (Iminus) over the portion of the cell that each eigenvalue
       ! can reach.  Do the same from the right interface in the cell,
       ! defining Iplus.  See Almgren et al. 2010 (Eq. 30) or Colella
       ! & Sekora (2008), or Miller & Colella (2002), Eq. 90.
       do m = 1, nwaves
          sigma_m = abs(eval_m(m))*dtdx
          sigma_p = abs(eval_p(m))*dtdx
          
          do n = 1, nprim

             ! only integrate if the wave is moving toward the interface
             ! (see Miller & Colella, Eg. 90).  This may not be necessary.
             if (eval_p(m) >= ZERO) then
                Iplus(m,n) = Qplus_p%data(i,n) - HALF*sigma_p* &
                     (Qplus_p%data(i,n) - Qminus_p%data(i,n) - &
                     (ONE - (TWO/THREE)*sigma_p)*Q6_p%data(i,n))
             else
                Iplus(m,n) = Q%data(i,n)

                if (n == iqxvel) then
                   Iplus(m,n) = Iplus(m,n) - vf % data(i+1,1)
                endif
             endif

             if (eval_m(m) <= ZERO) then
                Iminus(m,n) = Qminus_m%data(i,n) + HALF*sigma_m* &
                     (Qplus_m%data(i,n) - Qminus_m%data(i,n) + &
                     (ONE - (TWO/THREE)*sigma_m)*Q6_m%data(i,n))
             else
                Iminus(m,n) = Q%data(i,n)

                if (n == iqxvel) then
                   Iminus(m,n) = Iminus(m,n) - vf % data(i,1)
                endif
             endif

          enddo
       enddo


       ! the basic idea here is that we do a characteristic
       ! decomposition.  The jump in primitive variables (Q) can be
       ! transformed to a jump in characteristic variables using the
       ! left and right eigenvectors.  Then each wave tells us how
       ! much of each characteristic quantity reaches the interface
       ! over dt/2.  We only add the quantity if it moves toward the
       ! interface.
       !
       ! See Miller & Colella for a good discussion of the
       ! characteristic form.  The basic form is:
       !
       !
       !  n+1/2      ~               ~
       ! q         = q  -  sum l . ( q   - I  ) r
       !  i+1/2,L     L     i   i     L     +    i
       ! 
       !       ~
       ! Where q is the reference state.

       
       ! define the reference states -- Miller & Colella (2002) argue
       ! picking the fastest wave speed.  We follow the convention from
       ! Colella & Glaz, and Colella (1990) -- this is intended to 
       ! minimize the size of the jump that the projection operates on.
       if (.false.) then
          ! CASTRO method
          Qref_xm(:) = Q%data(i,:)
          Qref_xp(:) = Q%data(i,:)
       else
          ! Miller and Colella method
          if (eval_p(3) >= ZERO) then
             Qref_xp(:) = Iplus(3,:)
          else
             Qref_xp(:) = Q%data(i,:)
          endif

          if (eval_m(1) <= ZERO) then
             Qref_xm(:) = Iminus(1,:)
          else
             Qref_xm(:) = Q%data(i,:)
          endif             
       endif


       ! compute the dot product of each left eigenvector with (qref - I)
       do m = 1, nwaves    ! loop over waves
          beta_xm(m) = dot_product(lvec_m(m,:),Qref_xm(:) - Iminus(m,:))
          beta_xp(m) = dot_product(lvec_p(m,:),Qref_xp(:) - Iplus(m,:) )
       enddo

       ! finally, sum up all the jumps

       ! density
       Q_l%data(i+1,iqdens) = ZERO
       Q_r%data(i  ,iqdens) = ZERO

       do n = 1, nwaves
          if (eval_p(n) >= ZERO) then
             Q_l%data(i+1,iqdens) = Q_l%data(i+1,iqdens) + &
                  beta_xp(n)*rvec_p(n,iqdens)
          endif

          if (eval_m(n) <= ZERO) then
             Q_r%data(i,iqdens) = Q_r%data(i,iqdens) + &
                  beta_xm(n)*rvec_m(n,iqdens)
          endif
       enddo

       Q_l%data(i+1,iqdens) = Qref_xp(iqdens) - Q_l%data(i+1,iqdens)
       Q_r%data(i  ,iqdens) = Qref_xm(iqdens) - Q_r%data(i  ,iqdens)


       ! velocity
       Q_l%data(i+1,iqxvel) = ZERO
       Q_r%data(i  ,iqxvel) = ZERO

       do n = 1, nwaves
          if (eval_p(n) >= ZERO) then
             Q_l%data(i+1,iqxvel) = Q_l%data(i+1,iqxvel) + &
                  beta_xp(n)*rvec_p(n,iqxvel)
          endif

          if (eval_m(n) <= ZERO) then
             Q_r%data(i,iqxvel) = Q_r%data(i,iqxvel) + &
                  beta_xm(n)*rvec_m(n,iqxvel)
          endif
       enddo

       Q_l%data(i+1,iqxvel) = Qref_xp(iqxvel) - Q_l%data(i+1,iqxvel)
       Q_r%data(i  ,iqxvel) = Qref_xm(iqxvel) - Q_r%data(i  ,iqxvel)


       ! pressure
       Q_l%data(i+1,iqpres) = ZERO
       Q_r%data(i  ,iqpres) = ZERO

       do n = 1, nwaves
          if (eval_p(n) >= ZERO) then
             Q_l%data(i+1,iqpres) = Q_l%data(i+1,iqpres) + &
                  beta_xp(n)*rvec_p(n,iqpres)
          endif

          if (eval_m(n) <= ZERO) then
             Q_r%data(i,iqpres) = Q_r%data(i,iqpres) + &
                  beta_xm(n)*rvec_m(n,iqpres)
          endif
       enddo

       Q_l%data(i+1,iqpres) = Qref_xp(iqpres) - Q_l%data(i+1,iqpres)
       Q_r%data(i,iqpres)   = Qref_xm(iqpres) - Q_r%data(i,iqpres)
       
       ! flatten
       if (use_flattening) then
          Q_l%data(i+1,:) = (ONE - xi_p%data(i,1))*Q%data(i,:) + xi_p%data(i,1)*Q_l%data(i+1,:)
          Q_r%data(i  ,:) = (ONE - xi_m%data(i,1))*Q%data(i,:) + xi_m%data(i,1)*Q_r%data(i,:)
       endif
    enddo

    ! clean-up
    call destroy(Qminus_m)
    call destroy(Qplus_m)
    call destroy(Qminus_p)
    call destroy(Qplus_p)
    call destroy(Q6_p)
    call destroy(Q6_m)

    if (use_flattening) then
       call destroy(xi_m)
       call destroy(xi_p)
    endif


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

  end subroutine make_interface_states_ppm

end module interface_states_ppm_module
