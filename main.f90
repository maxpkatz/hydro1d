program hydro1d

  use datatypes_module
  use grid_module
  use params_module
  use probparams_module, only: init_probparams
  use init_module
  use variables_module
  use output_module
  use runtime_diag_module
  use bcs_module
  use dt_module
  use interface_states_godunov_module
  use interface_states_plm_module
  use interface_states_ppm_module
  use interface_states_ppm_temp_module
  use interface_velocities_module
  use riemann_module
  use update_module

  implicit none

  type(grid_t) :: grid
  type(gridvar_t) :: U
  type(gridedgevar_t) :: U_l, U_r, vf, fluxes

  integer :: n, ng
  real (kind=dp_t) :: t, dt, dt_new


  ! parse the inputs file
  call init_params()
  call init_probparams()

  ! set the number of ghostcells
  if (godunov_type == 0) then
     ng = 1
  else
     ng = 4
  endif


  ! build the grid and storage for grid variables, interface states,
  ! interface velocities, and fluxes
  call build(grid, nx, ng, xmin, xmax, xlboundary, xrboundary)

  call build(U, grid, ncons)
  call build(U_l, grid, ncons)
  call build(U_r, grid, ncons)

  call build(vf, grid, 1)
  
  call build(fluxes, grid, ncons)


  ! set the initial conditions
  n = 0
  t = 0.0_dp_t
  call init_data(U)

  ! set the boundary conditions
  call fillBCs(U)


  ! output the initial model
  call output(U, t, n)


  ! main evolution loop
  do while (t < tmax)

     ! set the boundary conditions
     call fillBCs(U)

     ! construct the face velocities
     call make_interface_velocities(U, vf)         
     
     ! compute the timestep
     call compute_dt(U, vf, n, dt_new)
     if (n > 0) then
        dt = min(dt_change*dt, dt_new)
     else
        dt = dt_new
     endif

     if (t + dt > tmax) dt = tmax - t

     
     ! construct the interface states
     if (godunov_type == 0) then
        call make_interface_states_godunov(U, U_l, U_r, vf, dt)
     else if (godunov_type == 1) then
        call make_interface_states_plm(U, U_l, U_r, vf, dt)
     else if (godunov_type == 2) then
        if (ppm_temp) then
           call make_interface_states_ppm_temp(U, U_l, U_r, dt)
        else
           call make_interface_states_ppm(U, U_l, U_r, dt)
        endif
     endif


     ! compute the fluxes
     call solve_riemann(U_l, U_r, vf, fluxes)


     ! do the conservative update
     call update(U, vf, fluxes, dt)

     t = t + dt
     n = n + 1


     ! output (if necessary)
     print *, n, t, dt

     if (plot_dt > 0.0_dp_t .and. &
          mod(t - dt, plot_dt) > mod(t, plot_dt)) then
        call output(U, t, n)
     endif


     ! do any runtime diagnostics
     call run_diag(U, t, n)

  enddo


  ! final output
  call output(U, t, n)


  ! clean-up
  call destroy(U)
  call destroy(grid)

end program hydro1d
