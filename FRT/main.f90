!     Multi Layer Shallow Water Equation Model with Free Surfaces.
!     Nonlinear Initial Value Problem is Solved.
!     2nd Order Central Difference Scheme is Used.
!     Arakawa Enstrophy Conservation Scheme is Used for Momentum Advection.

!     Unit in MKS
!     OCT 2007 Y.Yoshikawa
!     FEB 2009 Y.Yoshikawa
!     APR 2022 Y.Yoshikawa

program main

  use common
!  use euler_asselin
  use runge_kutta2
!  use harmonic_diffusion
  use biharmonic_diffusion
  use user_module


  implicit none

  character(*),parameter:: &
       !expid='flat2case2'
       !expid='seamount3case2'
       expid='seamount2case2'

  real(4):: &
       cputim1, cputim2

  ! initialization processes
  call cpu_time(cputim1)
  call configuration (expid)
  call configuration_check
  call user_configuration (expid)
  call useful_constants

  call initialize_integration_scheme
  call initialize_sgs_parameter
  call initialize_output (expid)
  call initialize_monitoring (expid)

  call initial_condition (expid)
  call user_initial_condition
  call user_output (expid)
  call output (expid)
  
  ! main processes
  do
     ! calculate diagnostic variables
     call diagnostic
     call sgs_parameter
     call time_interval

     ! calculate RHS of governng equations
     call inertial_forcings ! this must be called at first
     call user_forcings
     call sgs_forcings ! this must be called at last

     ! user's action
     call user_action
     ! calculate next step variables
     call time_integration
     call boundary_condition
     call user_boundary_condition

     ! post process
     call user_output (expid)
     call output (expid)
     call monitoring

     call swap_time_level
     
     if( istep/=0 ) cycle
     if( time >= time_to_end ) exit
  end do

  ! finalization processes
  call finalize_sgs_parameter
  call finalize_integration_scheme
  call finalize_output
  call cpu_time (cputim2)
  call finalize_monitoring (cputim2-cputim1)

end program main
