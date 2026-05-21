subroutine initial_condition (expid)

  use common

  implicit none

  character(*),intent(in):: &
       expid

  if( time_to_start == 0. ) then
     nstep=0
     istep=0
     time =0.
     dt=dtmin

!     if( inityp == 'user' ) then
!        call user_initial_condition
!     else
     ub(:,:,:)=inival_u
     vb(:,:,:)=inival_v
     eb(:,:,:)=inival_e
!     end if

     call boundary_condition

     ua(:,:,:)=ub(:,:,:)
     va(:,:,:)=vb(:,:,:)
     ea(:,:,:)=eb(:,:,:)

     time_to_output =0.
     time_to_monitor=0.

     sumdt=0.
     sumup=0
     us(:,:,:)=0.
     vs(:,:,:)=0.
     es(:,:,:)=0.

!     call user_output (expid)
!     call output (expid)
  else
     read(backup_unit) nstep,time,dt,ua,va,ea,ub,vb,eb
     rewind (backup_unit)

     time_to_output =dble(idint(time/ output_interval)+1)* output_interval
     time_to_monitor=dble(idint(time/monitor_interval)+1)*monitor_interval

     sumdt=0.
     sumup=0
  endif

end subroutine initial_condition
