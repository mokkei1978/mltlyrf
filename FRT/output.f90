subroutine initialize_output (expid)

  use common

  implicit none

  character(*),intent(in):: &
       expid

  !call available_file_unit(backup_unit)
  open(backup_unit,file='../OUT/'//trim(expid)//'.cnt',form='unformatted')

end subroutine initialize_output


subroutine output (expid)

  use common

  implicit none

  character(*),intent(in):: &
       expid

  real(8):: &
       ca,cb
  integer:: &
       m,n
  character(80):: &
       buff

  if( istep /= 0 ) return

!  ! sum for output of averaged variable
!
!  !us(:,:,:)=us(:,:,:)+sngl(ub(:,:,:)*dt)
!  !vs(:,:,:)=vs(:,:,:)+sngl(vb(:,:,:)*dt)
!  !hs(:,:,:)=hs(:,:,:)+sngl(hb(:,:,:)*dt)
!  !ps(:,:)  =ps(:,:)  +sngl(pr(:,:)  *dt)
!  !
!  !sumdt=sumdt+sngl(dt)
!  !sumup=sumup+1
!
!  us(:,:,:)=us(:,:,:)+ub(:,:,:)*dt
!  vs(:,:,:)=vs(:,:,:)+vb(:,:,:)*dt
!  es(:,:,:)=es(:,:,:)+eb(:,:,:)*dt
!
!  sumdt=sumdt+dt
!  sumup=sumup+1

  if( time < time_to_output ) return


  m=idnint(time_to_output/output_interval)

!  us(:,:,:)=us(:,:,:)/sumdt
!  vs(:,:,:)=vs(:,:,:)/sumdt
!  es(:,:,:)=es(:,:,:)/sumdt
!
!
!  write(buff,'(a,i4.4)') '../OUT/'//trim(expid)//'.av',m
!  call available_file_unit(n)
!  open(n,file=buff,form='unformatted',access='stream')
!  write(n) nstep,sngl(time_to_output),us,vs,es,sumdt,sumup
!  close(n)

  ! time-dt < time_to_output < time
  cb=1./dt*(time_to_output-time+dt)
  ca=1./dt*(time-time_to_output)

  !us(:,:,:)=sngl(ca*ub(:,:,:)+cb*ua(:,:,:))
  !vs(:,:,:)=sngl(ca*vb(:,:,:)+cb*va(:,:,:))
  !hs(:,:,:)=sngl(ca*hb(:,:,:)+cb*ha(:,:,:))
  !ps(:,:)  =sngl(pr(:,:))

  us(:,:,:)=ca*ub(:,:,:)+cb*ua(:,:,:)
  vs(:,:,:)=ca*vb(:,:,:)+cb*va(:,:,:)
  es(:,:,:)=ca*eb(:,:,:)+cb*ea(:,:,:)

  write(buff,'(a,i4.4)') '../OUT/'//trim(expid)//'.ss',m
  open(output_unit,file=buff,form='unformatted',access='stream')
  write(output_unit) nstep,sngl(time_to_output),us,vs,es
  close(output_unit)


!  ! initialize
!  sumdt=0.
!  sumup=0
!  us(:,:,:)=0.
!  vs(:,:,:)=0.
!  es(:,:,:)=0.

!  time_to_output=time_to_output+output_interval
  time_to_output=output_interval*dble(m+1)

  if( time > time_to_output ) then
    write(6,*) 'ERR:: too short output interval'
    stop
  end if


  ! backup double precision data
  write(backup_unit) nstep,time,dt,ua,va,ea,ub,vb,eb
  rewind(backup_unit)

end subroutine output


subroutine finalize_output

  use common

  implicit none

  close(backup_unit)

end subroutine finalize_output
