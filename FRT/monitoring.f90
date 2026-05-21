subroutine initialize_monitoring (expid)

  use common

  implicit none

  character(*),intent(in):: &
       expid

  !call available_file_unit(monitor_unit)

  if( time == 0. ) then
     open(monitor_unit,file='../LOG/'//expid//'.log',form='formatted')
  else
     open(monitor_unit,file='../LOG/'//expid//'.log',form='formatted',&
          status='old',position='append')     
  end if

  write(monitor_unit,'(a)') '  time        step    dt'// &
       '   umax      i   j   k'// &
       '   vmax      i   j   k'// &
       '   emax      i   j   k'!// &
!       '   Hmax            i   j'

end subroutine initialize_monitoring


subroutine monitoring

  use common

  implicit none

  integer,dimension(:):: &
       lmax(3)

  real(8):: &
       vmax1,vmax2 

  character(*),parameter:: &
       fmt='(a,es12.3,a,3i4)'
  character(196):: &
       line

  if( istep /= 0 ) return
  if( time < time_to_monitor ) return

  !write(line,'(f10.2,i8,f6.1)') time/86400., nstep, dt
  write(line,'(es10.2,i8,f6.1)') time, nstep, dt

  ! u
  vmax1=maxval(abs(ub(ibi:iei,jbi:jei,1:km)))
  lmax =maxloc(abs(ub(ibi:iei,jbi:jei,1:km)))
  vmax2=maxval(    ub(ibi:iei,jbi:jei,1:km))
  if( vmax1 /= vmax2 ) vmax1=-vmax1

  write(line,'(a,es10.1,3i4)') trim(line), vmax1, lmax

  ! v
  vmax1=maxval(abs(vb(ibi:iei,jbi:jei,1:km)))
  lmax =maxloc(abs(vb(ibi:iei,jbi:jei,1:km)))
  vmax2=maxval(    vb(ibi:iei,jbi:jei,1:km))
  if( vmax1 /= vmax2 ) vmax1=-vmax1

  write(line,'(a,es10.1,3i4)') trim(line), vmax1, lmax

  ! e
  vmax1=minval(eb(ibi:iei,jbi:jei,1:km))
  lmax =minloc(eb(ibi:iei,jbi:jei,1:km))

  write(line,'(a,es10.1,3i4)') trim(line), vmax1, lmax

  write(monitor_unit,'(a)') trim(line)
  do while ( time_to_monitor <= time )
     time_to_monitor=time_to_monitor+monitor_interval
  end do


end subroutine monitoring


subroutine finalize_monitoring (cputim)

  use common

  implicit none

  real(4):: &
       cputim

  write(monitor_unit,*) 'MSG:: total CPU time ',cputim
  close(monitor_unit)

end subroutine finalize_monitoring
