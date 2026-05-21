module runge_kutta2

  use common

  implicit none

  private
  public initialize_integration_scheme, &
       time_interval, &
       time_integration, &
       swap_time_level, &
       finalize_integration_scheme

  real(8),parameter:: &
       cfl=0.1

contains


  subroutine initialize_integration_scheme

  end subroutine initialize_integration_scheme


  subroutine time_interval

    real(8):: &
         dtnew, umax, vmax, dmax, cmax
    integer:: &
         i, j, k

    select case (istep)
    case (0)

       ! update dt
       !dtnew = 2.*dt

       ! surface gravity wave
       cmax=sqrt(gr(1)*depmax)
       dtnew=cfl*min(dx,dy)/cmax

       ! current
       umax=maxval(abs(ub(ibi:iei,jbi:jei,1:km)))
       if( umax > 0. ) dtnew=min(dtnew,cfl*dx/umax)
       vmax=maxval(abs(vb(ibi:iei,jbi:jei,1:km)))
       if( vmax > 0. ) dtnew=min(dtnew,cfl*dy/vmax)

       dtnew=dble(idint(dtnew/dtmin))*dtmin
       dtnew=min(dtnew,dtmax)

       dtnew=min(dtnew,1.1d0*dt) 
       dtnew=max(dtnew,0.9d0*dt) 

       if( dtnew < dtmin ) then
          write(6,*) 'ERR:: too small time interval.'
          write(6,*) '      time, dtnew ', time, dtnew
          stop
       end if
       dt=dtnew
    case (1)
       ! dt is unchanged
    end select

  end subroutine time_interval


  subroutine time_integration

    integer:: &
         i, j, k

    select case (istep)
    case (0)

       do k=1,km
          do j=jbi,jei
             do i=ibi,iei
                ub(i,j,k)=ua(i,j,k)+0.5*dt*gu(i,j,k)
                vb(i,j,k)=va(i,j,k)+0.5*dt*gv(i,j,k)
                eb(i,j,k)=ea(i,j,k)+0.5*dt*ge(i,j,k)
             end do
          end do
       end do

       istep=1
       time=time+0.5*dt

    case (1)

       do k=1,km
          do j=jbi,jei
             do i=ibi,iei
                ub(i,j,k)=ua(i,j,k)+dt*gu(i,j,k)
                vb(i,j,k)=va(i,j,k)+dt*gv(i,j,k)
                eb(i,j,k)=ea(i,j,k)+dt*ge(i,j,k)
             end do
          end do
       end do

       nstep=nstep+1
       istep=0
       time=time+0.5*dt
    end select

!    if( maxval(abs(wb)) > 1.d2 ) then
!       write(6,*) 'ERR:: overflow'
!       stop
!    end if

  end subroutine time_integration


  subroutine swap_time_level

    if( istep /=0 ) return

    ua(:,:,:)=ub(:,:,:)
    va(:,:,:)=vb(:,:,:)
    ea(:,:,:)=eb(:,:,:)

  end subroutine swap_time_level


  subroutine finalize_integration_scheme

  end subroutine finalize_integration_scheme
end module runge_kutta2
