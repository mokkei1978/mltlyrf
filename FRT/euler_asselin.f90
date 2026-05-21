module euler_asselin

  use common

  implicit none

  private
  public initialize_integration_scheme, &
       time_interval, &
       time_integration, &
       swap_time_level, &
       finalize_integration_scheme

  real(8):: &
       asf ! coefficient for asselin filter

contains
  subroutine initialize_integration_scheme

    asf=0.5d0

  end subroutine initialize_integration_scheme


  subroutine time_interval

    ! dt should be constant in this Euler scheme

  end subroutine time_interval


  subroutine time_integration

    integer:: &
         i, j, k

    ! update GU  <- n+1 step variable
    ! GU=u^(n+1)
    do k=1,km
       do j=jbi,jei
          do i=ibi,iei
             GU(i,j,k)=ua(i,j,k)+2.*dt*gu(i,j,k)
             GV(i,j,k)=va(i,j,k)+2.*dt*gv(i,j,k)
             GE(i,j,k)=ea(i,j,k)+2.*dt*ge(i,j,k)
          end do
       end do
    end do

    ! asselin filter is
    ! u^(n)=u^(n)+0.5*asf*(u^(n-1)-2*u^(n)+u^(n+1)).
    ! this filter is divided into two processes
    ! u^*=u^(n)+0.5*asf*(u^(n-1)-2*u^(n))
    ! u^(n)=u^*+0.5*asf*(u^(n+1))
    ! the first process is performed in this subroutine
    ! the second process is performed in swap_time_level


    ! update ua  <- n step variable (with asselin filter (first process))
    ! ua=u^*
    UA(:,:,:)=ub(:,:,:)+0.5*asf*(ua(:,:,:)-2.*ub(:,:,:))
    VA(:,:,:)=vb(:,:,:)+0.5*asf*(va(:,:,:)-2.*vb(:,:,:))
    EA(:,:,:)=eb(:,:,:)+0.5*asf*(ea(:,:,:)-2.*eb(:,:,:))

    ! update ub  <- n+1 step variable
    ! ub=u^(n+1)
    ub(:,:,:)=gu(:,:,:)
    vb(:,:,:)=gv(:,:,:)
    eb(:,:,:)=ge(:,:,:)

    nstep=nstep+1
    time =time+dt

  end subroutine time_integration


  subroutine swap_time_level
    ! ua=u^(n)
    ! ua=u^*+0.5*asf*(u^(n+1)) (asselin filter (second half process))
    ua(:,:,:)=UA(:,:,:)+0.5*asf*ub(:,:,:)
    va(:,:,:)=VA(:,:,:)+0.5*asf*vb(:,:,:)
    ea(:,:,:)=EA(:,:,:)+0.5*asf*eb(:,:,:)

  end subroutine swap_time_level


  subroutine finalize_integration_scheme

  end subroutine finalize_integration_scheme

end module euler_asselin
