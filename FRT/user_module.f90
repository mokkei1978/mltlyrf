module user_module

  use common
  implicit none

  private
  public &
       user_configuration, &
       user_initial_condition, &
       user_forcings, &
       user_action, &
       user_boundary_condition, &
       user_output


  real(8),allocatable:: &
       tx(:,:), ty(:,:) ! wind stress

  real(8):: &
       rho0, &
       tx_yamp, &
       tx_ylen, &
       tx_tamp, &
       tx_tper

  integer:: &
       user_sumup 

  integer:: &
       i, j, k, n

  real(8):: &
       user_sumdt

  logical:: &
       linear


contains
  subroutine user_configuration (expid)

    character(*),intent(in):: &
         expid
    real(8):: &
         x, y, xc, yc, lx, ly
    integer:: &
         i,j
    
    namelist &
         /user_parameter/ &
         rho0, &
         tx_yamp, &
         tx_ylen, &
         tx_tamp, &
         tx_tper


    ! allocation
    allocate( &
         tx(ib:ie,jb:je), &
         ty(ib:ie,jb:je) &
         )

    ! set default values

    ! read conficuration
    call available_file_unit(n)
    open(n,file='../CNF/'//trim(expid)//'.cnf')
    read(n,user_parameter)
    close(n)

    ty(:,:)=0.
    
    xc=100.*dx
    yc=dble(jmi)*dy*0.5

    lx=10.*dx
    ly=10.*dx

    select case (expid(1:8))
    case("seamount")
       ! topography
       do j=ib,je
          y=dy*(dble(j)-0.5)
          do i=ib,ie
             x=dx*(dble(i)-0.5)

             dp(i,j)=sum(hh(:)) &
                  -2000.*dexp(-(x-xc)*(x-xc)/lx/lx-(y-yc)*(y-yc)/ly/ly)
          end do
       end do
    end select

  end subroutine user_configuration


  subroutine user_initial_condition

    user_sumdt=0.
    user_sumup=0

  end subroutine user_initial_condition

  
  subroutine user_forcings

    real(8):: &
         y, uv

    ! wind stress
    do j=jb,je
       y=(dble(j-jb)-0.5)*dy
       tx(:,j)=-tx_yamp*dcos(2.*dpi*y/tx_ylen)+tx_tamp*dsin(2.*dpi*time/tx_tper)
       !write(6,*) j,tx(1,j)
    end do
!    stop

    do j=jbi,jei
      do i=ibi,iei
        gu(i,j,1)=gu(i,j,1) &
             +tx(i,j)/rho0/hh(1)
        gv(i,j,1)=gv(i,j,1) &
             +ty(i,j)/rho0/hh(1)
      end do
    end do

  end subroutine user_forcings


  subroutine user_action

  end subroutine user_action


  subroutine user_boundary_condition

  end subroutine user_boundary_condition


  subroutine user_output (expid)

    character(*),intent(in):: &
         expid

  end subroutine user_output


end module user_module
