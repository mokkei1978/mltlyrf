module biharmonic_diffusion

  use common

  implicit none

  private
  public &
       initialize_sgs_parameter, &
       sgs_parameter, &
       sgs_forcings, &
       finalize_sgs_parameter
       
  real(8),allocatable,public:: &
       dm1(:,:,:), dm2(:,:,:)

  real(8) :: &
       vsx, vsy, vsz, dfx, dfy

contains
  subroutine initialize_sgs_parameter

    allocate( dm1(ib:ie,jb:je,1:km), dm2(ib:ie,jb:je,1:km) )

    vsx = -xd*xd*xd*xd*vis
    vsy = -yd*yd*yd*yd*vis
    vsz = -bst

    dfx = -xd*xd*xd*xd*dif
    dfy = -yd*yd*yd*yd*dif

    dm1(:,:,:)=0.
    dm2(:,:,:)=0.

  end subroutine initialize_sgs_parameter

  subroutine sgs_parameter

  end subroutine sgs_parameter

  subroutine sgs_forcings

    real(8):: &
         dxp, dxm, dyp, dym
    integer:: &
         i, j, k

    ! gu
    ! d^2 u/dx^2
    ! laplacian = 0 at the boundary (symmetric profile)
    ! laplacian u at the boundary is used for laplacian laplacian u
    do k=1,km
      do j=jbi,jei
        do i=ibi,iei
          dxp=ua(i+1,j,k)-ua(i,j,k)
          dxm=ua(i-1,j,k)-ua(i,j,k)
          dm1(i,j,k)=(dxp+dxm)*du(i,j)
        end do
      end do
    end do
    ! boundary condition
    select case (bcx)
    case (cyclic) ! periodic
      dm1(ibi-1,:,:)=dm1(iei,:,:)
      dm1(iei+1,:,:)=dm1(ibi,:,:)
    case (wall) ! laplacian = 0 at the boundary (symmetric profile)
      dm1(ibi-1,:,:)=0.
      dm1(iei+1,:,:)=0.
    end select

    ! d^4 u/dx^4
    do k=1,km
      do j=jbi,jei
        do i=ibi,iei
          dxp=dm1(i+1,j,k)-dm1(i,j,k)
          dxm=dm1(i-1,j,k)-dm1(i,j,k)
          dm2(i,j,k)=(dxp+dxm)*du(i,j)
        end do
      end do
    end do
    ! boundary condition
    select case (bcx)
    case (cyclic) ! periodic
      dm2(ibi-1,:,:)=dm2(iei,:,:)
      dm2(iei+1,:,:)=dm2(ibi,:,:)
    case (wall) ! unused
      dm2(ibi-1,:,:)=0.
      dm2(iei+1,:,:)=0.
    end select

    gu(:,:,:)=gu(:,:,:)+vsx*dm2(:,:,:)


    select case(bcm)
    case (frslip) ! free slip
      ! d^2 u/dy^2
      do k=1,km
        do j=jbi,jei
          do i=ibi,iei
            dyp=(ua(i,j+1,k)-ua(i,j,k))*du2(i,j+1)
            dym=(ua(i,j-1,k)-ua(i,j,k))*du2(i,j-1)
            !if ( di(i,j+1)+di(i+1,j+1)==0. ) dyp=0.
            !if ( di(i,j-1)+di(i+1,j-1)==0. ) dym=0.
            dm1(i,j,k)=dyp+dym
          end do
        end do
      end do
      ! boundary condition
      select case (bcy)
      case (cyclic) ! periodic
        dm1(:,jbi-1,:)= dm1(:,jei,:)
        dm1(:,jei+1,:)= dm1(:,jbi,:)
      case (wall) ! symmetric profile
        dm1(:,jbi-1,:)= dm1(:,jbi,:)
        dm1(:,jei+1,:)= dm1(:,jei,:)
      end select
      ! d^4 u/dy^4
      do k=1,km
        do j=jbi,jei
          do i=ibi,iei
            dyp=(dm1(i,j+1,k)-dm1(i,j,k))*du2(i,j+1)
            dym=(dm1(i,j-1,k)-dm1(i,j,k))*du2(i,j-1)
            !if ( di(i,j+1)+di(i+1,j+1)==0. ) dyp=0.
            !if ( di(i,j-1)+di(i+1,j-1)==0. ) dym=0.
            dm2(i,j,k)=dyp+dym
          end do
        end do
      end do
      ! boundary condition
      select case (bcy)
      case (cyclic) ! periodic
        dm2(:,jbi-1,:)=dm2(:,jei,:)
        dm2(:,jei+1,:)=dm2(:,jbi,:)
      case (wall) ! unused
        dm2(:,jbi-1,:)=0.
        dm2(:,jei+1,:)=0.
      end select
    case (noslip) ! no slip
      ! d^2 u/dy^2
      do k=1,km
        do j=jbi,jei
          do i=ibi,iei
            dyp=ua(i,j+1,k)*du2(i,j+1)+ua(i,j,k)*du2(i,j+1)-2.*ua(i,j,k)
            dym=ua(i,j-1,k)*du2(i,j-1)+ua(i,j,k)*du2(i,j-1)-2.*ua(i,j,k)
            !if ( di(i,j+1)+di(i+1,j+1)==0. ) dyp=-2.*ua(i,j,k)
            !if ( di(i,j-1)+di(i+1,j-1)==0. ) dym=-2.*ua(i,j,k)
            dm1(i,j,k)=dyp+dym
          end do
        end do
      end do
      ! boundary condition
      select case (bcy)
      case (cyclic) ! periodic
        dm1(:,jbi-1,:)= dm1(:,jei,:)
        dm1(:,jei+1,:)= dm1(:,jbi,:)
      case (wall) ! anti-symmetric profile
        dm1(:,jbi-1,:)=-dm1(:,jbi,:)
        dm1(:,jei+1,:)=-dm1(:,jei,:)
      end select
      ! d^4 u/dy^4
      do k=1,km
        do j=jbi,jei
          do i=ibi,iei
            dyp=dm1(i,j+1,k)*du2(i,j+1)+dm1(i,j,k)*du2(i,j+1)-2.*dm1(i,j,k)
            dym=dm1(i,j-1,k)*du2(i,j-1)+dm1(i,j,k)*du2(i,j-1)-2.*dm1(i,j,k)
            !if ( di(i,j+1)+di(i+1,j+1)==0. ) dyp=-2.*dm1(i,j,k)
            !if ( di(i,j-1)+di(i+1,j-1)==0. ) dym=-2.*dm1(i,j,k)
            dm2(i,j,k)=dyp+dym
          end do
        end do
      end do
      ! boundary condition
      select case (bcy)
      case (cyclic) ! periodic
        dm2(:,jbi-1,:)=dm2(:,jei,:)
        dm2(:,jei+1,:)=dm2(:,jbi,:)
      case (wall) ! unused
        dm2(:,jbi-1,:)=0.
        dm2(:,jei+1,:)=0.
      end select
    end select
    gu(:,:,:)=gu(:,:,:)+vsy*dm2(:,:,:)

    ! bottom stress
    gu(:,:,km)=gu(:,:,km)+vsz*ua(:,:,km)


    ! gv
    select case (bcm)
    case (frslip) ! free slip
      ! d^2 v/dx^2
      do k=1,km
        do j=jbi,jei
          do i=ibi,iei
            dxp=(va(i+1,j,k)-va(i,j,k))*dv2(i+1,j)
            dxm=(va(i-1,j,k)-va(i,j,k))*dv2(i-1,j)
            !if ( di(i+1,j)+di(i+1,j+1)==0. ) dxp=0.
            !if ( di(i-1,j)+di(i-1,j+1)==0. ) dxm=0.
            dm1(i,j,k)=dxp+dxm
          end do
        end do
      end do
      ! boundary condition
      select case (bcx)
      case (cyclic) ! periodic
        dm1(ibi-1,:,:)= dm1(iei,:,:)
        dm1(iei+1,:,:)= dm1(ibi,:,:)
      case (wall)
        dm1(ibi-1,:,:)= dm1(ibi,:,:)
        dm1(iei+1,:,:)= dm1(iei,:,:)
      end select
      ! d^4 v/dx^4
      do k=1,km
        do j=jbi,jei
          do i=ibi,iei
            dxp=(dm1(i+1,j,k)-dm1(i,j,k))*dv2(i+1,j)
            dxm=(dm1(i-1,j,k)-dm1(i,j,k))*dv2(i-1,j)
            !if ( di(i+1,j)+di(i+1,j+1)==0. ) dxp=0.
            !if ( di(i-1,j)+di(i-1,j+1)==0. ) dxm=0.
            dm2(i,j,k)=dxp+dxm
          end do
        end do
      end do
      ! boundary condition
      select case (bcx)
      case (cyclic) ! periodic
        dm2(ibi-1,:,:)=dm2(iei,:,:)
        dm2(iei+1,:,:)=dm2(ibi,:,:)
      case (wall) ! unused
        dm2(ibi-1,:,:)=0.
        dm2(iei+1,:,:)=0.
      end select
    case (noslip)
      ! d^2 v/dx^2
      do k=1,km
        do j=jbi,jei
          do i=ibi,iei
            dxp=va(i+1,j,k)*dv2(i+1,j)+va(i,j,k)*dv2(i+1,j)-2.*va(i,j,k)
            dxm=va(i-1,j,k)*dv2(i-1,j)+va(i,j,k)*dv2(i-1,j)-2.*va(i,j,k)
            !if ( di(i+1,j)+di(i+1,j+1)==0. ) dxp=-2.*va(i,j,k)
            !if ( di(i-1,j)+di(i-1,j+1)==0. ) dxm=-2.*va(i,j,k)
            dm1(i,j,k)=dxp+dxm
          end do
        end do
      end do
      ! boundary condition
      select case (bcx)
      case (cyclic) ! periodic
        dm1(ibi-1,:,:)= dm1(iei,:,:)
        dm1(iei+1,:,:)= dm1(ibi,:,:)
      case (wall)
        dm1(ibi-1,:,:)=-dm1(ibi,:,:)
        dm1(iei+1,:,:)=-dm1(iei,:,:)
      end select
      ! d^4 v/dx^4
      do k=1,km
        do j=jbi,jei
          do i=ibi,iei
            dxp=dm1(i+1,j,k)*dv2(i+1,j)+dm1(i,j,k)*dv2(i+1,j)-2.*dm1(i,j,k)
            dxm=dm1(i-1,j,k)*dv2(i-1,j)+dm1(i,j,k)*dv2(i-1,j)-2.*dm1(i,j,k)
            !if ( di(i+1,j)+di(i+1,j+1)==0. ) dxp=-2.*dm1(i,j,k)
            !if ( di(i-1,j)+di(i-1,j+1)==0. ) dxm=-2.*dm1(i,j,k)
            dm2(i,j,k)=dxp+dxm
          end do
        end do
      end do
      ! boundary condition
      select case (bcx)
      case (cyclic) ! periodic
        dm2(ibi-1,:,:)=dm2(iei,:,:)
        dm2(iei+1,:,:)=dm2(ibi,:,:)
      case (wall) ! unused
        dm2(ibi-1,:,:)=0.
        dm2(iei+1,:,:)=0.
      end select
    end select
    gv(:,:,:)=gv(:,:,:)+vsx*dm2(:,:,:)


    ! d^2 v/dy^2
    do k=1,km
      do j=jbi,jei
        do i=ibi,iei
          dyp=va(i,j+1,k)-va(i,j,k)
          dym=va(i,j-1,k)-va(i,j,k)
          dm1(i,j,k)=(dyp+dym)*dv(i,j)
        end do
      end do
    end do
    ! boundary condition
    select case (bcy)
    case (cyclic) ! periodic
      dm1(:,jbi-1,:)=dm1(:,jei,:)
      dm1(:,jei+1,:)=dm1(:,jbi,:)
    case (wall) ! laplacian = 0 at the boundary (symmetric profile)
      dm1(:,jbi-1,:)=0.
      dm1(:,jei+1,:)=0.
    end select

    ! d^4 v/dy^4
    do k=1,km
      do j=jbi,jei
        do i=ibi,iei
          dyp=dm1(i,j+1,k)-dm1(i,j,k)
          dym=dm1(i,j-1,k)-dm1(i,j,k)
          dm2(i,j,k)=(dyp+dym)*dv(i,j)
        end do
      end do
    end do
    ! boundary condition
    select case (bcy)
    case (cyclic) ! periodic
      dm2(:,jbi-1,:)=dm2(:,jei,:)
      dm2(:,jei+1,:)=dm2(:,jbi,:)
    case (wall) ! unused
      dm2(:,jbi-1,:)=0.
      dm2(:,jei+1,:)=0.
    end select

    gv(:,:,:)=gv(:,:,:)+vsy*dm2(:,:,:)

    ! bottom stress
    gv(:,:,km)=gv(:,:,km)+vsz*va(:,:,km)



!    !ge
!    select case(bce)
!    case (insult) ! free flux
!      ! d^2 h/dx^2
!      if ( km > 1 ) then
!        do k=1,km-1
!          do j=jbi,jei
!            do i=ibi,iei
!              dxp=(ea(i+1,j,k)-ea(i,j,k))*di(i+1,j)
!              dxm=(ea(i-1,j,k)-ea(i,j,k))*di(i-1,j)
!              !if ( di(i+1,j)==0. ) dxp=0.
!              !if ( di(i-1,j)==0. ) dxm=0.
!              dm1(i,j,k)=dxp+dxm
!            end do
!          end do
!        end do
!      end if
!    
!      do k=km,km
!        do j=jbi,jei
!          do i=ibi,iei
!            dxp=(ea(i+1,j,k)-ea(i,j,k))*di(i+1,j)
!            dxm=(ea(i-1,j,k)-ea(i,j,k))*di(i-1,j)
!            !if ( di(i+1,j)==0. ) dxp=0.
!            !if ( di(i-1,j)==0. ) dxm=0.
!            dm1(i,j,k)=dxp+dxm
!          end do
!        end do
!      end do
!      ! boundary condition
!      select case (bcx)
!      case (cyclic) ! periodic
!        dm1(ibi-1,:,:)=dm1(iei,:,:)
!        dm1(iei+1,:,:)=dm1(ibi,:,:)
!      case (imperm) ! unused
!        dm1(ibi-1,:,:)= dm1(ibi,:,:)
!        dm1(iei+1,:,:)= dm1(iei,:,:)
!      end select
!
!      ! d^4h/dx^4
!      do k=1,km
!        do j=jbi,jei
!          do i=ibi,iei
!            dxp=(dm1(i+1,j,k)-dm1(i,j,k))*di(i+1,j)
!            dxm=(dm1(i-1,j,k)-dm1(i,j,k))*di(i-1,j)
!            !if ( di(i+1,j)==0. ) dxp=0.
!            !if ( di(i-1,j)==0. ) dxm=0.
!            dm2(i,j,k)=dxp+dxm
!          end do
!        end do
!      end do
!      ! boundary condition
!      select case (bcx)
!      case (cyclic) ! periodic
!        dm2(ibi-1,:,:)=dm2(iei,:,:)
!        dm2(iei+1,:,:)=dm2(ibi,:,:)
!      case (imperm) ! unused
!        dm2(ibi-1,:,:)=0.
!        dm2(iei+1,:,:)=0.
!      end select
!    case (condct)
!      ! d^2 h/dx^2
!      if ( km > 1 ) then
!        do k=1,km-1
!          do j=jbi,jei
!            do i=ibi,iei
!              dxp=(ea(i+1,j,k)-ea(i,j,k))*di(i+1,j) &
!                  +2.*(hh(k)-ea(i,j,k))*(1.-di(i+1,j))
!              dxm=(ea(i-1,j,k)-ea(i,j,k))*di(i-1,j) &
!                  +2.*(hh(k)-ea(i,j,k))*(1.-di(i-1,j))
!              !if ( di(i+1,j)==0. ) dxp=2.*hh(k)-2.*ea(i,j,k)
!              !if ( di(i-1,j)==0. ) dxm=2.*hh(k)-2.*ea(i,j,k)
!              dm1(i,j,k)=dxp+dxm
!            end do
!          end do
!        end do
!      end if
!      do k=km,km
!        do j=jbi,jei
!          do i=ibi,iei
!            dxp=(ea(i+1,j,k)-ea(i,j,k))*di(i+1,j) &
!               +(2.*hh(k)-2.*ea(i,j,k))*(1.-di(i+1,j))
!            dxm=(ea(i-1,j,k)-ea(i,j,k))*di(i-1,j) &
!               +(2.*hh(k)-2.*ea(i,j,k))*(1.-di(i-1,j))
!            !if ( di(i+1,j)==0. ) dxp=-2.*ea(i,j,k)
!            !if ( di(i-1,j)==0. ) dxm=-2.*ea(i,j,k)
!            dm1(i,j,k)=dxp+dxm
!          end do
!        end do
!      end do
!    ! boundary condition
!      select case (bcx)
!      case (cyclic) ! periodic
!        dm1(ibi-1,:,:)=dm1(iei,:,:)
!        dm1(iei+1,:,:)=dm1(ibi,:,:)
!      case (imperm) ! unused
!        dm1(ibi-1,:,:)=-dm1(ibi,:,:)
!        dm1(iei+1,:,:)=-dm1(iei,:,:)
!      end select
!      ! d^4h/dx^4
!      do k=1,km
!        do j=jbi,jei
!          do i=ibi,iei
!            dxp=dm1(i+1,j,k)*di(i+1,j)+dm1(i,j,k)*di(i+1,j)-2.*dm1(i,j,k)
!            dxm=dm1(i-1,j,k)*di(i-1,j)+dm1(i,j,k)*di(i-1,j)-2.*dm1(i,j,k)
!            !if ( di(i+1,j)==0. ) dxp=-2.*dm1(i,j,k)
!            !if ( di(i-1,j)==0. ) dxm=-2.*dm1(i,j,k)
!            dm2(i,j,k)=dxp+dxm
!          end do
!        end do
!      end do
!      ! boundary condition
!      select case (bcx)
!      case (cyclic) ! periodic
!        dm2(ibi-1,:,:)=dm2(iei,:,:)
!        dm2(iei+1,:,:)=dm2(ibi,:,:)
!      case (imperm) ! unused
!        dm2(ibi-1,:,:)=0.
!        dm2(iei+1,:,:)=0.
!      end select
!    end select
!
!    ge(:,:,:)=ge(:,:,:)+dfx*dm2(:,:,:)
!    
!    
!
!    select case(bce)
!    case (insult) ! free flux
!      ! d^2 h/dy^2
!      if ( km > 1 ) then
!        do k=1,km-1
!          do j=jbi,jei
!            do i=ibi,iei
!              dyp=(ea(i,j+1,k)-ea(i,j,k))*di(i,j+1)
!              dym=(ea(i,j-1,k)-ea(i,j,k))*di(i,j-1)
!              !if ( di(i,j+1)==0. ) dyp=0.
!              !if ( di(i,j-1)==0. ) dym=0.
!              dm1(i,j,k)=dyp+dym
!            end do
!          end do
!        end do
!      end if
!    
!      do k=km,km
!        do j=jbi,jei
!          do i=ibi,iei
!            dyp=(ea(i,j+1,k)-ea(i,j,k))*di(i,j+1)
!            dym=(ea(i,j-1,k)-ea(i,j,k))*di(i,j-1)
!            !if ( di(i,j+1)==0. ) dyp=0.
!            !if ( di(i,j-1)==0. ) dym=0.
!            dm1(i,j,k)=dyp+dym
!          end do
!        end do
!      end do
!      ! boundary condition
!      select case (bcy)
!      case (cyclic) ! periodic
!        dm1(:,jbi-1,:)=dm1(:,jei,:)
!        dm1(:,jei+1,:)=dm1(:,jbi,:)
!      case (imperm) ! symmetric
!        dm1(:,jbi-1,:)= dm1(:,jbi,:)
!        dm1(:,jei+1,:)= dm1(:,jei,:)
!      end select
!      ! d^4h/dy^4
!      do k=1,km
!        do j=jbi,jei
!          do i=ibi,iei
!            dyp=(dm1(i,j+1,k)-dm1(i,j,k))*di(i,j+1)
!            dym=(dm1(i,j-1,k)-dm1(i,j,k))*di(i,j-1)
!            !if ( di(i,j+1)==0. ) dyp=0.
!            !if ( di(i,j-1)==0. ) dym=0.
!            dm2(i,j,k)=dyp+dym
!          end do
!        end do
!      end do
!      ! boundary condition
!      select case (bcy)
!      case (cyclic) ! periodic
!        dm2(:,jbi-1,:)=dm2(:,jei,:)
!        dm2(:,jei+1,:)=dm2(:,jbi,:)
!      case (imperm) ! unused
!        dm2(:,jbi-1,:)=0.
!        dm2(:,jei+1,:)=0.
!      end select
!    case (condct)
!      ! d^2 h/dy^2
!      if ( km > 1 ) then
!        do k=1,km-1
!          do j=jbi,jei
!            do i=ibi,iei
!              dyp=(ea(i,j+1,k)-ea(i,j,k))*di(i,j+1) &
!                 +2.*(hh(k)-ea(i,j,k))*(1.-di(i,j+1))
!              dym=(ea(i,j-1,k)-ea(i,j,k))*di(i,j-1) &
!                 +2.*(hh(k)-ea(i,j,k))*(1.-di(i,j-1))
!              !if ( di(i,j+1)==0. ) dyp=2.*hh(k)-2.*ea(i,j,k)
!              !if ( di(i,j-1)==0. ) dym=2.*hh(k)-2.*ea(i,j,k)
!              dm1(i,j,k)=dyp+dym
!            end do
!          end do
!        end do
!      end if
!      do k=km,km
!        do j=jbi,jei
!          do i=ibi,iei
!            dyp=(ea(i,j+1,k)-ea(i,j,k))*di(i,j+1) &
!               +(2.*hh(k)-2.*ea(i,j,k))*(1.-di(i,j+1))
!            dym=(ea(i,j-1,k)-ea(i,j,k))*di(i,j-1) &
!               +(2.*hh(k)-2.*ea(i,j,k))*(1.-di(i,j-1))
!            !if ( di(i,j+1)==0. ) dxp=-2.*ea(i,j,k)
!            !if ( di(i,j-1)==0. ) dxm=-2.*ea(i,j,k)
!            dm1(i,j,k)=dyp+dym
!          end do
!        end do
!      end do
!      ! boundary condition
!      select case (bcy)
!      case (cyclic) ! periodic
!        dm1(:,jbi-1,:)=dm1(:,jei,:)
!        dm1(:,jei+1,:)=dm1(:,jbi,:)
!      case (imperm) ! unused
!        dm1(:,jbi-1,:)=-dm1(:,jbi,:)
!        dm1(:,jei+1,:)=-dm1(:,jei,:)
!      end select
!      ! d^4h/dy^4
!      do k=1,km
!        do j=jbi,jei
!          do i=ibi,iei
!            dyp=dm1(i,j+1,k)*di(i,j+1)+dm1(i,j,k)*di(i,j+1)-2.*dm1(i,j,k)
!            dym=dm1(i,j-1,k)*di(i,j-1)+dm1(i,j,k)*di(i,j-1)-2.*dm1(i,j,k)
!            !if ( di(i,j+1)==0. ) dyp=-2.*dm1(i,j,k)
!            !if ( di(i,j-1)==0. ) dym=-2.*dm1(i,j,k)
!            dm2(i,j,k)=dyp+dym
!          end do
!        end do
!      end do
!      ! boundary condition
!      select case (bcy)
!      case (cyclic) ! periodic
!        dm2(:,jbi-1,:)=dm2(:,jei,:)
!        dm2(:,jei+1,:)=dm2(:,jbi,:)
!      case (imperm) ! unused
!        dm2(:,jbi-1,:)=0.
!        dm2(:,jei+1,:)=0.
!      end select
!    end select
!    ge(:,:,:)=ge(:,:,:)+dfy*dm2(:,:,:)


    ! mask
    do k=1,km
      gu(:,:,k)=gu(:,:,k)*du(:,:)
      gv(:,:,k)=gv(:,:,k)*dv(:,:)
      ge(:,:,k)=ge(:,:,k)*di(:,:)
    end do

    ! boundary condition for periodic condition
    select case (bcx)
    case (cyclic) ! periodic
      gu(ibi-1,:,:)=gu(iei,:,:)
      gu(iei+1,:,:)=gu(ibi,:,:)
      gv(ibi-1,:,:)=gv(iei,:,:)
      gv(iei+1,:,:)=gv(ibi,:,:)
      ge(ibi-1,:,:)=ge(iei,:,:)
      ge(iei+1,:,:)=ge(ibi,:,:)
    case (wall)
      !gu(ie,:,:)=0.
    end select
    select case (bcy)
    case (cyclic) ! periodic
      gu(:,jbi-1,:)=gu(:,jei,:)
      gu(:,jei+1,:)=gu(:,jbi,:)
      gv(:,jbi-1,:)=gv(:,jei,:)
      gv(:,jei+1,:)=gv(:,jbi,:)
      ge(:,jbi-1,:)=ge(:,jei,:)
      ge(:,jei+1,:)=ge(:,jbi,:)
    case (wall)
      !gv(:,je,:)=0.
    end select

  end subroutine sgs_forcings


  subroutine finalize_sgs_parameter

    deallocate( dm1, dm2 )

  end subroutine finalize_sgs_parameter

end module biharmonic_diffusion
