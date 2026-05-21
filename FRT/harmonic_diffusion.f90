module harmonic_diffusion

  use common

  implicit none

  private
  public &
       initialize_sgs_parameter, &
       sgs_parameter, &
       sgs_forcings, &
       finalize_sgs_parameter

  real(8):: &
       vsx, vsy, vsz, dfx, dfy

contains
  subroutine initialize_sgs_parameter

    vsx = xd*xd*vis
    vsy = yd*yd*vis
    vsz = -bst

    dfx = xd*xd*dif
    dfy = yd*yd*dif

  end subroutine initialize_sgs_parameter


  subroutine sgs_parameter

  end subroutine sgs_parameter


  subroutine sgs_forcings

    real(8):: &
         dxp,dxm,dyp,dym
    integer:: &
         i, j, k

    ! gu
    select case(bcm)
    case (frslip) ! free slip
      do k=1,km
        do j=jbi,jei
          do i=ibi,iei
            dxp= ua(i+1,j,k)-ua(i,j,k)
            dxm= ua(i-1,j,k)-ua(i,j,k)
            dyp=(ua(i,j+1,k)-ua(i,j,k))*du2(i,j+1)
            dym=(ua(i,j-1,k)-ua(i,j,k))*du2(i,j-1)
            !if ( di(i,j+1)+di(i+1,j+1)==0. ) dyp=0.
            !if ( di(i,j-1)+di(i+1,j-1)==0. ) dym=0.
            gu(i,j,k)=gu(i,j,k)+vsx*(dxp+dxm)+vsy*(dyp+dym)
          end do
        end do
      end do
    case (noslip) ! free slip
      do k=1,km
        do j=jbi,jei
          do i=ibi,iei
            dxp=ua(i+1,j,k)-ua(i,j,k)
            dxm=ua(i-1,j,k)-ua(i,j,k)
            dyp=ua(i,j+1,k)*du2(i,j+1)+ua(i,j,k)*du2(i,j+1)-2.*ua(i,j,k)
            dym=ua(i,j-1,k)*du2(i,j-1)+ua(i,j,k)*du2(i,j-1)-2.*ua(i,j,k)
            !if ( di(i,j+1)+di(i+1,j+1)==0. ) dyp=-2.*ua(i,j,k)
            !if ( di(i,j-1)+di(i+1,j-1)==0. ) dym=-2.*ua(i,j,k)
            gu(i,j,k)=gu(i,j,k)+vsx*(dxp+dxm)+vsy*(dyp+dym)
          end do
        end do
      end do
    end select
    gu(:,:,km)=gu(:,:,km)+vsz*ua(:,:,km)


    ! gv
    select case(bcm)
    case (frslip) ! free slip
      do k=1,km
        do j=jbi,jei
          do i=ibi,iei
            dxp=(va(i+1,j,k)-va(i,j,k))*dv2(i+1,j)
            dxm=(va(i-1,j,k)-va(i,j,k))*dv2(i-1,j)
            dyp=va(i,j+1,k)-va(i,j,k)
            dym=va(i,j-1,k)-va(i,j,k)
            !if ( di(i+1,j)+di(i+1,j+1)==0. ) dxp=0.
            !if ( di(i-1,j)+di(i-1,j+1)==0. ) dxm=0.
            gv(i,j,k)=gv(i,j,k)+vsx*(dxp+dxm)+vsy*(dyp+dym)
          end do
        end do
      end do
    case (noslip) ! no slip
      do k=1,km
        do j=jbi,jei
          do i=ibi,iei
            dxp=va(i+1,j,k)*dv2(i+1,j)+va(i,j,k)*dv2(i+1,j)-2.*va(i,j,k)
            dxm=va(i-1,j,k)*dv2(i-1,j)+va(i,j,k)*dv2(i-1,j)-2.*va(i,j,k)
            dyp=va(i,j+1,k)-va(i,j,k)
            dym=va(i,j-1,k)-va(i,j,k)
            !if ( di(i+1,j)+di(i+1,j+1)==0. ) dxp=-2.*va(i,j,k)
            !if ( di(i-1,j)+di(i-1,j+1)==0. ) dxm=-2.*va(i,j,k)
            gv(i,j,k)=gv(i,j,k)+vsx*(dxp+dxm)+vsy*(dyp+dym)
          end do
        end do
      end do
    end select
    gv(:,:,km)=gv(:,:,km)+vsz*va(:,:,km)


    !ge
    select case(bce)
    case (insult) ! free flux
      if ( km > 1 ) then
        do k=1,km-1
          do j=jbi,jei
            do i=ibi,iei
              dxp=(ea(i+1,j,k)-ea(i,j,k))*di(i+1,j)
              dxm=(ea(i-1,j,k)-ea(i,j,k))*di(i-1,j)
              dyp=(ea(i,j+1,k)-ea(i,j,k))*di(i,j+1)
              dym=(ea(i,j-1,k)-ea(i,j,k))*di(i,j-1)
              ge(i,j,k)=ge(i,j,k)+dfx*(dxp+dxm)+dfy*(dyp+dym)
            end do
          end do
        end do
      end if

      do k=km,km
        do j=jbi,jei
          do i=ibi,iei
            dxp=(ea(i+1,j,k)-ea(i,j,k)-dp(i+1,j)+dp(i,j))*di(i+1,j)
            dxm=(ea(i-1,j,k)-ea(i,j,k)-dp(i-1,j)+dp(i,j))*di(i-1,j)
            dyp=(ea(i,j+1,k)-ea(i,j,k)-dp(i,j+1)+dp(i,j))*di(i,j+1)
            dym=(ea(i,j-1,k)-ea(i,j,k)-dp(i,j-1)+dp(i,j))*di(i,j-1)
            ge(i,j,k)=ge(i,j,k)+dfx*(dxp+dxm)+dfy*(dyp+dym)
          end do
        end do
      end do
    case (condct)
      if ( km > 1 ) then
        do k=1,km-1
          do j=jbi,jei
            do i=ibi,iei
              dxp=(ea(i+1,j,k)-ea(i,j,k))*di(i+1,j) &
                   +2.*(hh(k)-ea(i,j,k))*(1.-di(i+1,J))
              dxm=(ea(i-1,j,k)-ea(i,j,k))*di(i-1,j) &
                   +2.*(hh(k)-ea(i,j,k))*(1.-di(i-1,j))
              dyp=(ea(i,j+1,k)-ea(i,j,k))*di(i,j+1) &
                   +2.*(hh(k)-ea(i,j,k))*(1.-di(i,j+1))
              dym=(ea(i,j-1,k)-ea(i,j,k))*di(i,j-1) &
                   +2.*(hh(k)-ea(i,j,k))*(1.-di(i,j-1))
              !if ( di(i+1,j)==0. ) dxp=2.*hh(k)-2.*ea(i,j,k)
              !if ( di(i-1,j)==0. ) dxm=2.*hh(k)-2.*ea(i,j,k)
              !if ( di(i,j+1)==0. ) dyp=2.*hh(k)-2.*ea(i,j,k)
              !if ( di(i,j-1)==0. ) dym=2.*hh(k)-2.*ea(i,j,k)
              ge(i,j,k)=ge(i,j,k)+dfx*(dxp+dxm)+dfy*(dyp+dym)
            end do
          end do
        end do
      end if
      do k=km,km
        do j=jbi,jei
          do i=ibi,iei
            dxp=(ea(i+1,j,k)-ea(i,j,k)-dp(i+1,j)+dp(i,j))*di(i+1,j) &
                 +(2.*hh(k)-2.*ea(i,j,k)-dp(i+1,j)+dp(i,j))*(1.-di(i+1,j))
            dxm=(ea(i-1,j,k)-ea(i,j,k)-dp(i-1,j)+dp(i,j))*di(i-1,j) &
                 +(2.*hh(k)-2.*ea(i,j,k)-dp(i-1,j)+dp(i,j))*(1.-di(i-1,j))
            dyp=(ea(i,j+1,k)-ea(i,j,k)-dp(i,j+1)+dp(i,j))*di(i,j+1) &
                 +(2.*hh(k)-2.*ea(i,j,k)-dp(i,j+1)+dp(i,j))*(1.-di(i,j+1))
            dym=(ea(i,j-1,k)-ea(i,j,k)-dp(i,j-1)+dp(i,j))*di(i,j-1) &
                 +(2.*hh(k)-2.*ea(i,j,k)-dp(i,j-1)+dp(i,j))*(1.-di(i,j-1))
            !if ( di(i+1,j)==0. ) dxp=2.*hh(k)-2.*ea(i,j,k)-dp(i+1,j)+dp(i,j)
            !if ( di(i-1,j)==0. ) dxm=2.*hh(k)-2.*ea(i,j,k)-dp(i-1,j)+dp(i,j)
            !if ( di(i,j+1)==0. ) dyp=2.*hh(k)-2.*ea(i,j,k)-dp(i,j+1)+dp(i,j)
            !if ( di(i,j-1)==0. ) dym=2.*hh(k)-2.*ea(i,j,k)-dp(i,j-1)+dp(i,j)
            ge(i,j,k)=ge(i,j,k)+dfx*(dxp+dxm)+dfy*(dyp+dym)
          end do
        end do
      end do
    end select


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

  end subroutine finalize_sgs_parameter


end module harmonic_diffusion
