subroutine diagnostic

  use common

  implicit none

  real(8):: &
!       uy, vx, &
       hi!, findx
  integer:: &
       i, j, k

  fu(:,:,:)=0.
  fv(:,:,:)=0.
  pv(:,:,:)=0.
  ke(:,:,:)=0.

  !     ---  Momentum Flux  ---
  do k=1,km
    fu(ib:ie-1,:,k)=0.5*ub(ib:ie-1,:,k)*(hd(ib+1:ie,:,k)+hd(ib:ie-1,:,k))
    fv(:,jb:je-1,k)=0.5*vb(:,jb:je-1,k)*(hd(:,jb+1:je,k)+hd(:,jb:je-1,k))
  end do


  ! potential vorticity
  select case (bcm)
  case (frslip)            ! free flux (free slip)
    do k=1,km
      do j=jb,je-1
        do i=ib,ie-1
!          vx = xd*(vb(i+1,j,k)-vb(i,j,k))*dv2(i,j)*dv2(i+1,j)
!          uy = yd*(ub(i,j,k)-ub(i,j+1,k))*du2(i,j)*du2(i,j+1)
!
!          !findx=(di(i,j)+di(i,j+1))*(di(i+1,j)+di(i+1,j+1))
!          !if( findx==0. ) vx = 0.
!          !findx=(di(i,j)+di(i+1,j))*(di(i,j+1)+di(i+1,j+1))
!          !if( findx==0. ) uy = 0.

          hi =  hd(i+1,j+1,k)*di(i+1,j+1)+hd(i  ,j+1,k)*di(i  ,j+1) &
               +hd(i  ,j  ,k)*di(i  ,j  )+hd(i+1,j  ,k)*di(i+1,j  )
          if( hi/=0. ) &
               hi = (di(i+1,j+1)+di(i,j+1)+di(i,j)+di(i+1,j))/hi
!          pv(i,j,k) = ( cf(j)+vx+uy )*hi
          pv(i,j,k) = ( cf(j) )*hi
        end do
      end do
    end do
  case (noslip)            ! fixed (non slip)
    do k=1,km
      do j=jb,je-1
        do i=ib,ie-1
!          vx = xd*(vb(i+1,j,k)*dv2(i+1,j)+vb(i,j,k)*dv2(i+1,j)-2.*vb(i,j,k))
!          uy = yd*(2.*ub(i,j,k)-ub(i,j,k)*du2(i,j+1)-ub(i,j+1,k)*du2(i,j+1))
!
!          !findx=di(i+1,j)+di(i+1,j+1)
!          !if( findx==0. ) vx = -2.*xd*vb(i,j,k)
!          !findx=di(i,j+1)+di(i+1,j+1)
!          !if( findx==0. ) uy =  2.*yd*ub(i,j,k)

          hi =  hd(i+1,j+1,k)*di(i+1,j+1)+hd(i  ,j+1,k)*di(i  ,j+1) &
               +hd(i  ,j  ,k)*di(i  ,j  )+hd(i+1,j  ,k)*di(i+1,j  )
          if( hi/=0. ) &
               hi = (di(i+1,j+1)+di(i,j+1)+di(i,j)+di(i+1,j))/hi
!          pv(i,j,k) = ( cf(j)+vx+uy )*hi
          pv(i,j,k) = ( cf(j) )*hi
        end do
      end do
    end do
  end select


  ! kinetic energy
!  ke(ib+1:ie,jb+1:je,:)=0.25* &
!       (ub(ib+1:ie,jb+1:je,:)**2+ub(ib:ie-1,jb+1:je,:)**2 &
!       +vb(ib+1:ie,jb+1:je,:)**2+vb(ib+1:ie,jb:je-1,:)**2 )
  ke(:,:,:)=0.

  !     boundary conditions
  select case (bcx)
  case (cyclic) ! periodic
    fu(iei+1,:,:)=fu(ibi,:,:)
    pv(iei+1,:,:)=pv(ibi,:,:)
    ke(ibi-1,:,:)=ke(iei,:,:)
  case (wall)
    ! NOTE in case of solid wall, these variables are not used
  end select
  select case (bcy)
  case (cyclic) ! periodic
    fv(:,jei+1,:)=fv(:,jbi,:)
    pv(:,jei+1,:)=pv(:,jbi,:)
    ke(:,ibi-1,:)=ke(:,jei,:)
  case (wall)
    ! NOTE in case of solid wall, these variables are not used
  end select

end subroutine diagnostic
