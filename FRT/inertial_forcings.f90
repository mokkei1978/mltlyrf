subroutine inertial_forcings

  use common

  implicit none

  real(8):: &
       alph, &
       beta, &
       gamm, &
       delt, &
       phiu, &
       phil, &
       epsu, &
       epsl, &
       g1,   &
       g2

  integer:: &
       i, j, k

  gu(:,:,:)=0.
  gv(:,:,:)=0.
  ge(:,:,:)=0.


  ! pressure gradient
  ! surface pressure
  gu(ib:ie-1,:,1)=gr(1)*xd*(eb(ib:ie-1,:,1)-eb(ib+1:ie,:,1))
  gv(:,jb:je-1,1)=gr(1)*yd*(eb(:,jb:je-1,1)-eb(:,jb+1:je,1))
  ! internal pressure
  do k=2,km
     gu(ib:ie-1,:,k)=gu(ib:ie-1,:,k-1)+gr(k)*xd*(eb(ib:ie-1,:,k)-eb(ib+1:ie,:,k))
     gv(:,jb:je-1,k)=gv(:,jb:je-1,k-1)+gr(k)*yd*(eb(:,jb:je-1,k)-eb(:,jb+1:je,k))
  end do

  ! gu
  do k=1,km
    do j=jbi,jei
      do i=ibi,iei
        alph = 2.*pv(i+1,j,k)+pv(i,j,k)+2.*pv(i,j-1,k)+pv(i+1,j-1,k)
        beta = pv(i,j,k)+2.*pv(i-1,j,k)+pv(i-1,j-1,k)+2.*pv(i,j-1,k)
        gamm = 2.*pv(i,j,k)+pv(i-1,j,k)+2.*pv(i-1,j-1,k)+pv(i,j-1,k)
        delt = pv(i+1,j,k)+2.*pv(i,j,k)+pv(i,j-1,k)+2.*pv(i+1,j-1,k)
        epsu = pv(i+1,j,k)+pv(i,j,k)-pv(i,j-1,k)-pv(i+1,j-1,k)
        epsl = pv(i,j,k)+pv(i-1,j,k)-pv(i-1,j-1,k)-pv(i,j-1,k)

        g1 = -alph*fv(i+1,j  ,k)-beta*fv(i  ,j  ,k) &
             -gamm*fv(i  ,j-1,k)-delt*fv(i+1,j-1,k) &
             +epsu*fu(i+1,j  ,k)-epsl*fu(i-1,j  ,k)

        g2 = ke(i,j,k)-ke(i+1,j,k) ! energy pressure

        gu(i,j,k)=gu(i,j,k) & ! pressure term
             +adv*g1+xd*g2 ! nonlinear term
      end do
    end do
  end do



  ! gv
  do k=1,km
    do j=jbi,jei
      do i=ibi,iei
        alph = 2.*pv(i,j,k)+pv(i-1,j,k)+2.*pv(i-1,j-1,k)+pv(i,j-1,k)
        beta = pv(i,j,k)+2.*pv(i-1,j,k)+pv(i-1,j-1,k)+2.*pv(i,j-1,k)
        gamm = 2.*pv(i,j+1,k)+pv(i-1,j+1,k)+2.*pv(i-1,j,k)+pv(i,j,k)
        delt = pv(i,j+1,k)+2.*pv(i-1,j+1,k)+pv(i-1,j,k)+2.*pv(i,j,k)
        phiu = -pv(i,j+1,k)+pv(i-1,j+1,k)+pv(i-1,j,k)-pv(i,j,k)
        phil = -pv(i,j,k)+pv(i-1,j,k)+pv(i-1,j-1,k)-pv(i,j-1,k)

        g1 = +gamm*fu(i  ,j+1,k)+delt*fu(i-1,j+1,k) &
             +alph*fu(i-1,j  ,k)+beta*fu(i  ,j  ,k) &
             +phiu*fv(i  ,j+1,k)-phil*fv(i  ,j-1,k)

        g2 = ke(i,j,k)-ke(i,j+1,k)

        gv(i,j,k)=gv(i,j,k) & ! pressure term
             +adv*g1+yd*g2 ! nonlinear term
      end do
    end do
  end do


  ! ge
  do k=1,km
    do j=jbi,jei
      do i=ibi,iei
        g1 = fu(i-1,j,k)-fu(i,j,k)
        g2 = fv(i,j-1,k)-fv(i,j,k)

        ge(i,j,k)=ge(i,j,k) &
             +xd*g1+yd*g2
      end do
    end do
  end do
  do k=km-1,1,-1
     ge(:,:,k)=ge(:,:,k)+ge(:,:,k+1)
  end do

 end subroutine inertial_forcings
