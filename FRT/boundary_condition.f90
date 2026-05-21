subroutine boundary_condition

  use common

  implicit none

  integer:: &
       k

  select case (bcx)
  case (cyclic )
    ub(ibi-1,:,:)= ub(iei,:,:)
    ub(iei+1,:,:)= ub(ibi,:,:)
    vb(ibi-1,:,:)= vb(iei,:,:)
    vb(iei+1,:,:)= vb(ibi,:,:)
    eb(ibi-1,:,:)= eb(iei,:,:)
    eb(iei+1,:,:)= eb(ibi,:,:)
  case (wall)
    ! eastern boundaries
    ub(iei:ie,:,:)=0.
    select case (bcm)
    case (frslip)
      vb(iei+1,:,:)= vb(iei,:,:)
    case (noslip)
      vb(iei+1,:,:)=-vb(iei,:,:)
    end select
    select case (bce)
    case (insult)
      eb(iei+1,:,:)= eb(iei,:,:)
    case (condct)
      eb(iei+1,:,:)=-eb(iei,:,:)
    end select

    ! western boundaries
    ub(ib:ibi-1,:,:)= 0.
    select case (bcm)
    case (frslip)
      vb(ibi-1,:,:)= vb(ibi,:,:)
    case (noslip)
      vb(ibi-1,:,:)=-vb(ibi,:,:)
    end select
    select case (bce)
    case (insult)
      eb(ibi-1,:,:)= eb(ibi,:,:)
    case (condct)
      eb(ibi-1,:,:)=-eb(ibi,:,:)
    end select
  end select


  select case (bcy)
  case (cyclic)
    ub(:,jbi-1,:)= ub(:,jei,:)
    ub(:,jei+1,:)= ub(:,jbi,:)
    vb(:,jbi-1,:)= vb(:,jei,:)
    vb(:,jei+1,:)= vb(:,jbi,:)
    eb(:,jbi-1,:)= eb(:,jei,:)
    eb(:,jei+1,:)= eb(:,jbi,:)
  case (wall)
    ! northern boundaries
    vb(:,jei:je,:)= 0.
    select case (bcm)
    case (frslip)
      ub(:,jei+1,:)= ub(:,jei,:)
    case (noslip)
      ub(:,jei+1,:)=-ub(:,jei,:)
    end select
    select case (bce)
    case (insult)
      eb(:,jei+1,:)= eb(:,jei,:)
    case (condct)
      eb(:,jei+1,:)=-eb(:,jei,:)
    end select


    ! southern boundaries
    vb(:,jb:jbi-1,:)= 0.
    select case (bcm)
    case (frslip)
      ub(:,jbi-1,:)= ub(:,jbi,:)
    case (noslip)
      ub(:,jbi-1,:)=-ub(:,jbi,:)
    end select
    select case (bce)
    case (insult)
      eb(:,jbi-1,:)= eb(:,jbi,:)
    case (condct)
      eb(:,jbi-1,:)=-eb(:,jbi,:)
    end select
  end select

end subroutine boundary_condition
