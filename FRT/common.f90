module common

  implicit none

  real(8),parameter:: &
       undef=-999.

  integer:: &
       im, ib, ie, imi, ibi, iei, &
       jm, jb, je, jmi, jbi, jei, &
       km, kmi

  real(8):: &
       dx,   &
       dy,   &
       dt,   & ! time interval
       dtmin,& ! smallest time interval
       dtmax,& ! largest time interval
       f0,   & ! Coriolis coef.
       b0,   & ! beta (gradient of Coriolis coef.)
       gr01, & ! reduced gravity between 0-1 layer
               ! (should be 1 except for the 1.5 layer model)
       gr12, & ! reduced gravity between 1-2 layer 
       gr23, & ! reduced gravity between 2-3 layer 
       hh1,  & ! layer thickness (1st layer)
       hh2,  & ! layer thickness (2nd layer)
       hh3,  & ! layer thickness (3rd layer)
       vis,  & ! eddy viscosity
       dif,  & ! eddy diffusivity
       bst     ! bottom stress

  real(8):: &
       time_to_start,   & ! start time
       time_to_end,     & ! end time
       output_interval, & ! output time interval
       monitor_interval,& ! monitor interval
       time,            & ! current time
       time_to_output,  & ! output time
       time_to_monitor    ! monitor_time


  character(10):: &
       bndtyp_x, bndtyp_y, &
       bndtyp_m, bndtyp_e

  integer:: &
       bcx, bcy, bcm, bce
  integer,parameter:: &
       cyclic=0, &
       wall=1, &
       frslip=0, &
       noslip=1, &
       insult=0, &
       condct=1

  real(8):: &
       inival_u, &
       inival_v, &
       inival_e

  character(10):: &
       inityp
  
  integer:: &
       nstep, & ! current time step
       istep, & ! internal time step
       sumup    ! counter for averaging
  
  integer,parameter:: &
       monitor_unit=10, &
       output_unit=11, &
       backup_unit=12


  ! basic arrays

  real(8),allocatable:: &
       ua(:,:,:), ub(:,:,:), gu(:,:,:), uf(:,:,:), &
       va(:,:,:), vb(:,:,:), gv(:,:,:), vf(:,:,:), &
       ea(:,:,:), eb(:,:,:), ge(:,:,:), ef(:,:,:), &
       fu(:,:,:), fv(:,:,:), pv(:,:,:), ke(:,:,:), &
       df(:,:,:), hd(:,:,:), &
       dp(:,:)  , di(:,:)  , &
       du(:,:)  , dv(:,:)  , du2(:,:) , dv2(:,:) , &
       gr(:)    , hh(:)    , cf(:)


!  real(4),allocatable:: &
!       us(:,:,:), vs(:,:,:), hs(:,:,:), ps(:,:)
!  real(4):: &
!       sumdt
  real(8),allocatable:: &
       us(:,:,:), vs(:,:,:), es(:,:,:)
  real(8):: &
       sumdt

  ! frequently or globally used constants

  real(8):: &
       xd, yd, adv, aera, depmax, dpi

contains
  subroutine configuration (expid)

    character(*),intent(in):: &
         expid
    integer:: &
         n

    namelist &
         /basic_parameter/ &
         imi, jmi, kmi, &
         dx, dy, dt, dtmin, dtmax, &
         f0, b0, vis, dif, bst, gr01, gr12, gr23, hh1, hh2, hh3, &
         time_to_start, time_to_end, output_interval, monitor_interval, &
         /inicon_parameter/ &
         inival_u, inival_v, inival_e, inityp, &
         /bndcon_parameter/ &
         bndtyp_x, bndtyp_y, bndtyp_m, bndtyp_e


    ! set default values
    im=0
    jm=0
    imi=0
    jmi=0

    dx   =undef
    dy   =undef
    dtmin=undef
    dtmax=undef
    dt   =undef

    f0=undef
    b0=undef
    gr01=9.8
    gr12=undef
    gr23=undef
    hh1=undef
    hh2=undef
    hh3=undef
    vis=undef
    dif=undef
    bst=undef

    time_to_start   =undef
    time_to_end     =undef
    output_interval =undef
    monitor_interval=undef

    inival_u=undef
    inival_v=undef
    inival_e=undef
    inityp='undef'

    bndtyp_x='undef' ! 'cyclic' or 'wall'
    bndtyp_y='undef' ! 'cyclic' or 'wall'

    bndtyp_m='undef' ! 'frslip' or 'noslip'
    bndtyp_e='undef' ! 'insulating' or 'conducting'

    call available_file_unit(n)
    open(n,file='../CNF/'//trim(expid)//'.cnf')
    read(n,basic_parameter)
    read(n,inicon_parameter)
    read(n,bndcon_parameter)
    close(n)

    ibi=1 ; iei=imi
    jbi=1 ; jei=jmi

    im=imi+2 ; ib=ibi-1 ; ie=iei+1
    jm=jmi+2 ; jb=jbi-1 ; je=jei+1
    km=kmi


    ! array allocation
    allocate( &
         ua(ib:ie,jb:je,1:km), &
         ub(ib:ie,jb:je,1:km), &
         gu(ib:ie,jb:je,1:km), &
         uf(ib:ie,jb:je,1:km), &
         va(ib:ie,jb:je,1:km), &
         vb(ib:ie,jb:je,1:km), &
         gv(ib:ie,jb:je,1:km), &
         vf(ib:ie,jb:je,1:km), &
         ea(ib:ie,jb:je,1:km), &
         eb(ib:ie,jb:je,1:km), &
         ge(ib:ie,jb:je,1:km), &
         ef(ib:ie,jb:je,1:km), &
         fu(ib:ie,jb:je,1:km), &
         fv(ib:ie,jb:je,1:km), &
         pv(ib:ie,jb:je,1:km), &
         ke(ib:ie,jb:je,1:km), &
         df(ib:ie,jb:je,1:km), &
         hd(ib:ie,jb:je,1:km) &
    )

    allocate( &
         dp(ib:ie,jb:je), &
         di(ib:ie,jb:je), &
         du(ib:ie,jb:je), &
         dv(ib:ie,jb:je), &
         du2(ib:ie,jb:je),&
         dv2(ib:ie,jb:je) &
         )
    
    allocate( &
         gr(1:km), &
         hh(1:km), &
         cf(jb:je)  &
    )

    allocate( &
         us(ib:ie,jb:je,1:km), &
         vs(ib:ie,jb:je,1:km), &
         es(ib:ie,jb:je,1:km) &
         )


    ! set default values
    select case (km)
    case (1)
      gr(1)=gr01
      hh(1)=hh1
    case(2)
      gr(1)=gr01
      gr(2)=gr12
      hh(1)=hh1
      hh(2)=hh2
    case (3)
      gr(1)=gr01
      gr(2)=gr12
      gr(3)=gr23
      hh(1)=hh1
      hh(2)=hh2
      hh(3)=hh3
    end select


    ! default ocean depth
    dp(:,:)=sum(hh(:))



    ! set boundary parameters
    select case (bndtyp_x)
    case ('cyclic')
      bcx=cyclic
    case ('wall')
      bcx=wall
    case default
      stop 'ERR:: invalid bndtyp_x'
    end select
    select case (bndtyp_y)
    case ('cyclic')
      bcy=cyclic
    case ('wall')
      bcy=wall
    case default
      stop 'ERR:: invalid bndtyp_y'
    end select
    select case (bndtyp_m)
    case ('frslip')
      bcm=frslip
    case ('noslip')
      bcm=noslip
    case default
      stop 'ERR:: invalid bndtyp_m'
    end select
    select case (bndtyp_e)
    case ('insulating')
      bce=insult
    case ('conducting')
      bce=condct
    case default
      stop 'ERR:: invalid bndtyp_e'
    end select

  end subroutine configuration


  subroutine configuration_check

    if( km >= 4 ) stop 'ERR:: km >= 4 is not expected in the present configuration.'

    if( imi==0 ) stop 'ERR:: invalid imi'
    if( jmi==0 ) stop 'ERR:: invalid imi'

    if( dx==undef ) stop 'ERR:: invalid dx'
    if( dy==undef ) stop 'ERR:: invalid dy'
    if( dtmin==undef ) stop 'ERR:: invalid dtmin'
    if( dtmax==undef ) stop 'ERR:: invalid dtmax'
    if( dtmax==undef ) stop 'ERR:: invalid dtmax'

    if( f0==undef ) stop 'ERR:: invalid f0'
    if( b0==undef ) stop 'ERR:: invalid b0'
    if( gr12==undef .and. km >=2 ) stop 'ERR:: invalid gr12'
    if( gr23==undef .and. km >=3 ) stop 'ERR:: invalid gr23'
    if( hh1==undef .and. km >=1 ) stop 'ERR:: invalid hh1'
    if( hh2==undef .and. km >=2 ) stop 'ERR:: invalid hh2'
    if( hh3==undef .and. km >=3 ) stop 'ERR:: invalid hh3'
    if( vis==undef ) stop 'ERR:: invalid vis'
    if( dif==undef ) stop 'ERR:: invalid dif'
    if( bst==undef ) stop 'ERR:: invalid bst'

    if( time_to_start   ==undef ) stop 'ERR:: invalid time_to_start'
    if( time_to_end     ==undef ) stop 'ERR:: invalid time_to_end'
    if( output_interval ==undef ) stop 'ERR:: invalid output_interval'
    if( monitor_interval==undef ) stop 'ERR:: invalid monitor_interval'

    if( inival_u==undef ) stop 'ERR:: invalid inival_u'
    if( inival_v==undef ) stop 'ERR:: invalid inival_v'
    if( inival_e==undef ) stop 'ERR:: invalid inival_e'


  end subroutine configuration_check


  subroutine useful_constants

    integer:: &
         i,j

    dpi =4.d0*datan(1.d0)

    xd=1./dx
    yd=1./dy

    aera=1./dble(imi)/dble(jmi)

    adv = -1./24.


!    select case (km)
!    case (1)
!      gr(1)=gr01
!      hh(1)=hh1
!    case(2)
!      gr(1)=gr01
!      gr(2)=gr12
!      hh(1)=hh1
!      hh(2)=hh2
!    case (3)
!      gr(1)=gr01
!      gr(2)=gr12
!      gr(3)=gr23
!      hh(1)=hh1
!      hh(2)=hh2
!      hh(3)=hh3
!    end select

    do j=jb,je
      cf(j) = b0*dble(j)*dy
    end do
    cf(:) = cf(:)-sum(cf(:))/dble(jm)
    cf(:) = cf(:)+f0

    depmax=sum(hh)

    ! default ocean depth
    !dp(:,:)=depmax
    where ( dp(:,:) > 0. ) 
      di(:,:) = 1.
    else where
      di(:,:) = 0.
    end where


    du(ib:ie-1,:)=di(ib:ie-1,:)*di(ib+1:ie,:)
    select case (bcx)
    case (cyclic)
      du(ie,:)=du(ibi,:)
    case (wall)
      du(ie,:)=0.
    end select

    dv(:,jb:je-1)=di(:,jb:je-1)*di(:,jb+1:je)
    select case (bcy)
    case (cyclic)
      dv(:,je)=dv(:,jbi)
    case (wall)
      dv(:,je)=0.
    end select

    do j=jb,je
      do i=ib,ie-1
        if( di(i,j)+di(i+1,j)==0. ) then
          du2(i,j)=0.
        else
          du2(i,j)=1.
        end if
      end do
    end do
    select case (bcx)
    case (cyclic)
      du2(ie,:)=du2(ibi,:)
    case (wall)
      du2(ie,:)=0.
    end select
    do j=jb,je-1
      do i=ib,ie
        if( di(i,j)+di(i,j+1)==0. ) then
          dv2(i,j)=0.
        else
          dv2(i,j)=1.
        end if
      end do
    end do
    select case (bcy)
    case (cyclic)
      dv2(:,je)=dv2(:,jbi)
    case (wall)
      dv2(:,je)=0.
    end select

  end subroutine useful_constants

end module common
