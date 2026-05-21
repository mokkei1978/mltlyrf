subroutine available_file_unit(n)

  implicit none

  integer:: &
       m,n
  logical:: &
       file_status

  n=0
  do m=7,99
     inquire( m,opened=file_status )
     if( file_status ) cycle
     n=m
     exit
  end do

  if( n==0 ) then
     write(6,*) 'ERR:: no file unit is available'
     stop
  end if

end subroutine available_file_unit
