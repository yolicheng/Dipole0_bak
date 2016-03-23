      subroutine prt_eigval( n, ftag, wr, wi )
! print the eigenvalue 
      integer          n, ftag
      double precision wr( * ), wi( * )
!
      double precision zero
      parameter        ( zero = 0.0 )
      integer          j
!
!     write(ftag,*)
      do j = 1, n
         if( wi( j ).eq.zero ) then
            write(ftag,9998,advance='no') wr(j)
         else
            write(ftag,9999,advance='no') wr( j ), wi( j )
         end if
      end do
      write(ftag,*)
!
! 9998 format( 11(:,1x,f12.8) )
! 9999 format( 11(:,1x,'(',f14.8,',',f14.8,')') )
 9998 format( 11(f14.8, 1x) )
 9999 format( 11(f14.8, 2x, f14.8) )
      return
      end

