subroutine fqext(fft, nsm, lt, f_fft, nfq,  fqmx)
!  this subroutine is extract the frequency and amplitude from the data returned by twofft, which is complex numbers
!  and save them to file ftag, and return the first nfq dominant frequencies, to be written to file later

! 	Input
! nsm 	:	number of sampling points
! fft 	:	complex variable from twofft
! lt	: 	length of signal, used to compute frequency
! f_fft :	file tag to save the frequency and amplitude for fft1.dat
! nfq	:	number of  frequencies to keep, we only consider the number less than 10(included) 
!                       actually we take nfq = 4

!      Onput
! fqmx(nfq, 2):  the selected frequency(the first column), and amplitude(the second column) 


implicit none 
integer, parameter:: dp = kind(1.d0) 

integer, intent(in) ::  nsm, f_fft, nfq    
double complex, intent(in)  :: fft(nsm)
real(kind=dp), intent(in)   :: lt 
real(kind=dp), intent(out)  ::  fqmx(nfq,2) 

! local variables
integer :: i, j
real(kind=dp) :: f, amp, iamp, ramp 
real(kind=dp), parameter :: pi = 4.d0*datan(1.d0)
 CHARACTER(LEN=*), PARAMETER  :: fmt  = "(4f22.12)" ! the format for magnitude output


! the initial value should be zero
fqmx = 0.d0  

do i = 1, nsm/2
! radian per unit time ! without 2*pi, it is cycle per unit time
   f = (i-1) / lt  * pi * 2 ! ! for twofft, amp = amp/(nsm/2), f = f/2, but if t is the real time, f = f*2*pi/2 = f*pi

! is this the right way to get the frequency?! Yes! frequency is the number of sampling per unit time 
 
! the amplitude, the real and imaginary part, which are the coefficient of cos and sin respectively.     
  amp = cdabs( fft(i) )/ nsm * 2
  iamp = dimag(fft(i))/ nsm * 2
  ramp = dble(fft(i))/ nsm * 2 
      
  if(i == 1) then  ! for f=0, the constant shift should be devided by 2
    amp  = amp  / 2
    ramp = ramp/2  
  endif
    
  write(f_fft, fmt)  f, amp,  ramp,  iamp 
!  write(*, fmt)  f, amp,  ramp, iamp  ! frequency domain + phase domain(polar radius + phase angle)
   
  if (f .lt. 1.d-5) cycle ! the constant term, corresponding to f = 0 
  
! subroutine fqmax(fq, amp, xin, nfq, xout) 
  call fqmax(f, amp, nfq, fqmx)
  
!  print*, 'fqmx' 
!  do j = 1, nfq
!    write(*,'(3f20.10)') dx, fqmx(j, :)
!  enddo      
!  read(*,*)

!  if (dabs(amp) > 1.d-3) then  ! discard this part, we are not interested in all the 1.d-3 part...
!!    print*,  f, amp 
!    write(f_fft, fmt )  f ,amp ,  ramp, iamp 
!  endif
  
enddo 

!print*,  'finish frequency extraction of all the points'

! add a blank line to seperate orbit as block
write(f_fft, *) 
 
return
end    
    
    
