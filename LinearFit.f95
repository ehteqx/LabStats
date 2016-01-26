! LinearFit - v. 1.1
! A simple CLI software for computing the best linear fit (and associated uncertainties) of a data sample using the Least Squares method.
! (C) Emanuele Ballarin (ehteqx@gmail.com) -- 20/10/2015
!#######################################################################################################################################

program LinearFit

	implicit none

	integer, parameter 	:: ik = selected_int_kind(38)		! MAX: 38
	integer, parameter 	:: rk = selected_real_kind(33)		! MAX: 33

	real (kind = rk) :: s0, s1, s2, s3, s4, s5, m, q, d, sm, sq, cov, ro
	integer (kind = ik) :: arraydim = 7   ! MODIFY AS NEEDED THE DIMENSION OF THE FOLLOWING ARRAY
	real (kind = rk), dimension(7) :: x = [1.0,2.0,3.0,4.0,5.0,6.0,7.0], y= [3.2,5.3,6.5,9.0,10.3,13.4,14.5], vary, retta ! MODIFY AS NEEDED
	integer (kind = ik) :: i

	! # SOME COSMETICS #

	print*, ' '
	print*, '#######################################################################'
	print*, '                         LinearFit - v. 1.1                           '
	print*, '                                                                      '
	print*, '                 Copyright (C) 2015 Emanuele Ballarin                 '
	print*, '                                                                      '
	print*, 'LinearFit is free software, covered by the GNU General Public License v3,'
	print*, 'and comes with ABSOLUTELY NO WARRANTY WHATSOEVER.                     '
	print*, 'For more information about the license: https://www.gnu.org/licenses/ '
	print*, '#######################################################################'
	print*, ' '

	! The following variables are initialized with a default value
	s0 = 0.0_rk
	s1 = 0.0_rk
	s2 = 0.0_rk
	s3 = 0.0_rk
	s4 = 0.0_rk

	vary = (0.05_rk*y)**2 		! MODIFY THE VALUE OF VARIANCE AS NEEDED

	! The sums are computed within a do-cycle
	do i=1,arraydim
		s0 = s0+1_rk/vary(i)
		s1 = s1+y(i)/vary(i)
		s2 = s2+x(i)/vary(i)
		s3 = s3+x(i)*y(i)/vary(i)
		s4 = s4+x(i)*x(i)/vary(i)
	end do

	! Some additional temporary values and the final ones are computed
	d = s0*s4-s2**2
	m = (s0*s3-s2*s1)/d
	q = (s4*s1-s2*s3)/d
	sm = s0/d
	sq = s4/d
	cov = -s2/d
	ro = cov/sqrt(sm*sq)

	! The required values are then printed
	print*,'Coefficiente angolare retta di fit (m): ',m
	print*,'Intercetta retta di fit (q): ',q
	print*,'Deviazione standard su m: ',sm
	print*,'Deviazione standard su q: ',sq
	print*,'Covarianza (di m, q): ',cov
	print*,'Coefficiente Correlazione (tra m, q)',ro

	! The values of the fit function (linear) are computed and written to the array
	do i=1,arraydim
		retta(i)=m*x(i)+q
	end do

	! The values are written to a file, in a structured form
	do i=1,arraydim
		write(unit=1,fmt=*)x(i),y(i),retta(i),vary(i)
	end do

end program LinearFit
