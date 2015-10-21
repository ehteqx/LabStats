! RANDIST - v. 1.2 ENH {1}  ||  STATOOLS - v. 2.1 [Fortran Edition] {2}
! {1} A pseudo-ramdomic generator of flat, linear and exponential real-number distributions.
! {2} An essential collection of tools for statistical data analysis written in Fortran 95.

! (C) Emanuele Ballarin (ehteqx@gmail.com) -- 20/06/2015
!###############################################################################################

MODULE GLPREC			! Module for global precision constants
	
	implicit none
	
	integer, parameter 	:: ik = selected_int_kind(38)		! MAX: 38
	integer, parameter 	:: rk = selected_real_kind(33)		! MAX: 33
	
END MODULE GLPREC

MODULE STATOOLS			! Module for STATOOLS (Tools for statistical data analysis)
	
	use GLPREC
	implicit none
	
	CONTAINS
	
	subroutine binrange(dataarray, n, division, histmin, histmax, histep, bincount)		! Subroutine for adaptive bin-range definition
	
		implicit none
		
		real (kind = rk), dimension(n), intent(in) :: dataarray
		integer (kind = ik), intent(in) :: division
		real (kind = rk) :: minvalue, maxvalue, step
		integer (kind = ik), intent(in) :: n
		integer (kind = ik), intent(out) :: bincount
		real (kind = rk), intent(out) :: histmin, histmax, histep
		
		minvalue = 0.0_rk  				! MINVAL(dataarray) # UNCOMMENT FOR ADAPTIVE BEHAVIOUR (lower limit)
		maxvalue = MAXVAL(dataarray)	! The adaptive behaviour for the upper limit is default
		
		step = ((maxvalue - minvalue)/division) + (MINVAL(dataarray)/(2*division)) 	! Action required to avoid unfilled leftovers
		
		histmin = (minvalue - step)
		histmax = (maxvalue + step)
		histep = step
		bincount = division + 1_ik 			! The last bin MUST BE LEFT EMPTY. This verifies the presence of eventual leftovers.
			
	end subroutine binrange
	
	
	subroutine binfill(dataarray, n, histmin, histmax, histep, bincount, histogram)		! Subroutine for adaptive bin filling
	
		implicit none
		
		real(kind=rk), dimension(n), intent(in) :: dataarray
		real (kind = rk), intent(in) :: histmin, histmax, histep
		integer (kind = ik), intent(in) :: n, bincount
		integer (kind = ik), allocatable, intent(out) :: histogram(:)
		integer (kind = ik) :: i, j, upbound, dowbound
		
						
		allocate(histogram(bincount))
	
		do i=1_ik, bincount	
			histogram(i) = 0_rk
		end do
	
		do i = 1_ik, n
		
			checker: do j = 1_ik, bincount
			
				if (dataarray(i) >= (j-1_ik)*histep .AND. dataarray(i) < (j)*histep) then
					histogram(j) = histogram(j) + 1_ik
					
					exit checker
									
				end if
			
			end do checker
		
		end do
			
	end subroutine binfill
	
END MODULE STATOOLS


PROGRAM RANDIST
	
	use GLPREC
	use STATOOLS
	implicit none
	
	integer (kind = ik) :: i, j, n, division
	integer (kind = 4) :: rseed=4395432
	real (kind = rk) :: tau
	real (kind = rk) ::	histmin_flat, histmin_linear, histmin_expon
	real (kind = rk) ::	histmax_flat, histmax_linear, histmax_expon
	real (kind = rk) ::	histep_flat, histep_linear, histep_expon
	integer (kind = ik) ::	bincount_flat, bincount_linear, bincount_expon
	real (kind = rk), allocatable :: flat(:), linear(:), expon(:)
	integer (kind = ik), allocatable :: histogram_flat(:), histogram_linear(:), histogram_expon(:)
		
	
	tau = 0.5_rk	! MODIFY AS NEEDED
	
	print*, 'Inserisci il numero di valori da generare... '
	read*, n
	
	print*, 'Inserisci il numero di bin dello istogramma ... '
	read*, division
	
	allocate(flat(n))
	allocate(linear(n))
	allocate(expon(n))
	
	CALL RANDOM_SEED(rseed)
	CALL RANDOM_NUMBER(flat)
	
	do i = 1_ik, n
	
		linear(i) = sqrt(flat(i))
		expon(i) = (-1.0_rk)*tau*log(1.0_rk-flat(i))
		
	end do
	
	! Invoking the statistical functions and copying the data
	
	CALL binrange(flat, n, division, histmin_flat, histmax_flat, histep_flat, bincount_flat)
	CALL binrange(linear, n, division, histmin_linear, histmax_linear, histep_linear, bincount_linear)
	CALL binrange(expon, n, division, histmin_expon, histmax_expon, histep_expon, bincount_expon)
	
	allocate(histogram_flat(bincount_flat))
	allocate(histogram_linear(bincount_linear))
	allocate(histogram_expon(bincount_expon))
	
	CALL binfill(flat, n, histmin_flat, histmax_flat, histep_flat, bincount_flat, histogram_flat)
	CALL binfill(linear, n, histmin_linear, histmax_linear, histep_linear, bincount_linear, histogram_linear)
	CALL binfill(expon, n, histmin_expon, histmax_expon, histep_expon, bincount_expon, histogram_expon)	
	
	! Writing to file all the data...
	
	open(unit=1, file='distributions.dat')
	open(unit=2, file='h-flat.dat')
	open(unit=3, file='h-linear.dat')
	open(unit=4, file='h-exponential.dat')

	do i=1_ik, n	
		write(unit=1,fmt=*)i, flat(i), expon(i), linear(i)		
	end do
	
	do j = 1_ik, bincount_flat
		write(unit=2,fmt=*)j, histogram_flat(j)		
	end do
	
	do j = 1_ik, bincount_linear
		write(unit=3,fmt=*)j, histogram_linear(j)		
	end do
	
	do j = 1_ik, bincount_expon
		write(unit=4,fmt=*)j, histogram_expon(j)		
	end do
	
END PROGRAM RANDIST
