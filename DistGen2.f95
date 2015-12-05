! DISTGEN2 - v. 1.2 NGD
! A pseudo-ramdomic generator of real-number distributions, derived from the Standard Gaussian Distribution

! It relies heavily on George Marsaglia's KISS linear RNG (Fortran built-in) and G.M.'s Box-Muller Polar Method (reimplemented).
! It uses the GLPREC, HISTOG and POLARGAUSS modules/subroutines by Emanuele Ballarin.

! (C) Emanuele Ballarin (ehteqx@gmail.com)
!###############################################################################################################################

! ### MODULES IMPORT & DEFINITION ###

MODULE GLPREC			                 ! Module for global precision constants

	implicit none

	integer, parameter 	:: ik = selected_int_kind(38)		! MAX: 38
	integer, parameter 	:: rk = selected_real_kind(33)		! MAX: 33

END MODULE GLPREC

MODULE GLCONST			! Module for global and universal constants

        use GLPREC
        implicit none

    	real (kind = rk), parameter 	:: twopi = 2*acos(-1.0_rk)
    	!integer (kind = ik), parameter 	::

    END MODULE GLCONST

MODULE HISTOOLS			        ! Module for histogram manipulation and plotting

	use GLPREC
	implicit none

	CONTAINS

	subroutine binrange(dataarray, n, division, histmin, histmax, bincount, binwidth) ! Subroutine for adaptive bin setup

		implicit none

		real (kind = rk), dimension(n), intent(in) :: dataarray       ! Data input
        integer (kind = ik), intent(in) :: n                          ! Dimension of the array
		real (kind = rk), intent(in) :: division                      ! Bin width (nearest to...)

		real (kind = rk) :: histspan                                  ! Working variable

		integer (kind = ik), intent(out) :: bincount                  ! Number of bins
		real (kind = rk), intent(out) :: histmin, histmax, binwidth   ! Some output variables

        histmin = real((floor(MINVAL(dataarray))), rk)              !  Histogram minimum (NOT OPTIMIZED - MODIFY FOR AD-HOC BEHAVIOUR)
        histmax = real((ceiling(MAXVAL(dataarray)) + 1), rk)        !  Histogram maximum (NOT OPTIMIZED - MODIFY FOR AD-HOC BEHAVIOUR)
        histspan = histmax - histmin                                ! Histogram spanning space
        bincount = anint((histspan/division), ik)                   ! Number of bins
        binwidth = histspan/real((bincount), rk)                    ! Bin amplitude

	end subroutine binrange


	subroutine binfill(dataarray, n, histmin, histmax, bincount, binwidth, histogram)		! Subroutine for bin filling

		implicit none

		real(kind=rk), dimension(n), intent(in) :: dataarray                      ! Data input
        integer (kind = ik), intent(in) :: n, bincount                            ! Array dimension // Number of bins
		real (kind = rk), intent(in) :: histmin, histmax, binwidth                ! As above (binrange)
		integer (kind = ik), allocatable, intent(out) :: histogram(:)             ! Output array
		integer (kind = ik) :: i, j                                               ! Indexes
        real (kind = ik) :: upbound, downbound                                     ! Working variables

		allocate(histogram(bincount))

		do i=1_ik, bincount
			histogram(i) = 0_ik
		end do


        do i = 1_ik, n

            downbound = histmin
            upbound = histmin + binwidth

			filler: do j = 1_ik, bincount

				if ((dataarray(i) >= downbound) .AND. (dataarray(i) < upbound)) then
					histogram(j) = histogram(j) + 1_ik

					exit filler

                else
                    downbound = downbound + binwidth
                    upbound = upbound + binwidth

				end if

			end do filler

		end do

	end subroutine binfill

END MODULE HISTOOLS

MODULE GAUSSIAN         ! Module for the general-purpose generation of Gaussian Distributions

    use GLPREC
    use GLCONST
    implicit none

    CONTAINS

        SUBROUTINE POLARGAUSS(output, n, mean, stdev)		! Subroutine for the generation of Gausssian Distributions (Marsaglia Polar Method)

            implicit none

                real (kind = rk), dimension(n), intent(out) :: output       ! Output array (contains the distribution)
                integer (kind = ik), intent(in) :: n                        ! Dimension of the output array (must be user-passed)
                real (kind = rk), intent(in) :: mean, stdev                ! Mean and std deviation of the distribution

                integer (kind = 4) :: rseed_a = 5931041, rseed_b = 6327490, rseed_c = 6301731  ! The random seeds used with the subroutine (change for different results)
                real (kind = ik), dimension(n) :: randomarr_a, randomarr_b, randomarr_c        ! Random-filled arrays
                integer (kind = ik) :: i, memo = 0_ik

                CALL RANDOM_SEED(rseed_a)
                CALL RANDOM_NUMBER(randomarr_a)

                CALL RANDOM_SEED(rseed_b)
                CALL RANDOM_NUMBER(randomarr_b)

                CALL RANDOM_SEED(rseed_c)
                CALL RANDOM_NUMBER(randomarr_c)

                do i = 1, n, 1

                    do while (randomarr_a(i) .EQ. 0.0_rk)

                        randomarr_a(i) = randomarr_c(memo+1_ik)
                        memo = memo + 1

                        end do

                    output(i) = ((sqrt(-2.0_rk * log(randomarr_a(i))) * cos(twopi * randomarr_b(i))) * stdev) + mean

                    end do

            END SUBROUTINE

    END MODULE GAUSSIAN

! ### MAIN PROGRAM ###

PROGRAM RANDIST

	use GLPREC
	use GLCONST
	use HISTOOLS
	use GAUSSIAN
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

	! # SOME COSMETICS #

	print*, ' '
	print*, '#######################################################################'
	print*, '               RANDIST - v. 1.2 ENH || STATOOLS - v. 2.1              '
	print*, '                                                                      '
	print*, '                  Copyright (C) 2015 Emanuele Ballarin                '
	print*, '                                                                      '
	print*, 'All this is free software, covered by the GNU General Public License v3,'
	print*, 'and comes with ABSOLUTELY NO WARRANTY WHATSOEVER.                     '
	print*, 'For more information about the license: https://www.gnu.org/licenses/ '
	print*, '#######################################################################'
	print*, ' '

	print*, 'Inserisci il numero di valori da generare... '
	read*, n

	print*, 'Inserisci il numero di bin dello istogramma ... '
	read*, division

	print*, 'Computation started...'		! Some info for the user

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

	print*, 'Computation completed!'		! Some info for the user

END PROGRAM RANDIST

! ### END OF THE FILE
