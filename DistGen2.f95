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

PROGRAM DISTGEN2

	use GLPREC				! All the required modules are imported
	use GLCONST
	use HISTOOLS
	use GAUSSIAN
	implicit none			! As it is standard-prescribed

	! ## VARIABLES DEFINITION ##

	integer (kind = ik) :: dim		! Dimension of the default distribution array
	real (kind = rk), allocatable :: stdnorm(:), square(:), antilin(:), antisq(:), antisqsh(:)		! Distribution arrays
	integer (kind = ik), allocatable :: h_stdnorm(:), h_square(:), h_antilin(:), h_antisq(:), h_antisqsh(:)	! Histogram counters
	real (kind = rk) :: division

	real (kind = rk) :: h_min, h_max, binwth		! Working variables for subroutines
	integer (kind = ik) :: bincnt

	integer (kind = 4) :: rseed = 6423389		! For the random number generator (change for different results)

	integer (kind = ik) :: i, j 		! Indexes (some can be unnecessary)

	! ## INTRO SCREEN ##

	print*, ' '
	print*, '#######################################################################'
	print*, '                         DISTGEN2 - v. 1.2 NGD                        '
	print*, '                                                                      '
	print*, '                  Copyright (C) 2015 Emanuele Ballarin                '
	print*, '                                                                      '
	print*, 'This is free software, covered by the GNU General Public License v3,  '
	print*, 'and comes with ABSOLUTELY NO WARRANTY WHATSOEVER.                     '
	print*, 'For more information about the license: https://www.gnu.org/licenses/ '
	print*, '#######################################################################'
	print*, ' '

	! ## ... ##

	print*, 'Inserisci il numero di valori da generare... '
	read*, dim

	print*, 'Inserisci la ampiezza desiderata per i bin dello istogramma ... '
	read*, division				! User-proposed bin width (may be adapted by the program)

	print*, 'Computation started...'		! Some info for the user

	allocate(stdnorm(dim))
	allocate(square(dim))
	allocate(antilin(dim))
	allocate(antisq(dim))
	allocate(antisqsh(dim))

	! Generating the Standard Gaussian Distribution

	CALL POLARGAUSS(stdnorm, dim, 0.0_rk, 1.0_rk)

	! Generating all the other distributions

	do i = 1, dim, 1

		square(i) = ((stdnorm(i))*(stdnorm(i)))
		antilin(i) = ((1.0_rk)/(stdnorm(i)))
		antisq(i) = ((1.0_rk)/((stdnorm(i))*(stdnorm(i))))
		antisqsh(i) = ((1.0_rk)/((stdnorm(i) + 10.0_rk)*(stdnorm(i) + 10.0_rk)))

	end do

	! Invoking the plotting functions and copying the data to file

	CALL binrange(stdnorm, dim, division, h_min, h_max, bincnt, binwth)
	allocate(h_stdnorm(bincnt))
	CALL binfill(stdnorm, dim, h_min, h_max, bincnt, binwth, h_stdnorm)

	open(unit=2, file='stdnorm.dat')
	write(unit=2,fmt=*)"bin number ", " times", " left bound "
	do j = 1_ik, bincnt
		write(unit=2,fmt=*)j, h_stdnorm(j), (h_min + (j*(binwth)))
	end do

	CALL binrange(square, dim, division, h_min, h_max, bincnt, binwth)
	allocate(h_square(bincnt))
	CALL binfill(square, dim, h_min, h_max, bincnt, binwth, h_square)

	open(unit=3, file='square.dat')
	write(unit=3,fmt=*)"bin number ", " times", " left bound "
	do j = 1_ik, bincnt
		write(unit=3,fmt=*)j, h_square(j), (h_min + (j*(binwth)))
	end do

	CALL binrange(antilin, dim, division, h_min, h_max, bincnt, binwth)
	allocate(h_antilin(bincnt))
	CALL binfill(antilin, dim, h_min, h_max, bincnt, binwth, h_antilin)

	open(unit=4, file='antilin.dat')
	write(unit=4,fmt=*)"bin number ", " times", " left bound "
	do j = 1_ik, bincnt
		write(unit=4,fmt=*)j, h_antilin(j), (h_min + (j*(binwth)))
	end do

	CALL binrange(antisq, dim, division, h_min, h_max, bincnt, binwth)
	allocate(h_antisq(bincnt))
	CALL binfill(antisq, dim, h_min, h_max, bincnt, binwth, h_antisq)

	open(unit=5, file='antisq.dat')
	write(unit=5,fmt=*)"bin number ", " times", " left bound "
	do j = 1_ik, bincnt
		write(unit=5,fmt=*)j, h_antisq(j), (h_min + (j*(binwth)))
	end do

	CALL binrange(antisqsh, dim, division, h_min, h_max, bincnt, binwth)
	allocate(h_antisqsh(bincnt))
	CALL binfill(antisqsh, dim, h_min, h_max, bincnt, binwth, h_antisqsh)

	open(unit=6, file='antisqsh.dat')
	write(unit=6,fmt=*)"bin number ", " times", " left bound "
	do j = 1_ik, bincnt
		write(unit=6,fmt=*)j, h_antisqsh(j), (h_min + (j*(binwth)))
	end do

	! Writing to file all the remaining data...

	open(unit=1, file='distributions.dat')

	write(unit=1,fmt=*)"index ", " stdnorm ", " square ", " antilin ", " antisq ", " antisqsh "
	do j=1_ik, dim
		write(unit=1,fmt=*)j, stdnorm(j), square(j), antilin(j), antisq(j), antisqsh(j)
	end do

	print*, 'Computation completed!'		! Some info for the user

END PROGRAM DISTGEN2

! ### END OF THE FILE
