! EXAM EXERCISE N° 5.2 -- SIMULATION -- Version 2.0
! A straightforward, automated simulator for large scale iterated measurements
! (C) Emanuele Ballarin (ehteqx@gmail.com)
!##############################################################################

! § MODULES (Declaration and Definition)

MODULE GLPREC			! Module for global precision constants

	implicit none

	integer, parameter 	:: ik = selected_int_kind(38)		! MAX: 38
	integer, parameter 	:: rk = selected_real_kind(33)		! MAX: 33

    END MODULE GLPREC

MODULE GLCONST			! Module for global and universal constants

        use GLPREC
        implicit none

    	real (kind = rk), parameter 	:: twopi = 2*acos(-1.0_rk)
        real (kind = rk), parameter 	:: pi = acos(-1.0_rk)

    END MODULE GLCONST

MODULE GAUSSIAN         ! Module for the general-purpose generation of Gaussian Distributions

    use GLPREC
    use GLCONST
    implicit none

    CONTAINS

        SUBROUTINE POLARGAUSS(output, n, mean, stdev, rseed_a, rseed_b, rseed_c)	! Subroutine for the generation of Gausssian Distributions (Marsaglia Polar Method)

            implicit none

                integer (kind = ik), intent(in) :: n                         ! Dimension of the output array (must be user-passed)
                real (kind = rk), dimension(n), intent(out) :: output        ! Output array (contains the distribution)
                real (kind = rk), intent(in) :: mean, stdev                  ! Mean and std deviation of the distribution
                integer (kind = 4), intent(inout) :: rseed_a, rseed_b, rseed_c  ! The random seeds used with the subroutine (change for different results)

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

MODULE LSLMFIT         ! Module for the general-purpose calculation of best linear fits and associated uncertainties.

    use GLPREC
    implicit none

    CONTAINS

        SUBROUTINE LINFIT(dimarr, arr_x, arr_y, errvar, out_m, out_q, out_cov, out_corr, out_eq, out_em, linefit)		! Subroutine for the calculation of linear fits (LSLM)

            implicit none

				integer (kind = ik), intent(in) :: dimarr
				real (kind = rk), dimension(dimarr), intent(in) :: arr_x, arr_y
				real (kind = rk), dimension(dimarr), intent(inout) :: errvar
				real (kind = rk), dimension(dimarr), intent(out) :: linefit
				real (kind = rk), intent(out) :: out_m, out_q, out_cov, out_corr, out_eq, out_em

				real (kind = rk) :: s0, s1, s2, s3, s4, s5, d
				integer (kind = ik) :: i, k

				! The following variables are initialized with a default value
				s0 = 0.0_rk
				s1 = 0.0_rk
				s2 = 0.0_rk
				s3 = 0.0_rk
				s4 = 0.0_rk

                ! Formula per LSLM (StdDev o StdDev * Values) FLAGGED !
                do k=1,dimarr
                    errvar(k) = (((errvar(k))*(arr_y(k))/(sqrt(3.0_rk)))**(2.0_rk))
                end do

				! The sums are computed within a do-cycle
				do i=1,dimarr
					s0 = s0+1_rk/errvar(i)
					s1 = s1+arr_y(i)/errvar(i)
					s2 = s2+arr_x(i)/errvar(i)
					s3 = s3+arr_x(i)*arr_y(i)/errvar(i)
					s4 = s4+arr_x(i)*arr_x(i)/errvar(i)
				end do

				! Some additional temporary values and the final ones are computed
				d = s0*s4-s2**2
				out_m = (s0*s3-s2*s1)/d
				out_q = (s4*s1-s2*s3)/d
				out_em = s0/d
				out_eq = s4/d
				out_cov = -s2/d
				out_corr = out_cov/sqrt(out_em*out_eq)

				! The values of the fit function (linear) are computed and written to the array
				do i=1,dimarr
					linefit(i)=out_m*arr_x(i)+out_q
				end do

            END SUBROUTINE

    END MODULE LSLMFIT

! ******************************************************************************

! § MAIN PROGRAM (Nonmodular, Interactive, Complete)

PROGRAM SIMULATION

	use GLPREC
	use GLCONST
	use GAUSSIAN
	use LSLMFIT
	implicit none

	! Indexes
	integer (kind = ik) :: z = 0, w = 0, memory = 0
	integer (kind = ik) :: k = 0

	! Values of m_0 and q_0 in [ y = (m_0)*x + (q_0) ]; Arbitrary, but fixed.
	real (kind = rk) :: m_zero = pi, q_zero = 1.35

	! Arrays of x-es and y-es
	real (kind = rk), dimension(7) :: arr_x = [1.0,2.0,3.0,4.0,5.0,6.0,7.0] ! Fixed sampling of x-es
	real (kind = rk), dimension(7) :: arr_y = [0.0,0.0,0.0,0.0,0.0,0.0,0.0] ! Initialized zeroes
	real (kind = rk), dimension(7) :: error_y = [0.0,0.0,0.0,0.0,0.0,0.0,0.0] ! Initialized zeroes

	! Things necessary to generate the random seeds for gaussian error generation
	integer (kind = ik) :: cntsim = 0 		             ! Number of simulations (user-defined)
	integer (kind = 4) :: foresee_d				             ! Seed of seeds (user-defined for randomness)
	real (kind = rk), allocatable :: aseeds_buffer = 0	     ! Array of seeds, buffered, real-valued
	integer (kind = 4), allocatable :: aseeds = 0            ! Array of seeds, positive integers

	! Result arrays (estimated m_0 and q_0 for each simulation)
	real (kind = rk), allocatable :: estimated_m = 0.0
	real (kind = rk), allocatable :: estimated_q = 0.0

	! Mock variables (end nowhere)
	real (kind = rk) :: out_m, out_q, out_cov, out_corr
	real (kind = rk), dimension(7) :: pit

    ! Bootstrapping message
    print*, ' '
    print*, '#######################################################################'
	print*, '                       SIMFIT - v. 2.1 GAUSSIAN                       '
	print*, '                                                                      '
	print*, '                  Copyright (C) 2015 Emanuele Ballarin                '
	print*, '                                                                      '
	print*, 'This is free software, covered by the GNU General Public License v3,  '
	print*, 'and comes with ABSOLUTELY NO WARRANTY WHATSOEVER.                     '
	print*, 'For more information about the license: https://www.gnu.org/licenses/ '
	print*, '#######################################################################'
    print*, ' '
    print*, 'Program loaded!'

    ! Some user interaction (for the random seed and the number of simulations)
    print*, 'Insert the number of simulations that you want to run... '
    read*, cntsim
    print*, ' '

    print*, 'Inserisci un numero di 7 (sette) cifre da usare come seme random... '
    read*, foresee_d
    print*, ' '

    ! User notification
    print*, 'Pre-computing random numbers for the simulation...'

    ! The seeds for gaussian error generation are finally generated
	allocate(aseeds_buffer(3_ik*cntsim))
	allocate(aseeds(3_ik*cntsim))

	CALL RANDOM_SEED(foresee_d)
	CALL RANDOM_NUMBER(aseeds_buffer)

    aseeds_buffer = ((aseeds_buffer)*(10000000.0_rk))

    do z = 1_ik,cntsim,1_ik
        aseeds(z) = INT((floor(aseeds_buffer(z))), ik)
    end do

    ! User notification
    print*, 'OK.'
    print*, ' '

    print*, 'Simulation started. Please wait, it may take some time to complete...'
    print*, ' '
    print*, 'Using', m_zero, 'as m_0 central value'
    print*, 'Using', q_zero, 'as q_0 central value'
    print*, ' '

    ! The simulation cycle is started and the necessary data is computed
    allocate(estimated_m(cntsim))
    allocate(estimated_q(cntsim))

    do w = 1, cntsim, 1 ! Cycle of simulations

        arr_y = (((m_zero)*(arr_x)) + (q_zero)) ! The y-es array is (re)initialized

        ! Generating errors...
        CALL POLARGAUSS(error_y, 7_ik, 0.0_rk, 0.05, aseeds(memory + 1_ik), aseeds(memory + 2_ik), aseeds(memory + 3_ik))

        arr_y = (arr_y + error_y) ! Adding errors to the y-es array

        CALL LINFIT(7_ik, arr_x, arr_y, ((0.05_rk)**2_ik), out_m, out_q, out_cov, out_corr, estimated_q(w), estimated_m(w), pit)

        memory = memory + 3_ik ! Info for the random numbers generator

    end do

    ! User notification
    print*, 'OK.'
    print*, ' '
    print*, ' '
    print*, 'Bye bye, have a nice day!'
    print*, ' '


END PROGRAM SIMULATION
