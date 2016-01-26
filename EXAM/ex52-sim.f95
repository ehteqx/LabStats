! EXAM EXERCISE N° 5.2 -- SIMULATION
! A straightforward, automated simulator for large scale iterated measurements
! (C) Emanuele Ballarin (ehteqx@gmail.com)
!##############################################################################

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

MODULE LSLMFIT         ! Module for the general-purpose calculation of best linear fits and associated uncertainties.

    use GLPREC
    implicit none

    CONTAINS

        SUBROUTINE LINFIT(dimarr, arr_x, arr_y, errvar, out_m, out_q, out_cov, out_corr, out_eq, out_em, linefit)		! Subroutine for the calculation of linear fits (LSLM)
							! PASS VARIANCE FROM THE PROGRAM
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

! ### END OF THE MODULES -- MAIN PROGRAM FOLLOWS ### !

PROGRAM SIMULATION

	use GLPREC
	use GLCONST
	use GAUSSIAN
	use LSLMFIT
	implicit none

	! Things necessary to generate the random seeds for gaussian error generation
	integer (kind = ik) :: cntsim = 0_ik 		! NUMBER OF SIMS
	integer (kind = 4) :: foresee_d				! SEED OF SEEDS
	real (kind = rk), allocatable :: aseeds_buffer = 0	! ARRAY OF RANDOM SEEDS (REAL BUFFER)
	integer (kind = 4), allocatable :: aseeds = 0	! ARRAY OF RANDOM SEEDS (INTEGER)

	! Values of m_0 and q_0 in [ y = (m_0)*x + (q_0) ]
	real (kind = rk) :: m_zero = pi, q_zero = 1.35_rk ! FIXED

	! Arrays of x-es and y-es
	real (kind = rk), dimension(7) :: arr_x = [1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk] ! FIXED
	real (kind = rk), dimension(7) :: arr_y = [0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk,0.0_rk] ! GENERATED

	! Array of variances and working array for variances (as in LINFIT)
	real (kind = rk), dimension(7) :: errvar = 0.05_rk, errvarprd = 0.0_rk ! FIXED --> GENERATED, ! FIXED

	! Things for calling LINFIT
	real (kind = rk), dimension(7) :: linefit		! LINEAR FIT
	integer (kind = ik) :: dimarr = 7				! LINEAR FIT
	real (kind = rk) :: out_m, out_q, out_cov, out_corr, out_eq, out_em ! LINEAR FIT

	! Some user interaction (for the random seed and the number of simulations)
	print*, 'Insert the number of simulations that you want to run... '
	read*, cntsim

	print*, 'Inserisci un numero di 7 (sette) cifre da usare come seme random... '
	read*, foresee_d

	! The seeds for gaussian error generation are in turn generated
	allocate(aseeds_buffer(cntsim))
	allocate(aseeds(cntsim))

	CALL RANDOM_SEED(foresee_d)
	CALL RANDOM_NUMBER(aseeds_buffer)

	aseeds_buffer = ((aseeds_buffer)*(1000000.0_rk))

	for ...

	! ### RANDOM ERROR GENERATION (MARSAGLIA'S KISS + BOX-MILLER TRANSFORM) ### !




END PROGRAM SIMULATION
