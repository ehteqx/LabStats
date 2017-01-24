! EXAM EXERCISE N° 5.2 New -- SIMULATION -- Version 3.1
! A straightforward, automated simulator for large scale iterated measurements
! (C) Emanuele Ballarin (ehteqx@gmail.com)
!*******************************************************************************

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

        SUBROUTINE POLARGAUSS(output, n, mean, stdev, rseed_a)	! Subroutine for the generation of Gausssian Distributions (Marsaglia Polar Method)

            implicit none

                integer (kind = ik), intent(in) :: n                         ! Dimension of the output array (must be user-passed)
                real (kind = rk), dimension(n), intent(out) :: output        ! Output array (contains the distribution)
                real (kind = rk), intent(in) :: mean, stdev                  ! Mean and std deviation of the distribution
                integer (kind = 4), intent(inout) :: rseed_a                 ! The random seeds used with the subroutine (change for different initial seed)

                integer (kind = 4) :: rseed_b, rseed_c
                real (kind = ik), dimension(n) :: randomarr_a, randomarr_b, randomarr_c ! Randomly-filled arrays
                integer (kind = ik) :: i, memo = 0_ik

                rseed_b = rseed_a + 1_ik
                rseed_c = rseed_a + 2_ik

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

MODULE CHISQUARED       ! Module for the computation of a Chi-Squared Distribution sampling

    use GLPREC
    implicit none

    CONTAINS

        FUNCTION XSQ(dglib, ics) result(ips)

            implicit none

                integer (kind = ik), intent(in) :: dglib
                real (kind = rk), intent(in) :: ics
                real (kind = rk) :: ips

                integer (kind = ik) :: w

				real (kind = rk) :: alpha
				real (kind = rk) :: beta
				real (kind = rk) :: gam

                ! DEFINITIONS
                alpha = ((real(dglib, rk))/2.0_rk)
                beta = (1.0_rk/2.0_rk)
                gam = GAMMA(alpha)

                ! RESULT (Chi-Squared distribution)
                ips = ((beta**alpha)*(ics**(alpha - 1.0_rk))*(EXP(-(beta*ics)))/(gam))

            END FUNCTION

    END MODULE CHISQUARED

! ******************************************************************************

! § MAIN PROGRAM (Nonmodular, Interactive, Complete)

PROGRAM SIMULATION

	use GLPREC
	use GLCONST
	use GAUSSIAN
	use LSLMFIT
	use CHISQUARED
	implicit none

    ! Indexes
    integer (kind = ik) :: p = 0, w = 0, k = 0, z = 0, h = 0, i = 0, j = 0, f = 0

    ! Values of m_0 and q_0 in [ y = (m_0)*x + (q_0) ]; Arbitrary, but fixed.
    real (kind = rk) :: m_zero = 1.0, q_zero = 8    ! DO NOT CHANGE DURING EXECUTION

    ! Gaussian seeds
    integer (kind = 4) :: fixseed1 = 5621478, fixseed2 =6953174, ofaseed = 5316985  ! OFA = One For All

    ! Chi-Squared/Gaussian testing arrays and variables
    real (kind = rk), dimension(:), allocatable :: true_xsq
    real (kind = rk), dimension(:), allocatable :: estm_xsq
    real (kind = rk), dimension(:), allocatable :: true_gss_m
    real (kind = rk), dimension(:), allocatable :: true_gss_q

    ! Result arrays (estimated m_0 and q_0 for each simulation)
    real (kind = rk), dimension(:), allocatable :: estimated_m
    real (kind = rk), dimension(:), allocatable :: estimated_q

	real (kind = rk), dimension(:), allocatable :: estimated_merr
	real (kind = rk), dimension(:), allocatable :: estimated_qerr
	real (kind = rk), dimension(:), allocatable :: correl

    ! Arrays of x-es and y-es
    real (kind = rk), dimension(7) :: arr_x = [1.0,2.0,3.0,4.0,5.0,6.0,7.0] ! Fixed sampling of x-es
    real (kind = rk), dimension(7) :: exact_y = [0.0,0.0,0.0,0.0,0.0,0.0,0.0] ! Initialized zeroes
    real (kind = rk), dimension(7) :: arr_y = [0.0,0.0,0.0,0.0,0.0,0.0,0.0] ! Initialized zeroes
    real (kind = rk), dimension(7) :: error_y = [0.0,0.0,0.0,0.0,0.0,0.0,0.0] ! Initialized zeroes
	real (kind = rk), dimension(7) :: stdevarr = [0.0,0.0,0.0,0.0,0.0,0.0,0.0] ! Initialized zeroes (STANDARD DEVS)
	real (kind = rk), dimension(7) :: variancevarr = [0.0,0.0,0.0,0.0,0.0,0.0,0.0] ! Initialized zeroes (VARIANCES)
	real (kind = rk), dimension(7) :: linearfitted = [0.0,0.0,0.0,0.0,0.0,0.0,0.0] ! Initialized zeroes

    ! Things necessary to generate the random seeds for gaussian error generation
    integer (kind = ik) :: cntsim = 0       ! Number of simulations (user-defined)
    integer (kind = 4) :: foresee_d		    ! Master seed (user-defined for randomness)

    ! Chi-Squared/Gaussian testing arrays and variables
    real (kind = rk) :: varp = 1.0, varm = 1.0, chisum = 0.0      ! DEVIAZIONE STANDARD? (Sono inseriti valori di test: 1.0)

	! Self-destroying variables
	real (kind = rk) :: dest1	!, dest2

    !...
    ! ################ INSERIRE LA LISTA DELLE VARIABILI DAL FILE A FIANCO ##############################################
    !...

    ! Bootstrapping message
    print*, ' '
    print*, '************************************************************************'
    print*, '                     SIMFIT - v. 3.1 GAUSSIAN ERRORS                   '
    print*, '                                                                       '
    print*, '                   Copyright (C) 2015 Emanuele Ballarin                '
    print*, '                                                                       '
    print*, ' This is free software, covered by the GNU General Public License v3,  '
    print*, ' and comes with ABSOLUTELY NO WARRANTY WHATSOEVER.                     '
    print*, ' For more information about the license: https://www.gnu.org/licenses/ '
    print*, '************************************************************************'
    print*, ' '
    print*, 'Program loaded!'
    print*, ' '

    ! Some user interaction (for the random seed and the number of simulations)
    print*, 'Insert the number of simulations that you want to run... '
    read*, cntsim
    print*, ' '

    print*, 'Inserisci un numero di 7 (sette) cifre da usare come seme random... '
    read*, foresee_d
    print*, ' '

    print*, 'Simulation started. Please wait, it may take some time to complete...'
    print*, ' '
    print*, 'Using', m_zero, 'as m_0 central value'
    print*, 'Using', q_zero, 'as q_0 central value'
    print*, ' '
    print*, 'Computing...'

    ! Simulation started...
    exact_y = (((m_zero)*(arr_x)) + (q_zero))		! Array initialization (exact y)
    ofaseed = foresee_d

	do k = 1, 7, 1
		stdevarr(k) = (0.05_rk)*exact_y(k)
		variancevarr(k) = stdevarr(k)**2
	end do

	!Allocating allocatable arrays
    allocate(estimated_m(cntsim))
        estimated_m = 0.0_rk
    allocate(estimated_q(cntsim))
        estimated_q = 0.0_rk
	allocate(estimated_merr(cntsim))
        estimated_m = 0.0_rk
    allocate(estimated_qerr(cntsim))
        estimated_q = 0.0_rk
	allocate(correl(cntsim))
        estimated_q = 0.0_rk
    allocate(estm_xsq(cntsim))
        estm_xsq = 0.0_rk
    allocate(true_xsq((10_ik)*cntsim))
        true_xsq = 0.0_rk
    allocate(true_gss_m((10_ik)*cntsim))
        true_gss_m = 0.0_rk
    allocate(true_gss_q((10_ik)*cntsim))
        true_gss_q = 0.0_rk

	! Outer cycle of simulations
    do w = 1, cntsim, 1

        ! Generating errors (of stated stdev) for simulated y
        do z = 1, 7, 1
            CALL POLARGAUSS(error_y(z), 1_ik, 0.0_rk, stdevarr(z), ofaseed)
            ofaseed = ofaseed + 1
        end do

        arr_y = exact_y + error_y   ! Computing simulations as EXACT + ERROR

        CALL LINFIT(7_ik, arr_x, arr_y, variancevarr, estimated_m(w), estimated_q(w),&		! On two lines because...
		& dest1, correl(w), estimated_qerr(w), estimated_merr(w), linearfitted)					! ...else the line is truncated

		!Chi-Squared testing (REDUCED)
		chisum = 0.0_rk
		do h = 1, 7, 1
		    chisum = chisum + (((linearfitted(h) - arr_y(h))**2)/(variancevarr(h)))
		end do
		estm_xsq(w) = chisum

	end do

    ! User notification
    print*, 'OK.'
    print*, ' '
    print*, 'Generating distributions and testing hypoteses...'

    ! Chi-Squared testing
    do p = 1, cntsim*10_ik, 1
        true_xsq(p) = XSQ(5_ik, (real(p, rk)))
    end do

    ! Gaussian testing (m)
    CALL POLARGAUSS(true_gss_m, ((10_ik)*cntsim), m_zero, varm, fixseed1)       ! DEVIAZIONE STANDARD?

    ! Gaussian testing (q)
    CALL POLARGAUSS(true_gss_q, ((10_ik)*cntsim), q_zero, varp, fixseed2)       ! DEVIAZIONE STANDARD?

    ! User notification
    print*, 'OK.'
    print*, ' '
    print*, 'Writing data to file...'

	! Writing data to file the necessary data
	open(unit=2, file='results.dat')
		do i = 1, cntsim, 1
			write(unit=2,fmt=*)i, estimated_m(i), estimated_q(i), estm_xsq(i)
		end do

	open(unit=3, file='distrib.dat')
		do j = 1, 10*cntsim, 1
			write(unit=3,fmt=*)j, true_gss_m(j), true_gss_q(j), true_xsq(j)
		end do

	! Chi-Squared internal ellipse (direct print) <---- ????
	open(unit=4, file='xq.dat')
		do f = 1, cntsim, 1
			if (((1.0_rk/(1.0_rk - (correl(f)**2)))*(((estimated_m(f)&
			& - m_zero)**2)/estimated_merr(f)) + (((estimated_q(f)&
			& - q_zero)**2)/estimated_qerr(f))&
			& - ((0.5_rk*correl(f)*(estimated_m(f) - m_zero)*(estimated_q(f)&
			& - q_zero))&
			& / ((sqrt(estimated_merr(f)*estimated_qerr(f))))))&
			& .LE. 1.0_rk ) then
				write(unit=4,fmt=*)f, estimated_m(f), estimated_q(f)
			end if

				!print*, correl(f), estimated_merr(f), estimated_qerr(f)

		end do

	! User notification
	print*, 'OK.'
	print*, ' '
	print*, ' '
	print*, 'Bye bye, have a nice day!'
	print*, ' '

END PROGRAM SIMULATION