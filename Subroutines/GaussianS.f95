! POLARGAUSS - v. 1.1 ENH (SubRoutinized Version)
! An optimized implementation of Box-Muller/George Marsaglia's Polar Method for Gaussian Distributions random sampling.
! (C) Emanuele Ballarin (ehteqx@gmail.com) -- 13/11/2015
!######################################################################################################################

MODULE GLPREC			! Module for global precision constants

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

! ### HERE THE MODULE FINISHES. THE PROGRAM PROVIDED BELOW MUST BE USED FOR TESTING PURPOSES ONLY ###

PROGRAM TEST

        use GLPREC
        use GLCONST
        use GAUSSIAN
        implicit none

        integer (kind = ik) :: n = 1000000              ! Modify the dimension as needed...
        real (kind = rk), dimension(1000000) :: output  ! The value of n must be the same as the dimension of the array.

        CALL POLARGAUSS(output, n, 0.0_rk, 1.0_rk)      ! The subroutine POLARGAUSS is called.

        !print*, output                                 ! Uncomment the line if you want to print the output array on screen.

    END PROGRAM TEST
