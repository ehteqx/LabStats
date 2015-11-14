! LinearFit - v. 1.2 ENH
! A simple CLI software for computing the best linear fit (and associated uncertainties) of a data sample using the Least Squares method.
! (SubRoutinized version)
! (C) Emanuele Ballarin (ehteqx@gmail.com) -- 13/11/2015
!########################################################################################################################################

MODULE GLPREC			! Module for global precision constants

	implicit none

	integer, parameter 	:: ik = selected_int_kind(38)		! MAX: 38
	integer, parameter 	:: rk = selected_real_kind(33)		! MAX: 33

    END MODULE GLPREC

MODULE LSLMFIT         ! Module for the general-purpose calculation of best linear fits and associated uncertainties.

    use GLPREC
    implicit none

    CONTAINS

        SUBROUTINE LINFIT(dimarr, arr_x, arr_y, errvar, out_m, out_q, out_cov, out_corr, out_eq, out_em, linefit)		! Subroutine for the calculation of linear fits (LSLM)
							! PASS VARIANCE FROM THE PROGRAM
            implicit none

				integer (kind = ik), intent(in) :: dimarr
				real (kind = rk), dimension(dimarr), intent(in) :: arr_x, arr_y, errvar
				real (kind = rk), dimension(dimarr), intent(out) :: linefit
				real (kind = rk), intent(out) :: out_m, out_q, out_cov, out_corr, out_eq, out_em

				real (kind = rk) :: s0, s1, s2, s3, s4, s5, d
				integer (kind = ik) :: i

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

! ### HERE THE MODULE FINISHES. THE PROGRAM PROVIDED BELOW MUST BE USED FOR TESTING PURPOSES ONLY ###

PROGRAM TEST

        use GLPREC
        use LSLMFIT
        implicit none

		real (kind = rk), dimension(5) :: arr_x = [1,2,3,4,5], arr_y = [1.1,2.1,3.1,3.9,4.8], errvar = 0.4
		real (kind = rk), dimension(5) :: linefit
		integer (kind = ik) :: dimarr = 5
		real (kind = rk) :: out_m, out_q, out_cov, out_corr, out_eq, out_em

		CALL LINFIT(dimarr, arr_x, arr_y, errvar, out_m, out_q, out_cov, out_corr, out_eq, out_em, linefit)

		print*, dimarr, arr_x, arr_y, errvar, out_m, out_q, out_cov, out_corr, out_eq, out_em, linefit

    END PROGRAM TEST
