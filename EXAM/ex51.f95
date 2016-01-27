! EXAM EXERCISE N° 5.1 --- (Using: LinearFit - v. 1.2 ENH)
! (C) Emanuele Ballarin (ehteqx@gmail.com)
!##########################################################

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
                    errvar(k) = (((errvar(k))*(arr_y(k))/(sqrt(3.0_rk)))**(2.0_ik))
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

PROGRAM EX51

        use GLPREC
        use LSLMFIT
        implicit none

		real (kind = rk), dimension(7) :: arr_x = [1.0_rk,2.0_rk,3.0_rk,4.0_rk,5.0_rk,6.0_rk,7.0_rk]
		real (kind = rk), dimension(7) :: arr_y = [3.2_rk,5.3_rk,6.5_rk,9.0_rk,10.3_rk,13.4_rk,14.5_rk]
		real (kind = rk), dimension(7) :: errvar = 0.05_rk, errvarprd = 0.0_rk
		real (kind = rk), dimension(7) :: linefit
		integer (kind = ik) :: dimarr = 7
		real (kind = rk) :: out_m, out_q, out_cov, out_corr, out_eq, out_em

        errvarprd = errvar

		CALL LINFIT(dimarr, arr_x, arr_y, errvarprd, out_m, out_q, out_cov, out_corr, out_eq, out_em, linefit)

		print*, dimarr, arr_x, arr_y, errvar, errvarprd, out_m, out_q, out_cov, out_corr, out_eq, out_em, linefit

    END PROGRAM EX51
