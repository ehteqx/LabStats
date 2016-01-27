MODULE CHISQUARED         ! Module for the computation of a Chi-Squared Distribution sampling

    use GLPREC
    implicit none

    CONTAINS

        FUNCTION XSQ(dglib, ics), result(ips)

            implicit none

                integer (kind = ik), intent(in) :: dglib
                real (kind = rk), intent(in) :: ics
                real (kind = rk), intent(out) :: ips

                integer (kind = ik) :: w

                ! DEFINITIONS
                real (kind = rk) :: alpha = ((real(dglib, rk))/2.0_rk)
                real (kind = rk) :: beta = (1.0_rk/2.0_rk)
                real (kind = rk) :: gam = GAMMA(alpha)

                ! RESULT (Chi-Squared distribution)
                ips = ((beta**alpha)*(ics**(alpha - 1.0_rk))*(EXP(-(beta*ics)))/(gam))

            END FUNCTION

    END MODULE CHISQUARED
