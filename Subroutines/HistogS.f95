! HISTOG - v. 3.1
! An enhanced toolset for manipulating and plotting histogram data in Fortran 95 (SubRoutinized Version)

! (C) Emanuele Ballarin (ehteqx@gmail.com)
!########################################################################################################

MODULE GLPREC			                 ! Module for global precision constants

	implicit none

	integer, parameter 	:: ik = selected_int_kind(38)		! MAX: 38
	integer, parameter 	:: rk = selected_real_kind(33)		! MAX: 33

END MODULE GLPREC

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
