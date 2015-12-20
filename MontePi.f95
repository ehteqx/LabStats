! MONTEPI - v. 1.2
! A multiprecision calculator of the Pi constant using a simple pseudo-ramdomic Monte Carlo algorithm.

! It relies heavily on George Marsaglia's KISS linear RNG (Fortran built-in).
! It uses the GLPREC and GLCONST modules/subroutines by Emanuele Ballarin.

! (C) Emanuele Ballarin (ehteqx@gmail.com)
!######################################################################################################

! ### MODULES IMPORT & DEFINITION ###

MODULE GLPREC			                 ! Module for global precision constants

	implicit none

	integer, parameter 	:: ik = selected_int_kind(38)		! MAX: 38
	integer, parameter 	:: rk = selected_real_kind(33)		! MAX: 33

END MODULE GLPREC

MODULE GLCONST			! Module for global and universal constants

        use GLPREC
        implicit none

    	real (kind = rk), parameter 	:: onepi = acos(-1.0_rk)  ! Pi; max. precision available

    END MODULE GLCONST


PROGRAM MONTEPI

	use GLPREC				! All the required modules are imported
	use GLCONST
	implicit none			! As it is standard-prescribed

	! ## VARIABLES DEFINITION ##

    integer (kind = ik) :: points		! Number of points used for the MC simulation

    real (kind = rk), allocatable :: xcoord(:), ycoord(:)     ! Arrays of coordinates (from 0.0 to 1.0)
    integer (kind = 4) :: rseed_x = 5931041, rseed_y = 6327490  ! The random seeds used to generate the coordinates

    real (kind = rk), allocatable :: values(:), errors(:)     ! Working arrays
    integer (kind = ik) :: inside = 0_ik    ! Counter

    real (kind = rk) :: pigreco = 0.0_rk

	integer (kind = ik) :: i, j   ! Indexes (some can be unnecessary)

	! ## INTRO SCREEN ##

	print*, ' '
	print*, '#######################################################################'
	print*, '                           MONTEPI - v. 1.2                           '
	print*, '                                                                      '
	print*, '                  Copyright (C) 2015 Emanuele Ballarin                '
	print*, '                                                                      '
	print*, 'This is free software, covered by the GNU General Public License v3,  '
	print*, 'and comes with ABSOLUTELY NO WARRANTY WHATSOEVER.                     '
	print*, 'For more information about the license: https://www.gnu.org/licenses/ '
	print*, '#######################################################################'
	print*, ' '

	! ## USER INTERACTION ##

	print*, 'Inserisci il numero di punti da utilizzare per la simulazione Monte-Carlo (the more, the more precise) '
	read*, points

    print*, 'Computation started...'		! Some info for the user

    ! ## PREPARING FOR RUN ##

    allocate(xcoord(points))        ! The coordinates arrays are allocated using the user-input number of points
    allocate(ycoord(points))

    CALL RANDOM_SEED(rseed_x)       ! The arrays are filled with random coordinates, between 0.0 and 1.0
    CALL RANDOM_NUMBER(xcoord)

    CALL RANDOM_SEED(rseed_y)
    CALL RANDOM_NUMBER(ycoord)

    allocate(values(points))
    allocate(errors(points))

    do i = 1_ik, points

        if ((((xcoord(i))**2) + ((ycoord(i))**2)) .LE. 1.0_rk) then

            inside = inside + 1_ik

        end if

        values(i) = ((4.0_rk)*(REAL(inside, rk))/(REAL(i, rk)))
        errors(i) = ABS((values(i)) - onepi)

    end do

    pigreco = values(points)

    print*, 'Computation completed!'		! Some info for the user
    print*, ' '

    print*, 'Copying to file...'		! Some info for the user

    open(unit=300, file='pi.dat')
    write(unit=300,fmt=*)pigreco

    print*, 'Do you want to copy to file all the intermediate values and errors? (Type 'yes' or 'no')'		! Some info for the user

    open(unit=100, file='values.dat')
    open(unit=200, file='errors.dat')
    write(unit=100,fmt=*)values
    write(unit=200,fmt=*)errors


    print*, 'DONE!'		! Some info for the user
    print*, ' '
    print*, ' '

    print*, 'The estimated value for Pi, using ', points, ' points is: ', pigreco	! Some info for the user

END PROGRAM MONTEPI

! ### END OF THE FILE
