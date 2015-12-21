! BUFFONEEDLE - v. 1.2 SIM
! A pseudo-ramdomic Monte Carlo simulator of Buffon's Needle Problem.

! It relies heavily on George Marsaglia's KISS linear RNG (Fortran built-in).
! It uses the GLPREC and GLCONST modules/subroutines by Emanuele Ballarin.

! (C) Emanuele Ballarin (ehteqx@gmail.com)
!###############################################################################

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


PROGRAM BUFFONEEDLE

	use GLPREC				! All the required modules are imported
	use GLCONST
	implicit none			! As it is standard-prescribed

	! ## VARIABLES DEFINITION ##

    integer (kind = ik) :: throws		! Number of throws used to simulate the event
	real (kind = rk) :: floor, needle, ratio	! Lenght of the needle and the distance detween tiles

    real (kind = rk), allocatable :: dist(:), angle(:)     ! Arrays of coordinates (distance, angle)
    integer (kind = 4) :: rseed_x = 6187321, rseed_y = 2336850  ! The random seeds used to generate the coordinates

    real (kind = rk), allocatable :: values(:), errors(:)     ! Working arrays
    integer (kind = ik) :: crosses = 0_ik    ! Counter

    real (kind = rk) :: pigreco = 0.0_rk

	integer (kind = ik) :: i, j   ! Indexes

	character :: check				! Checker for user interaction
	logical :: lequal = (.FALSE.)		! Boolean checkbox for hypothesis verification

	! ## INTRO SCREEN ##

	print*, ' '
	print*, '#######################################################################'
	print*, '                         BUFFONEEDLE - v. 1.2                         '
	print*, '                                                                      '
	print*, '                  Copyright (C) 2015 Emanuele Ballarin                '
	print*, '                                                                      '
	print*, 'This is free software, covered by the GNU General Public License v3,  '
	print*, 'and comes with ABSOLUTELY NO WARRANTY WHATSOEVER.                     '
	print*, 'For more information about the license: https://www.gnu.org/licenses/ '
	print*, '#######################################################################'
	print*, ' '

	! ## USER INTERACTION ##

	print*, 'Inserisci il numero di lanci da utilizzare per la simulazione Monte-Carlo (the more, the more precise) '
	read*, throws

	print*, 'Inserisci la distanza tra le fughe parallele del pavimento (valore reale arbitrario) '
	read*, floor

	do while (lequal .EQV. (.FALSE.))			! While cycle to have the needle length smaller than the tile distance
		print*, 'Inserisci la lunghezza dello ago (minore o uguale della distanza tra le fughe) '
		read*, needle

		if (needle .LE. floor) then

			lequal = .TRUE.

		else
			lequal = .FALSE.
			print*, 'Non hai inserito un valore consentito! Riprova...'

		end if

	end do

    print*, 'Computation started...'		! Some info for the user

	ratio = floor/needle

    ! ## PREPARING FOR RUN ##

    allocate(dist(throws))        ! The coordinates arrays are allocated using the user-input number of points
    allocate(angle(throws))

    CALL RANDOM_SEED(rseed_x)       ! The arrays are filled with random coordinates, between 0.0 and 1.0
    CALL RANDOM_NUMBER(dist)
	dist = (dist)*floor/2.0			! The array coordinates span is adjusted according to user input

    CALL RANDOM_SEED(rseed_y)
    CALL RANDOM_NUMBER(angle)
	angle = (angle)*onepi/2.0		! The array coordinates span is adjusted

    allocate(values(throws))
    allocate(errors(throws))

    do i = 1_ik, throws

        if ((dist(i)) .LE. ((needle)*(sin(angle(i)))/2.0)) then

            crosses = crosses + 1_ik

        end if

        values(i) = 2.0/(((REAL(crosses, rk))/(REAL(i, rk)))*ratio)
        errors(i) = ABS((values(i)) - onepi)

    end do

    pigreco = values(throws)

    print*, 'Computation completed!'		! Some info for the user
    print*, ' '

    print*, 'Copying to file...'		! Some info for the user

    open(unit=300, file='pi.dat')
    write(unit=300,fmt=*)pigreco

    print*, 'Do you want to copy to file the temp. values and errors? (Type y / n)'	! Some info for the user
    read*, check

    if (check .EQ. 'y') then

    	open(unit=100, file='values.dat')
    	open(unit=200, file='errors.dat')
    	write(unit=100,fmt=*)values
    	write(unit=200,fmt=*)errors

    end if

    print*, 'DONE!'		! Some info for the user
    print*, ' '
    print*, ' '

    print*, 'The estimated value for Pi, using ', throws, ' throws is: ', pigreco	! Some info for the user

END PROGRAM BUFFONEEDLE

! ### END OF THE FILE
