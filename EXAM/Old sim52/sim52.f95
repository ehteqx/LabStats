    do w = 1, cntsim, 1 ! Cycle of simulations

		do ciccio = 1, 7, 1
			arr_y(ciccio) = (((m_zero)*(arr_x(ciccio))) + (q_zero)) ! The y-es array is (re)initialized
		end do

        ! Generating errors...
        CALL POLARGAUSS(error_y, 7_ik, 0.0_rk, 0.05_rk, aseeds(memory + 1_ik), aseeds(memory + 2_ik), aseeds(memory + 3_ik))

		arr_y_bu = arr_y		  ! Copy to backup
		arr_y = (arr_y + error_y) ! Adding errors to the y-es array

		variarray = ((0.05_rk)**2_ik)
		do pippo = 1, 7, 1
			variarrayos(pippo) = variarray(pippo)*(arr_y_bu(pippo))**2
		end do

        CALL LINFIT(7_ik, arr_x, arr_y, variarrayos, out_m, out_q, out_cov, out_corr, pizza, pizzas, pit)

		estimated_m(w) = out_m
		estimated_q(w) = out_q

		do pippo = 1,7,1
		!variarrayos(pippo) = pizzas*((arr_x(pippo))**2) + pizza
		!variarrayos
		end do

		! Precomputations for Chi-Squared testing
		chisum = 0.0_rk

		do p = 1, 7, 1
		    chisum = chisum + ((pit(p) - arr_y_bu(p))**2)/(variarrayos(p))
		end do

		estm_xsq(w) = chisum

		! For the RNG
        memory = memory + 3_ik ! Info for the random numbers generator

    end do

!***************************************************************************************************************

	! Writing data to file (only the necessary data)
	open(unit=1, file='exacts.dat')
		write(unit=1,fmt=*)'VERTICALLY ALIGNED'
		write(unit=1,fmt=*)'m vero', 'q vero', 'domain; sampling', 'incertezza y'
		write(unit=1,fmt=*)m_zero, q_zero, '[0.0, 7.0]; 1.0', '0.05'

	open(unit=2, file='results.dat')
		write(unit=2,fmt=*)'HORIZONTALLY ALIGNED'
		write(unit=2,fmt=*)'Sim. n°', 'm stim.', 'q stim', 'X2 risult.'
		do i = 1, cntsim, 1
			write(unit=2,fmt=*)i, estimated_m(i), estimated_q(i), estm_xsq(i)
		end do

	open(unit=3, file='distrib.dat')
		write(unit=3,fmt=*)'HORIZONTALLY ALIGNED'
		write(unit=3,fmt=*)'Val. n°', 'gaussiana m', 'gaussiana q', 'X2 n = 6'
		do i = 1, 10*cntsim, 1
			write(unit=3,fmt=*)i, true_gss_m(i), true_gss_q(i), true_xsq(i)
		end do

	! User notification
    print*, 'OK.'
    print*, ' '
    print*, ' '
    print*, 'Bye bye, have a nice day!'
    print*, ' '

END PROGRAM SIMULATION
