program runModel

	use problem

	implicit none

	integer :: i,ct,values1(8),values2(8)
	character :: dt*8,tm*10,zn*5,fMesh*32,fRt*32,fFem*32

	ct = command_argument_count()
	if(ct .lt. 3) then
		write(*,'(a)')"Please provide all data file names."
		write(*,'(a)')"Retry. Exiting now..."
		stop
	end if
	if(ct .gt. 3) then
		write(*,'(a)')"Superfluous arguments at the end. Discarded."
		ct = 3
	end if
	call get_command_argument(1,fMesh)
	write (*,'(a)') "Mesh data file: ", trim(fMesh)
	call get_command_argument(2,fRt)
	write (*,'(a)') "MCRT data file: ", trim(fRt)
	call get_command_argument(3,fFem)
	write (*,'(a)') "FEM data file: ", trim(fFem)

	call date_and_time(dt,tm,zn,values1)
	write(*,'(a,(3i8,2x),i8)') "Start time: ", values1(5:8)

	call runCoupledModel(fMesh,fRt,fFem)

	call date_and_time(dt,tm,zn,values2)
	write(*,'(a,(3i8,2x),i8)') "End time: ",values2(5:8)

end program runModel
