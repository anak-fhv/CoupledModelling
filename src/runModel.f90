program runModel

	use problem

	implicit none

	integer :: i,ct,values1(8),values2(8)
	character :: simType*2,dt*8,tm*10,zn*5,fMesh*32,fRt*32,fFem*32

	call date_and_time(dt,tm,zn,values1)
	write(*,'(a,(3i8,2x),i8)') "Start time: ", values1(5:8)

	ct = command_argument_count()
	if(ct .lt. 3) then
		write(*,'(a)')"Please provide all data, including file names."
		write(*,'(a)')"Retry. Exiting now..."
		stop
	end if

	call get_command_argument(1,simType)

	if(simType .eq. "RT") then
		write(*,'(/a)') "Only Ray Tracing run chosen."
		call get_command_argument(2,fMesh)
		write (*,'(a,2x,a)') "Mesh data file: ", trim(fMesh)
		call get_command_argument(3,fRt)
		write (*,'(a,2x,a)') "MCRT data file: ", trim(fRt)
		call rtSimple(fMesh,fRt)
	elseif(simType .eq. "FE") then
		write(*,'(/a)') "Only FEM solution sought."
		call get_command_argument(2,fMesh)
		write (*,'(a,2x,a)') "Mesh data file: ", trim(fMesh)
		call get_command_argument(3,fFem)
		write (*,'(a,2x,a)') "FEM data file: ", trim(fFem)
		call femSimple(fMesh,fFem)
	elseif(simType .eq. "CM") then
		write(*,'(/a)') "Coupled RT-FEM simulation selected."
		if(ct .lt. 4) then
			write(*,'(a)')"Some data still missing, please retry."
			write(*,'(a)')"Exiting now..."
			stop
		elseif(ct .gt. 4) then
			write(*,'(a)')"Superfluous arguments discarded."
		else
			continue
		end if
		call get_command_argument(2,fMesh)
		write (*,'(a,2x,a)') "Mesh data file: ", trim(fMesh)
		call get_command_argument(3,fRt)
		write (*,'(a,2x,a)') "MCRT data file: ", trim(fRt)
		call get_command_argument(4,fFem)
		write (*,'(a,2x,a)') "FEM data file: ", trim(fFem)
		call runCoupledModel(fMesh,fRt,fFem)
	else
		write(*,'(a)')"Simulation type chosen not recognised."		
	end if

	call date_and_time(dt,tm,zn,values2)
	write(*,'(a,(3i8,2x),i8)') "End time: ",values2(5:8)

end program runModel
