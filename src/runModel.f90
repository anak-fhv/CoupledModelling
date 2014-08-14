program runModel

	use mesh, only:	meshFile,meshNumNodes,meshNumElems,meshNumDoms,		&
				nBins,meshNumSurfs,meshVerts,meshElems,meshSurfs,		&
				readMesh,getElementNeighbours
	use fem
	use rt

	implicit none

	integer,parameter :: pDatFilNum=111
	integer :: values1(8),values2(8)
	character(*),parameter :: pDatFile="../data/probData.dat"
	character :: date*8,time*10,zone*5

	call date_and_time(date,time,zone,values1)
	write(*,'(a,(3i8,2x),i8)') "Start time: ", values1(5:8)
	open(pDatFilNum,file=pDatFile)
	read(pDatFilNum,*)
	read(pDatFilNum,*) meshFile
	read(pDatFilNum,*)
	read(pDatFilNum,*) nBins
	close(pDatFilNum)
	call readMesh()
	call getElementNeighbours()
	call date_and_time(date,time,zone,values2)
	write(*,'(a,(3i8,2x),i8)') "End time: ",values2(5:8)
	write(*,'(a,(3i8,2x),i8)') "Runtime: ",values2(5:8)-values1(5:8)
end program runModel
