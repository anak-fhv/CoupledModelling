program runModel

	use mesh, only:	meshFile,meshNumNodes,meshNumElems,meshNumDoms,		&
				nBins,meshNumSurfs,meshVerts,meshElems,meshSurfs,		&
				readMesh,getElementNeighbours
	use fem
	use rt

	implicit none

	integer,parameter :: pDatFilNum=111,trOutFilNum=222
	integer :: i,values1(8),values2(8),nHood(4,2)
	character(*),parameter :: pDatFile="../data/probData.dat",			&
	trOutFile="../data/tetraData.out"
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
	open(trOutFilNum,file=trOutFile)
	do i=1,meshNumElems
		nHood = meshElems(i)%neighbours
		write(trOutFilNum,'(1x,4(i8,1x,i4))') transpose(nHood)
	end do
	call date_and_time(date,time,zone,values2)
	write(*,'(a,(3i8,2x),i8)') "End time: ",values2(5:8)
	write(*,'(a,(3i8,2x),i8)') "Runtime: ",values2(5:8)-values1(5:8)
end program runModel
