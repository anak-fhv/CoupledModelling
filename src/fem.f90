module fem

	use femHelper
	implicit none

	contains

	subroutine femSimple()
		integer,parameter :: probFilNum=102
		integer,allocatable :: mConstTSfIds(:)
		character(*),parameter :: probFile='../data/femData.dat'
		character(72) :: mFileName

		open(probFilNum,file=probFile)
		read(probFilNum,*)
		read(probFilNum,*) mFileName
		read(probFilNum,*)
		read(probFilNum,*) femByFile
		read(probFilNum,*)
		read(probFilNum,*) femResFile
		read(probFilNum,*)
		read(probFilNum,*) femTransient
		read(probFilNum,*)
		read(probFilNum,*) femSolverMaxIter
		read(probFilNum,*)
		read(probFilNum,*) femSolverType
		if(femTransient) then
			read(probFilNum,*)
			read(probFilNum,*) femTrScheme			
		end if

		call runFem(mFileName)
	end subroutine femSimple

end module fem
