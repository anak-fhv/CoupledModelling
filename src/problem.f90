module problem

	use fem
	use rt
	use mesh
	use utilities

	implicit none

	contains

	subroutine femSimple()
		integer,parameter :: femProbFilNum=102
		integer,allocatable :: mConstTSfIds(:)
		character(*),parameter :: probFile='../data/femData.dat'
		character(72) :: mFileName

		open(femProbFilNum,file=probFile)
		read(femProbFilNum,*)
		read(femProbFilNum,*) mFileName
		read(femProbFilNum,*)
		read(femProbFilNum,*) femByFile
		read(femProbFilNum,*)
		read(femProbFilNum,*) femResFile
		read(femProbFilNum,*)
		read(femProbFilNum,*) femTransient
		read(femProbFilNum,*)
		read(femProbFilNum,*) femSolverMaxIter
		read(femProbFilNum,*)
		read(femProbFilNum,*) femSolverType
		if(femTransient) then
			read(femProbFilNum,*)
			read(femProbFilNum,*) femTrScheme			
		end if

		call runFem(mFileName)
	end subroutine femSimple

	subroutine rtSimple()
		integer,parameter :: rtDatFileNum=102
		integer :: i,mainCtr,rtIterNum,nConstTSfs,mBinNum(3)
		integer,allocatable :: mConstTSfIds(:)
		real(8) :: totEmPow,rayPow,emSurfArea
		real(8),allocatable :: pRatio(:),mSfConstTs(:)
		character(*),parameter :: rtProbDatFile='../data/probData.dat'
		character(72) :: mFileName

!	Some hard coded values still exist
		rtIterNum = 16
		emSurfArea = 4.0d0
!	Need to rid ourselves of these pesky constants

		open(rtDatFileNum,file=rtProbDatFile)
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) mFileName
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) mBinNum
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) nConstTSfs
		allocate(rtEmSfIds(nConstTSfs))
		allocate(mSfConstTs(nConstTSfs))
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtEmSfIds
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) mSfConstTs
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumRays
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtKappa
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtSigma
		close(rtDatFileNum)

		do mainCtr = 1,rtIterNum
			if(mainCtr .eq. 1) then
				call rtInitMesh(mFileName,mBinNum,mConstTSfIds,			&
				mSfConstTs)
				nConstTSfs = size(mSfConstTs,1)
				allocate(pRatio(nConstTSfs))
				totEmPow = sum(mSfConstTs**4.0d0)
				do i=1,nConstTSfs-1
					pRatio(i) = sum(mSfConstTs(1:i)**4.0d0)/totEmPow
				end do
				pRatio(nConstTSfs) = 1.0d0
				write(*,*) "pRatio: ", pRatio
				rtRefRayPow = sigB*totEmPow*emSurfArea/rtNumRays
				write(*,*) "rayPow: ", rayPow
				call traceFromSurf(pRatio)
				rtElemSto = rtElemAbs
				rtWallInf = rtElemAbs
				rtElemAbs = 0
			else
				call traceFromVol()
				rtElemSto = rtElemAbs + rtWallInf
				rtElemAbs = 0
			end if
			write(*,'(a,2x,i4)')"Main iteration number: ", mainCtr
			write(*,'(a,2x,i8)')"Absorbed numbers: ", sum(rtElemSto)
		end do
	end subroutine rtSimple

	subroutine rtFemSimple()
		integer,parameter :: rtDatFileNum=101,femProbFilNum=102
		integer :: i,mainCtr,rtIterNum,nConstTSfs,mBinNum(3)
		integer,allocatable :: mConstTSfIds(:)
		real(8) :: totEmPow,rayPow,emSurfArea
		real(8),allocatable :: pRatio(:),mSfConstTs(:)
		character(*),parameter :: rtProbDatFile='../data/probData.dat',	&
		femProbFile='../data/femData.dat'
		character(72) :: mFileName,tempFileName,ctrString

!	Some hard coded values still exist
		rtIterNum = 40
		emSurfArea = 4.0d0
!	Need to rid ourselves of these pesky constants

!	Read data for the RT and FEM problems from separate files (temporary)
!	Read RT data
		open(rtDatFileNum,file=rtProbDatFile)
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) mFileName
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) mBinNum
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) nConstTSfs
		allocate(rtEmSfIds(nConstTSfs))
		allocate(mSfConstTs(nConstTSfs))
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtEmSfIds
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) mSfConstTs
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumRays
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtKappa
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtSigma
		close(rtDatFileNum)

!	Read FEM data
		open(femProbFilNum,file=femProbFile)
		read(femProbFilNum,*)
		read(femProbFilNum,*) mFileName
		read(femProbFilNum,*)
		read(femProbFilNum,*) femByFile
		read(femProbFilNum,*)
		read(femProbFilNum,*) femResFile
		read(femProbFilNum,*)
		read(femProbFilNum,*) femTransient
		read(femProbFilNum,*)
		read(femProbFilNum,*) femSolverMaxIter
		read(femProbFilNum,*)
		read(femProbFilNum,*) femSolverType
		if(femTransient) then
			read(femProbFilNum,*)
			read(femProbFilNum,*) femTrScheme			
		end if

		tempFileName = femResFile
		do mainCtr = 1,rtIterNum
			write(ctrString,'(i4.4)') mainCtr
			femResFile = trim(adjustl(tempFileName))//trim(adjustl(ctrString))
			if(mainCtr .eq. 1) then
				call rtInitMesh(mFileName,mBinNum,mConstTSfIds,			&
				mSfConstTs)
				nConstTSfs = size(mSfConstTs,1)
				allocate(pRatio(nConstTSfs))
				totEmPow = sum(mSfConstTs**4.0d0)
				do i=1,nConstTSfs-1
					pRatio(i) = sum(mSfConstTs(1:i)**4.0d0)/totEmPow
				end do
				pRatio(nConstTSfs) = 1.0d0
				write(*,*) "pRatio: ", pRatio
				rtRefRayPow = sigB*totEmPow*emSurfArea/rtNumRays
				write(*,*) "rayPow: ", rtRefRayPow
				call traceFromSurf(pRatio)
				rtElemSto = rtElemAbs
				rtWallInf = rtElemAbs
				rtElemAbs = 0
				call runFem(mFileName)
			else
				call traceVolPowerBased()
				rtElemSto = rtElemAbs + rtWallInf
				rtElemAbs = 0
				call runFem(mFileName)
			end if
			write(*,'(a,2x,i4)')"Main iteration number: ", mainCtr
			write(*,'(a,2x,i8)')"Absorbed numbers: ", sum(rtElemSto)
		end do

	end subroutine rtFemSimple

end module problem
