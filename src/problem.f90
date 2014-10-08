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
		rtIterNum = 120
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


	subroutine rtReEmission()
		integer,parameter :: rtDatFileNum=102
		integer :: i,mainCtr,rtIterNum,nConstTSfs,mBinNum(3)
		integer,allocatable :: mConstTSfIds(:)
		real(8) :: totEmPow,rayPow,emSurfArea
		real(8),allocatable :: pRatio(:),mSfConstTs(:)
		character(*),parameter :: rtProbDatFile='../data/probData.dat'
		character(72) :: mFileName

!	Some hard coded values still exist
		rtIterNum = 80
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
			end if
			call traceOut(pRatio)
			rtElemSto = rtElemAbs
			rtElemAbs = 0
			write(*,'(a,2x,i4)')"Main iteration number: ", mainCtr
			write(*,'(a,2x,i8)')"Absorbed numbers: ", sum(rtElemSto)
		end do

	end subroutine rtReEmission

	subroutine rtLED()
		integer :: i,mainCtr,rtIterNum,mBinNum(3)
		real(8),allocatable :: pRatio(:)
		character(72) :: mFileName

		call readRtData(mBinNum,mSfConstTs,mSfConstQs,mFileName)
		rtBeta = rtKappa + rtSigma
		rtIterNum = 1
		do mainCtr = 1,rtIterNum
			if(mainCtr .eq. 1) then
				call rtInit(mBinNum,mSfConstTs,mSfConstQs,mFileName)
				allocate(pRatio(1))
				pRatio = (/1.d0/)
				rtRefRayPow = 1e-6
				call traceFromSurfLED(pRatio)
			end if
		end do	
	end subroutine rtLED

	subroutine readRtData(mBinNum,mSfConstTs,mSfConstQs,mFileName)
		integer,parameter :: rtDatFileNum=102
		integer :: nEmSfs,nConstTSfs,nConstQSfs,nTrSfs,nNonPartSfs,		&
		mBinNum(3)
		real(8),allocatable :: mSfConstTs(:),mSfConstQs(:)
		character(*),parameter :: rtProbDatFile='../data/ledData.dat'		
		character(72) :: mFileName

		open(rtDatFileNum,file=rtProbDatFile)
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) mFileName
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) mBinNum
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumEmSfs
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumCTSfs
		if(nConstTSfs .gt. 0) then
			allocate(rtConstTSfIds(nConstTSfs))
			allocate(mSfConstTs(nConstTSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtConstTSfIds
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) mSfConstTs
		else
			do i=1,4
				read(rtDatFileNum,*)
			end do
		end if
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumCQSfs
		if(nConstQSfs .gt. 0) then
			allocate(rtConstQSfIds(nConstQSfs))
			allocate(mSfConstQs(nConstQSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtConstQSfIds
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) mSfConstQs
		else
			do i=1,4
				read(rtDatFileNum,*)
			end do
		end if
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumTrSfs
		if(nTrSfs .gt. 0) then
			allocate (rtTrSfIds(nTrSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtTrSfIds
		else
			do i=1,2
				read(rtDatFileNum,*)
			end do
		end if
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumNpSfs
		if(nNonPartSfs .gt. 0) then
			allocate (rtNpSfIds(nNonPartSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtNpSfIds
		else
			do i=1,2
				read(rtDatFileNum,*)
			end do
		end if
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtKappa
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtSigma
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtRefrInd
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumRays
		close(rtDatFileNum)
	end subroutine readRtData

	subroutine rtFemSimple()
		integer,parameter :: rtDatFileNum=101,femProbFilNum=102
		integer :: i,mainCtr,rtIterNum,nConstTSfs,mBinNum(3)
		integer,allocatable :: mConstTSfIds(:)
		real(8) :: totEmPow,rayPow,emSurfArea
		real(8),allocatable :: pRatio(:),mSfConstTs(:),oldSrc(:),		&
		avSrc(:)
		character(*),parameter :: rtProbDatFile='../data/probData.dat',	&
		femProbFile='../data/femData.dat'
		character(72) :: mFileName,tempFileName,ctrString

!	Some hard coded values still exist
		rtIterNum = 80
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
				allocate(oldSrc(meshNumNodes))
				oldSrc = 0.0d0
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
			else
				call traceVolPowerBased()
				rtElemSto = rtElemAbs + rtWallInf
				rtElemAbs = 0
			end if
			call getAveragedSource(oldSrc,avSrc)
			rtNodalSrc = 0.0d0
			call setMeshNodalValues(avSrc,"S")
			call runFem(mFileName)
			oldSrc = avSrc
			deallocate(avSrc)
			write(*,'(a,2x,i4)')"Main iteration number: ", mainCtr
			write(*,'(a,2x,i8)')"Absorbed numbers: ", sum(rtElemSto)
		end do

	end subroutine rtFemSimple

	subroutine getAveragedSource(oldSrc,avSrc)
		real(8),parameter :: expAvFact = 0.05d0
		real(8),intent(in) :: oldSrc(:)
		real(8),allocatable,intent(out) :: avSrc(:)

		allocate(avSrc(meshNumNodes))
		avSrc = rtNodalSrc*expAvFact + oldSrc*(1.0d0-expAvFact)
	end subroutine getAveragedSource

end module problem