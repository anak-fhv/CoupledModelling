module problem

	use fem
	use rt
	use mesh
	use utilities

	implicit none

	contains

	subroutine runCoupledModel(fMesh,fRt,fFem)
		integer :: domChk,emSfChk,bySfChk,chkNumSfs,mBinNum(3)
		character(32),intent(in) :: fMesh,fRt,fFem
		character(32) :: fMeshDat,fRtDat,fFemDat
		character(32) :: mFileName
		real(8),allocatable :: pRatio(:)

		fMeshDat = commDatDir//trim(adjustl(fMesh))//commDatExt
		fRtDat = commDatDir//trim(adjustl(fRt))//commDatExt
		fFemDat = commDatDir//trim(adjustl(fFem))//commDatExt
		
		call readMeshDataFile(fMeshDat)
		call readRtDataFile(fRtDat)
		call readFemDataFile(fFemDat)
		call rtInit()
!		allocate(pRatio(1))
		allocate(rtSfPowRatio(rtNumEmSfs))
!		pRatio = (/1.d0/)
		rtSfPowRatio = (/1.d0/)
!		rtRefRayPow = 0.36d0/real(rtMCNumRays,8)	! Hard coded, needs to change
		rtRefRayPow = rtSysRadPow/real(rtMCNumRays,8)
!		rtBeta = rtKappa + rtSigma
!		rtAbsThr = rtKappa/rtBeta
!		call traceFromSurfLED(pRatio)
		call traceFromSurfLED()
		call binScreenPolar(20)
		call setMeshNodalValues(rtNodalSrc,"S")
		open(963,file="../data/lmodNodalSrc.dat")
		write(963,'(f20.12)') rtNodalSrc
		close(963)
		rtNodalSrc = 0.d0
		femBCsKnown = .true.
		call runFem()
	end subroutine runCoupledModel

	subroutine readMeshDataFile(fMeshDat)
		integer,parameter :: nMeshDat=101
		character(32),intent(in) :: fMeshDat

		open(nMeshDat,file=fMeshDat)
		read(nMeshDat,*)
		read(nMeshDat,*) meshFilePre
		read(nMeshDat,*)
		read(nMeshDat,*) meshNBins
		read(nMeshDat,*)
		read(nMeshDat,*) meshNumNodes
		read(nMeshDat,*)
		read(nMeshDat,*) meshNumElems
		read(nMeshDat,*)
		read(nMeshDat,*) meshNumDoms
		read(nMeshDat,*)
		read(nMeshDat,*) meshNumSurfs
		read(nMeshDat,*)
		read(nMeshDat,*) meshNumBys
		if(.not.(allocated(meshBys))) then
			allocate(meshBys(meshNumBys))
		end if
		read(nMeshDat,*)
		read(nMeshDat,*) meshBys
		read(nMeshDat,*)
		read(nMeshDat,*) meshNumIntfcs
		if((meshNumBys+meshNumIntfcs) .ne. meshNumSurfs) then
			write(*,'(a)')"Error in specifying the boundaries."
			write(*,'(a)')"Number of boundaries specified differ &
			& from number of surfaces in mesh."
			stop
		end if
		if(meshNumIntfcs .gt. 0) then
			if(.not.(allocated(meshIntfcs))) then
				allocate(meshIntfcs(meshNumIntfcs))
			end if
			read(nMeshDat,*)
			read(nMeshDat,*) meshIntfcs
		else
			call skipReadLines(nMeshDat,2)
		end if
		close(nMeshDat)
	end subroutine readMeshDataFile

	subroutine readRtDataFile(fRtDat)
		integer,parameter :: nRtDat=102
		integer :: i,rtNumDoms
		character(32),intent(in) :: fRtDat

		open(nRtDat,file=fRtDat)
		call skipReadLines(nRtDat,2)
		read(nRtDat,*) rtMCNumRays
		read(nRtDat,*)
		read(nRtDat,*) rtSysRadPow
		read(nRtDat,*)
		read(nRtDat,*) rtLimOutPts
		read(nRtDat,*)
		read(nRtDat,*) rtLimReEmDrops
		read(nRtDat,*)
		read(nRtDat,*) rtSctThr
		read(nRtDat,*)
		read(nRtDat,*) rtfResPre
		read(nRtDat,*)
		read(nRtDat,*) rtLoggerMode

		call skipReadLines(nRtDat,3)
		read(nRtDat,*) rtNumDoms
		if(rtNumDoms .ne. meshNumDoms) then
			write(*,'(a)')"Domain info in RT data incorrect."
			write(*,'(a)')"Check RT data file."
			stop
		end if
		read(nRtDat,*)
		read(nRtDat,*) rtNumLambdas

!		if(.not.(allocated(rtKappa))) allocate(rtKappa(meshNumDoms))
!		if(.not.(allocated(rtSigma))) allocate(rtSigma(meshNumDoms))
!		if(.not.(allocated(rtRefrInd))) allocate(rtRefrInd(meshNumDoms))
!		if(.not.(allocated(rtReEmThr))) allocate(rtReEmThr(meshNumDoms))
!		if(.not.(allocated(rtBeta))) allocate(rtBeta(meshNumDoms))
!		if(.not.(allocated(rtAbsThr))) allocate(rtAbsThr(meshNumDoms))

		if(.not.(allocated(rtKappa))) then
			allocate(rtKappa(meshNumDoms,rtNumLambdas))
		end if
		if(.not.(allocated(rtSigma))) then
			allocate(rtSigma(meshNumDoms,rtNumLambdas))
		end if
		if(.not.(allocated(rtBeta))) then
			allocate(rtBeta(meshNumDoms,rtNumLambdas))
		end if
		if(.not.(allocated(rtRefrInd))) then
			allocate(rtRefrInd(meshNumDoms,rtNumLambdas))
		end if
		if(.not.(allocated(rtReEmThr))) then
			allocate(rtReEmThr(meshNumDoms,rtNumLambdas))
		end if
		if(.not.(allocated(rtAnisFac))) then
			allocate(rtAnisFac(meshNumDoms,rtNumLambdas))
		end if
		if(.not.(allocated(rtAbsThr))) then
			allocate(rtAbsThr(meshNumDoms,rtNumLambdas))
		end if		

		read(nRtDat,*)
		do i=1,meshNumDoms
			read(nRtDat,*) rtKappa(i,:)
		end do
		read(nRtDat,*)
		do i=1,meshNumDoms
			read(nRtDat,*) rtSigma(i,:)
		end do
		read(nRtDat,*)
		do i=1,meshNumDoms
			read(nRtDat,*) rtRefrInd(i,:)
		end do
		read(nRtDat,*)
		do i=1,meshNumDoms
			read(nRtDat,*) rtReEmThr(i,:)
		end do
		read(nRtDat,*)
		do i=1,meshNumDoms
			read(nRtDat,*) rtAnisFac(i,:)
		end do

!		Assign values of beta here itself
		do i=1,meshNumDoms
			rtBeta(i,:) = rtKappa(i,:) + rtSigma(i,:)
		end do

!		Assign the absorption threshold here itself
		do i=1,meshNumDoms
			rtAbsThr(i,:) = rtKappa(i,:)/rtBeta(i,:)
		end do

		call skipReadLines(nRtDat,3)
		read(nRtDat,*) rtNumEmSfs
		read(nRtDat,*)
		read(nRtDat,*) rtNumCTEmSfs
		if(rtNumCTEmSfs .gt. 0) then
			allocate(rtCTEmSfIds(rtNumCTEmSfs))
			allocate(rtCTEmSfVals(rtNumCTEmSfs))
			read(nRtDat,*)
			read(nRtDat,*) rtCTEmSfIds
			read(nRtDat,*)
			read(nRtDat,*) rtCTEmSfVals
		else
			call skipReadLines(nRtDat,4)
		end if
		read(nRtDat,*)
		read(nRtDat,*) rtNumCQEmSfs
		if(rtNumCQEmSfs .gt. 0) then
			allocate(rtCQEmSfIds(rtNumCQEmSfs))
			allocate(rtCQEmSfVals(rtNumCQEmSfs))
			read(nRtDat,*)
			read(nRtDat,*) rtCQEmSfIds
			read(nRtDat,*)
			read(nRtDat,*) rtCQEmSfVals
		else
			call skipReadLines(nRtDat,4)
		end if

		call skipReadLines(nRtDat,3)
		read(nRtDat,*) rtNumActBys
		read(nRtDat,*)
		read(nRtDat,*) rtNumBBSfs
		if(rtNumBBSfs .gt. 0) then
			allocate (rtBBSfIds(rtNumBBSfs))
			read(nRtDat,*)
			read(nRtDat,*) rtBBSfIds
		else
			call skipReadLines(nRtDat,2)
		end if
		read(nRtDat,*)
		read(nRtDat,*) rtNumTrSfs
		if(rtNumTrSfs .gt. 0) then
			allocate (rtTrSfIds(rtNumTrSfs))
			read(nRtDat,*)
			read(nRtDat,*) rtTrSfIds
		else
			call skipReadLines(nRtDat,2)
		end if
		read(nRtDat,*)
		read(nRtDat,*) rtNumNPSfs
		if(rtNumNpSfs .gt. 0) then
			allocate (rtNpSfIds(rtNumNpSfs))
			read(nRtDat,*)
			read(nRtDat,*) rtNpSfIds
		else
			call skipReadLines(nRtDat,2)
		end if
		read(nRtDat,*)
		read(nRtDat,*) rtNumParSfs
		read(nRtDat,*)
		read(nRtDat,*) rtNumParDiffSfs
		if(rtNumParDiffSfs .gt. 0) then
			allocate(rtParDiffSfIds(rtNumParDiffSfs))
			allocate(rtParDiffSfVals(rtNumParDiffSfs))
			read(nRtDat,*)
			read(nRtDat,*) rtParDiffSfIds
			read(nRtDat,*)
			read(nRtDat,*) rtParDiffSfVals
		else
			call skipReadLines(nRtDat,4)
		end if
		read(nRtDat,*)
		read(nRtDat,*) rtNumParSpecSfs
		if(rtNumParSpecSfs .gt. 0) then
			allocate(rtParSpecSfIds(rtNumParSpecSfs))
			allocate(rtParSpecSfVals(rtNumParSpecSfs))
			read(nRtDat,*)
			read(nRtDat,*) rtParSpecSfIds
			read(nRtDat,*)
			read(nRtDat,*) rtParSpecSfVals
		else
			call skipReadLines(nRtDat,4)
		end if
		close(nRtDat)
	end subroutine readRtDataFile

	subroutine readFemDataFile(fFemDat)
		integer,parameter :: nFemDat=103
		integer :: i,femNumDoms
		character(32),intent(in) :: fFemDat	

		open(nFemDat,file=fFemDat)
		call skipReadLines(nFemDat,2)
		read(nFemDat,*) femSolverType
		read(nFemDat,*)
		read(nFemDat,*) femSolverMaxIter
		read(nFemDat,*)
		read(nFemDat,*) femResFile
		read(nFemDat,*)
		read(nFemDat,*) femTransient
		if(femTransient) then
			read(nFemDat,*)
			read(nFemDat,*) femTrScheme
		else
			call skipReadLines(nFemDat,2)
		end if

		call skipReadLines(nFemDat,3)

		read(nFemDat,*) femNumDoms
		if(femNumDoms .ne. meshNumDoms) then
			write(*,'(a)')"Domain info in RT data incorrect."
			write(*,'(a)')"Check RT data file."
			stop
		end if
		allocate(femKs(meshNumDoms))
		allocate(femRhos(meshNumDoms))
		allocate(femCaps(meshNumDoms))
		allocate(femDomGens(meshNumDoms))
		allocate(femBCs(meshNumSurfs))
		allocate(femBVals(meshNumSurfs))
		read(nFemDat,*)
        do i=1,meshNumDoms
			read(nFemDat,*)
            read(nFemDat,*) femKs(i)
        end do
		read(nFemDat,*)
		do i=1,meshNumDoms
			read(nFemDat,*)
			read(nFemDat,*) femRhos(i)
		end do
		read(nFemDat,*)
		do i=1,meshNumDoms
			read(nFemDat,*)
			read(nFemDat,*) femCaps(i)
		end do
		read(nFemDat,*)
		do i=1,meshNumDoms
			read(nFemDat,*)
			read(nFemDat,*) femDomGens(i)
		end do
		read(nFemDat,*)
		read(nFemDat,*) femAmbT

		call skipReadLines(nFemDat,3)

        read(nFemDat,*) femBCs
        read(nFemDat,*)
        read(nFemDat,*) femBVals
		close(nFemDat)
	end subroutine readFemDataFile

!	subroutine femSimple()
!		integer,parameter :: femProbFilNum=102
!		integer,allocatable :: mConstTSfIds(:)
!		character(*),parameter :: probFile='../data/femData.dat'
!		character(72) :: mFileName

!		open(femProbFilNum,file=probFile)
!		read(femProbFilNum,*)
!		read(femProbFilNum,*) mFileName
!		read(femProbFilNum,*)
!		read(femProbFilNum,*) femByFile
!		read(femProbFilNum,*)
!		read(femProbFilNum,*) femResFile
!		read(femProbFilNum,*)
!		read(femProbFilNum,*) femTransient
!		read(femProbFilNum,*)
!		read(femProbFilNum,*) femSolverMaxIter
!		read(femProbFilNum,*)
!		read(femProbFilNum,*) femSolverType
!		if(femTransient) then
!			read(femProbFilNum,*)
!			read(femProbFilNum,*) femTrScheme			
!		end if

!		call runFem()
!	end subroutine femSimple

!	subroutine rtSimple()
!		integer,parameter :: rtDatFileNum=102
!		integer :: i,mainCtr,rtIterNum,nConstTSfs,mBinNum(3)
!		integer,allocatable :: mConstTSfIds(:)
!		real(8) :: totEmPow,rayPow,emSurfArea
!		real(8),allocatable :: pRatio(:),mSfConstTs(:)
!		character(*),parameter :: rtProbDatFile='../data/probData.dat'
!		character(72) :: mFileName

!!	Some hard coded values still exist
!		rtIterNum = 120
!		emSurfArea = 4.0d0
!!	Need to rid ourselves of these pesky constants

!		open(rtDatFileNum,file=rtProbDatFile)
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) mFileName
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) mBinNum
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) nConstTSfs
!		allocate(rtEmSfIds(nConstTSfs))
!		allocate(mSfConstTs(nConstTSfs))
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtEmSfIds
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) mSfConstTs
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtMCNumRays
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtKappa
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtSigma
!		close(rtDatFileNum)

!		do mainCtr = 1,rtIterNum
!			if(mainCtr .eq. 1) then
!!				call rtInitMesh(mFileName,mBinNum,mConstTSfIds,			&
!!				mSfConstTs)
!				nConstTSfs = size(mSfConstTs,1)
!				allocate(pRatio(nConstTSfs))
!				totEmPow = sum(mSfConstTs**4.0d0)
!				do i=1,nConstTSfs-1
!					pRatio(i) = sum(mSfConstTs(1:i)**4.0d0)/totEmPow
!				end do
!				pRatio(nConstTSfs) = 1.0d0
!				write(*,*) "pRatio: ", pRatio
!				rtRefRayPow = sigB*totEmPow*emSurfArea/rtMCNumRays
!				write(*,*) "rayPow: ", rayPow
!				call traceFromSurf(pRatio)
!				rtElemSto = rtElemAbs
!				rtWallInf = rtElemAbs
!				rtElemAbs = 0
!			else
!				call traceFromVol()
!				rtElemSto = rtElemAbs + rtWallInf
!				rtElemAbs = 0
!			end if
!			write(*,'(a,2x,i4)')"Main iteration number: ", mainCtr
!			write(*,'(a,2x,i8)')"Absorbed numbers: ", sum(rtElemSto)
!		end do
!	end subroutine rtSimple

!	subroutine rtLED()
!		integer :: i,mainCtr,rtIterNum,mBinNum(3)
!		real(8),allocatable :: pRatio(:),mSfConstTs(:),mSfConstQs(:)
!		character(72) :: mFileName

!		call readRtData(mFileName,mBinNum)
!		rtIterNum = 1
!		do mainCtr = 1,rtIterNum
!			if(mainCtr .eq. 1) then
!				call rtInit(mFileName,mBinNum)
!				allocate(pRatio(1))
!				pRatio = (/1.d0/)
!				rtRefRayPow = 0.355/rtMCNumRays	! Hard coded, needs to change
!				rtBeta = rtKappa + rtSigma
!				rtAbsThr = rtKappa/rtBeta
!				call traceFromSurfLED(pRatio)
!			end if
!		end do	
!	end subroutine rtLED

!	subroutine rtFEMLED()
!		integer :: i,mainCtr,rtIterNum,mBinNum(3)
!		real(8),allocatable :: pRatio(:),mSfConstTs(:),mSfConstQs(:)
!		character(72) :: mFileName

!		call readRtData(mFileName,mBinNum)
!		rtIterNum = 1
!		do mainCtr = 1,rtIterNum
!			if(mainCtr .eq. 1) then
!				call rtInit(mFileName,mBinNum)
!				allocate(pRatio(1))
!				pRatio = (/1.d0/)
!				rtRefRayPow = 0.355/rtMCNumRays	! Hard coded, needs to change
!				rtBeta = rtKappa + rtSigma
!				rtAbsThr = rtKappa/rtBeta
!				call traceFromSurfLED(pRatio)
!			end if
!			call setMeshNodalValues(rtNodalSrc,"S")
!			rtNodalSrc = 0.d0
!			call readFEMData(mFileName)
!			call runFem(mFileName)
!		end do
!	end subroutine rtFEMLED

!	subroutine prepRtGen()
!		integer :: domChk,emSfChk,bySfChk,chkNumSfs,mBinNum(3)
!		character(72) :: mFileName
!		real(8),allocatable :: pRatio(:)

!		call readRtData(mFileName,mBinNum)
!		call rtInit(mFileName,mBinNum)

!		emSfChk = rtNumCTEmSfs+rtNumCQEmSfs
!		if(emSfChk .ne. rtNumEmSfs)	then
!			write(*,'(a)') "Emitting surfaces not correctly specified."
!			write(*,'(a,2x,i2)')"Stated number of emitting surfaces: ",	&
!			rtNumEmSfs
!			write(*,'(a,2x,i2)')"The sum of CT and CQ surfaces: ",		&
!			emSfChk
!		end if
!		bySfChk = rtNumBBSfs+rtNumTrSfs+rtNumNPSfs+rtNumParSfs
!		if(meshNumDoms .eq. 1) then
!			chkNumSfs = meshNumSurfs
!		else
!			if(meshNumDoms .eq. 2) chkNumSfs = meshNumSurfs-2
!			if(meshNumDoms .eq. 3) chkNumSfs = meshNumSurfs-3
!			if(meshNumDoms .eq. 4) chkNumSfs = meshNumSurfs-6
!			if(meshNumDoms .gt. 4) write(*,'(a)') "Surfcheck skipped."
!		end if
!		if(bySfChk .ne. chkNumSfs) then
!			write(*,'(a)') "Boundary surface specification incorrect."
!			write(*,'(a,2x,i2)')"Mesh surface number: ",meshNumSurfs
!			write(*,'(a,2x,i2)')"The sum of all specified boundaries: ",&
!			bySfChk
!		end if

!		allocate(pRatio(1))
!		pRatio = (/1.d0/)
!		rtRefRayPow = 0.355/rtMCNumRays	! Hard coded, needs to change
!		rtBeta = rtKappa + rtSigma
!		rtAbsThr = rtKappa/rtBeta
!		call traceFromSurfLED(pRatio)
!		call setMeshNodalValues(rtNodalSrc,"S")
!		rtNodalSrc = 0.d0
!		call readFEMData(mFileName)
!		call runFem(mFileName)
!		call binScreenPolar(20)
!	end subroutine prepRtGen

!	subroutine readRtSingleDomData(mFileName,mBinNum)
!		integer,parameter :: rtDatFileNum=102
!		integer :: mBinNum(3)
!		character(*),parameter :: rtDatFile='../data/ledRtData.dat'		
!		character(72) :: mFileName

!		open(rtDatFileNum,file=rtDatFile)
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) mFileName
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) mBinNum
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtKappa
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtSigma
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtRefrInd
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtReEmThr
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtMCNumRays

!		call skipReadLines(rtDatFileNum,3)

!		read(rtDatFileNum,*) rtNumEmSfs
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtNumCTEmSfs
!		if(rtNumCTEmSfs .gt. 0) then
!			allocate(rtCTEmSfIds(rtNumCTEmSfs))
!			allocate(rtCTEmSfVals(rtNumCTEmSfs))
!			read(rtDatFileNum,*)
!			read(rtDatFileNum,*) rtCTEmSfIds
!			read(rtDatFileNum,*)
!			read(rtDatFileNum,*) rtCTEmSfVals
!		else
!			call skipReadLines(rtDatFileNum,4)
!		end if
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtNumCQEmSfs
!		if(rtNumCQEmSfs .gt. 0) then
!			allocate(rtCQEmSfIds(rtNumCQEmSfs))
!			allocate(rtCQEmSfVals(rtNumCQEmSfs))
!			read(rtDatFileNum,*)
!			read(rtDatFileNum,*) rtCQEmSfIds
!			read(rtDatFileNum,*)
!			read(rtDatFileNum,*) rtCQEmSfVals
!		else
!			call skipReadLines(rtDatFileNum,4)
!		end if

!		call skipReadLines(rtDatFileNum,3)

!		read(rtDatFileNum,*) rtNumBBSfs
!		if(rtNumBBSfs .gt. 0) then
!			allocate (rtBBSfIds(rtNumBBSfs))
!			read(rtDatFileNum,*)
!			read(rtDatFileNum,*) rtBBSfIds
!		else
!			call skipReadLines(rtDatFileNum,2)
!		end if
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtNumTrSfs
!		if(rtNumTrSfs .gt. 0) then
!			allocate (rtTrSfIds(rtNumTrSfs))
!			read(rtDatFileNum,*)
!			read(rtDatFileNum,*) rtTrSfIds
!		else
!			call skipReadLines(rtDatFileNum,2)
!		end if
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtNumNPSfs
!		if(rtNumNpSfs .gt. 0) then
!			allocate (rtNpSfIds(rtNumNpSfs))
!			read(rtDatFileNum,*)
!			read(rtDatFileNum,*) rtNpSfIds
!		else
!			call skipReadLines(rtDatFileNum,2)
!		end if

!		close(rtDatFileNum)
!	end subroutine readRtSingleDomData

!	subroutine rtFemSimple()
!		integer,parameter :: rtDatFileNum=101,femProbFilNum=102
!		integer :: i,mainCtr,rtIterNum,nConstTSfs,mBinNum(3)
!		integer,allocatable :: mConstTSfIds(:)
!		real(8) :: totEmPow,rayPow,emSurfArea
!		real(8),allocatable :: pRatio(:),mSfConstTs(:),oldSrc(:),		&
!		avSrc(:)
!		character(*),parameter :: rtProbDatFile='../data/probData.dat',	&
!		femProbFile='../data/femData.dat'
!		character(72) :: mFileName,tempFileName,ctrString

!!	Some hard coded values still exist
!		rtIterNum = 80
!		emSurfArea = 4.0d0
!!	Need to rid ourselves of these pesky constants

!!	Read data for the RT and FEM problems from separate files (temporary)
!!	Read RT data
!		open(rtDatFileNum,file=rtProbDatFile)
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) mFileName
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) mBinNum
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) nConstTSfs
!		allocate(rtEmSfIds(nConstTSfs))
!		allocate(mSfConstTs(nConstTSfs))
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtEmSfIds
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) mSfConstTs
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtMCNumRays
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtKappa
!		read(rtDatFileNum,*)
!		read(rtDatFileNum,*) rtSigma
!		close(rtDatFileNum)

!		tempFileName = femResFile

!		do mainCtr = 1,rtIterNum
!			write(ctrString,'(i4.4)') mainCtr
!			femResFile = trim(adjustl(tempFileName))//trim(adjustl(ctrString))
!			if(mainCtr .eq. 1) then
!!				call rtInitMesh(mFileName,mBinNum,mConstTSfIds,			&
!!				mSfConstTs)
!				allocate(oldSrc(meshNumNodes))
!				oldSrc = 0.0d0
!				nConstTSfs = size(mSfConstTs,1)
!				allocate(pRatio(nConstTSfs))
!				totEmPow = sum(mSfConstTs**4.0d0)
!				do i=1,nConstTSfs-1
!					pRatio(i) = sum(mSfConstTs(1:i)**4.0d0)/totEmPow
!				end do
!				pRatio(nConstTSfs) = 1.0d0
!				write(*,*) "pRatio: ", pRatio
!				rtRefRayPow = sigB*totEmPow*emSurfArea/rtMCNumRays
!				write(*,*) "rayPow: ", rtRefRayPow
!				call traceFromSurf(pRatio)
!				rtElemSto = rtElemAbs
!				rtWallInf = rtElemAbs
!				rtElemAbs = 0
!			else
!				call traceVolPowerBased()
!				rtElemSto = rtElemAbs + rtWallInf
!				rtElemAbs = 0
!			end if
!			call getAveragedSource(oldSrc,avSrc)
!			rtNodalSrc = 0.0d0
!			call setMeshNodalValues(avSrc,"S")
!			call runFem()
!			oldSrc = avSrc
!			deallocate(avSrc)
!			write(*,'(a,2x,i4)')"Main iteration number: ", mainCtr
!			write(*,'(a,2x,i8)')"Absorbed numbers: ", sum(rtElemSto)
!		end do

!	end subroutine rtFemSimple

	subroutine getAveragedSource(oldSrc,avSrc)
		real(8),parameter :: expAvFact = 0.05d0
		real(8),intent(in) :: oldSrc(:)
		real(8),allocatable,intent(out) :: avSrc(:)

		allocate(avSrc(meshNumNodes))
		avSrc = rtNodalSrc*expAvFact + oldSrc*(1.0d0-expAvFact)
	end subroutine getAveragedSource

	subroutine binScreenPolar(nBins)
		integer,parameter :: nScrPts=151,nScrBins=152,lcYellow=6,		&
		lcBlue=3
		integer,intent(in) :: nBins
		integer :: i,j,nLines,error,lambda,binTh,binPh
		integer,allocatable :: bY(:,:),bB(:,:)
		real(8) :: th,ph,stepTh,stepPh,tempD1,tempD2,pt(3),dir(3)
		real(8),allocatable :: thetas(:),phis(:)
!		character(*),parameter :: fScrPts=commResDir//"scrPts.out",	&
!								  fScrBins=commResDir//"scrPolBins.out"
		character :: fScrPts*72,fScrBins*72

		fScrPts=commResDir//trim(adjustl(rtfResPre))//"_scrPts.out"
		fScrBins=commResDir//trim(adjustl(rtfResPre))//"_scrPolBins.res"
		open(nScrPts,file=fScrPts)
		call getFileNumLines(nScrPts,nLines)
		close(nScrPts)

		write(*,*) "nTotalPoints: ", nLines
		stepTh = pi/(2.d0*real(nbins,8))
		stepPh = 2.d0*pi/real(nbins,8)
		allocate(thetas(nBins+1))
		allocate(phis(nBins+1))
		thetas(1) = 0.d0
		phis(1) = 0.d0
		do i=1,nBins
			thetas(i+1) = i*stepTh
			phis(i+1) = i*stepPh
		end do
		allocate(bY(nBins,nBins))
		allocate(bB(nBins,nBins))
		bY = 0
		bB = 0
		open(nScrPts,file=fScrPts)
		do i=1,nLines
			read(nScrPts,'(3(e16.9,2x),i2)') pt,lambda
			! Remember to change next line to use an arbitrary origin
			dir = pt/(norm2(pt))
			! Points going backward have to be ignored 
			! (although there shouldn't be any)
			if(dir(3).lt. 0.d0) cycle
			th = acos(dir(3))
			tempD1 = dir(1)/sin(th)
			tempD2 = dir(2)/sin(th)
			if((tempD1.gt. 0.d0).and.(tempD2.gt. 0.d0)) then
				ph = acos(tempD1)
			elseif((tempD1.lt. 0.d0).and.(tempD2.gt. 0.d0)) then
				ph = acos(tempD1)
			elseif((tempD1.lt. 0.d0).and.(tempD2.lt. 0.d0)) then
				ph = atan(tempD2/tempD1) + pi
			elseif((tempD1.gt. 0.d0).and.(tempD2.lt. 0.d0)) then
				ph = asin(tempD2) + 2.d0*pi
			else
				cycle
			end if
			do j=1,nBins
				if((th.gt.thetas(j)).and.(th.lt.thetas(j+1)))then
					binTh = j
				end if
				if((ph.gt.phis(j)).and.(ph.lt.phis(j+1))) then
					binPh = j
				end if
			end do
			if(lambda .eq. lcYellow) then
				bY(binTh,binPh) = bY(binTh,binPh) + 1
			else
				bB(binTh,binPh) = bB(binTh,binPh) + 1
			end if
		end do
		close(nScrPts)

		open(nScrBins,file=fScrBins)
		do i=1,nBins
			write(nScrBins,'(<nBins>(i8,2x))') bY(i,:)
			write(nScrBins,'(<nBins>(i8,2x))') bB(i,:)
		end do
		close(nScrBins)

	end subroutine binScreenPolar

	subroutine binWallPoints(nBins,scrMin,scrMax)
		integer,parameter :: nScrIncs=151,nScrBins=152,lcYellow=6,		&
		lcBlue=3
		integer,intent(in) :: nBins
		integer :: i,nLines,error,lambda,bLocs(2)
		integer,allocatable :: bY(:,:),bB(:,:)
		real(8) :: pt(3),edges(3)
		real(8),intent(in) :: scrMin(3),scrMax(3)
		character(*),parameter :: fSfEmPts=commResDir//"wallIncs.out",	&
								  fScPtsBin=commResDir//"wallBins.out"

		open(nScrIncs,file=fSfEmPts)
		call getFileNumLines(nScrIncs,nLines)
		write(*,*) "nLines: ", nLines
		edges = scrMax - scrMin
		allocate(bY(nBins,nBins))
		allocate(bB(nBins,nBins))
		bY = 0
		bB = 0
		open(nScrIncs,file=fSfEmPts)
		do i=1,nLines
			read(nScrIncs,'(3(e16.9,2x),i2)') pt,lambda
			bLocs = nint(((pt(1:2)-scrMin(1:2))/edges(1:2))*nBins)
			if((bLocs(1).gt.0).and.(bLocs(2).gt.0)) then
				if((bLocs(1).le.nBins).and.(bLocs(2).le.nBins)) then
					if(lambda .eq. lcYellow) bY(bLocs(1),bLocs(2)) = 	&
					bY(bLocs(1),bLocs(2)) + 1
					if(lambda .eq. lcBlue) bB(bLocs(1),bLocs(2)) = 		&
					bB(bLocs(1),bLocs(2)) + 1
				end if
			end if
		end do
		close(nScrIncs)

		open(nScrBins,file=fScPtsBin)
		do i=1,nBins
			write(nScrBins,'(<nBins>(i8,2x))') bY(i,:)
			write(nScrBins,'(<nBins>(i8,2x))') bB(i,:)
		end do
		close(nScrBins)

	end subroutine binWallPoints

end module problem
