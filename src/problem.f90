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
		read(rtDatFileNum,*) rtMCNumRays
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtKappa
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtSigma
		close(rtDatFileNum)

		do mainCtr = 1,rtIterNum
			if(mainCtr .eq. 1) then
!				call rtInitMesh(mFileName,mBinNum,mConstTSfIds,			&
!				mSfConstTs)
				nConstTSfs = size(mSfConstTs,1)
				allocate(pRatio(nConstTSfs))
				totEmPow = sum(mSfConstTs**4.0d0)
				do i=1,nConstTSfs-1
					pRatio(i) = sum(mSfConstTs(1:i)**4.0d0)/totEmPow
				end do
				pRatio(nConstTSfs) = 1.0d0
				write(*,*) "pRatio: ", pRatio
				rtRefRayPow = sigB*totEmPow*emSurfArea/rtMCNumRays
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


	subroutine rtLED()
		integer :: i,mainCtr,rtIterNum,mBinNum(3)
		real(8),allocatable :: pRatio(:),mSfConstTs(:),mSfConstQs(:)
		character(72) :: mFileName

		call readRtData(mFileName,mBinNum)
		rtIterNum = 1
		do mainCtr = 1,rtIterNum
			if(mainCtr .eq. 1) then
				call rtInit(mFileName,mBinNum)
				allocate(pRatio(1))
				pRatio = (/1.d0/)
				rtRefRayPow = 0.355/rtMCNumRays	! Hard coded, needs to change
				rtBeta = rtKappa + rtSigma
				rtAbsThr = rtKappa/rtBeta
				call traceFromSurfLED(pRatio)
			end if
		end do	
	end subroutine rtLED

	subroutine rtFEMLED()
		integer :: i,mainCtr,rtIterNum,mBinNum(3)
		real(8),allocatable :: pRatio(:),mSfConstTs(:),mSfConstQs(:)
		character(72) :: mFileName

		call readRtData(mFileName,mBinNum)
		rtIterNum = 1
		do mainCtr = 1,rtIterNum
			if(mainCtr .eq. 1) then
				call rtInit(mFileName,mBinNum)
				allocate(pRatio(1))
				pRatio = (/1.d0/)
				rtRefRayPow = 0.355/rtMCNumRays	! Hard coded, needs to change
				rtBeta = rtKappa + rtSigma
				rtAbsThr = rtKappa/rtBeta
				call traceFromSurfLED(pRatio)
			end if
			call setMeshNodalValues(rtNodalSrc,"S")
			rtNodalSrc = 0.d0
			call readFEMData(mFileName)
			call runFem(mFileName)
		end do
	end subroutine rtFEMLED

	subroutine prepRtGen()
		integer :: domChk,emSfChk,bySfChk,chkNumSfs,mBinNum(3)
		character(72) :: mFileName
		real(8),allocatable :: pRatio(:)

		call readRtData(mFileName,mBinNum)
		call rtInit(mFileName,mBinNum)

		emSfChk = rtNumCTEmSfs+rtNumCQEmSfs
		if(emSfChk .ne. rtNumEmSfs)	then
			write(*,'(a)') "Emitting surfaces not correctly specified."
			write(*,'(a,2x,i2)')"Stated number of emitting surfaces: ",	&
			rtNumEmSfs
			write(*,'(a,2x,i2)')"The sum of CT and CQ surfaces: ",		&
			emSfChk
		end if
		bySfChk = rtNumBBSfs+rtNumTrSfs+rtNumNPSfs+rtNumParSfs
		if(meshNumDoms .eq. 1) then
			chkNumSfs = meshNumSurfs
		else
			if(meshNumDoms .eq. 2) chkNumSfs = meshNumSurfs-2
			if(meshNumDoms .eq. 3) chkNumSfs = meshNumSurfs-3
			if(meshNumDoms .eq. 4) chkNumSfs = meshNumSurfs-6
			if(meshNumDoms .gt. 4) write(*,'(a)') "Surfcheck skipped."
		end if
		if(bySfChk .ne. chkNumSfs) then
			write(*,'(a)') "Boundary surface specification incorrect."
			write(*,'(a,2x,i2)')"Mesh surface number: ",meshNumSurfs
			write(*,'(a,2x,i2)')"The sum of all specified boundaries: ",&
			bySfChk
		end if

		allocate(pRatio(1))
		pRatio = (/1.d0/)
		rtRefRayPow = 0.355/rtMCNumRays	! Hard coded, needs to change
		rtBeta = rtKappa + rtSigma
		rtAbsThr = rtKappa/rtBeta
		call traceFromSurfLED(pRatio)
		call setMeshNodalValues(rtNodalSrc,"S")
		rtNodalSrc = 0.d0
		call readFEMData(mFileName)
		call runFem(mFileName)
	end subroutine prepRtGen

	subroutine readRtData(mFileName,mBinNum)
		integer,parameter :: rtDatFileNum=102
		integer :: mBinNum(3)
		character(*),parameter :: rtDatFile='../data/gridLedRtData.dat'		
		character(72) :: mFileName

		open(rtDatFileNum,file=rtDatFile)
		call skipReadLines(rtDatFileNum,2)
		read(rtDatFileNum,*) mFileName
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) mBinNum
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtMCNumRays
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtLimOutPts
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtLimReEmDrops
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtfResPre
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtLoggerMode

		call skipReadLines(rtDatFileNum,3)
		read(rtDatFileNum,*) meshNumDoms
		if(.not.(allocated(rtKappa))) allocate(rtKappa(meshNumDoms))
		if(.not.(allocated(rtSigma))) allocate(rtSigma(meshNumDoms))
		if(.not.(allocated(rtRefrInd))) allocate(rtRefrInd(meshNumDoms))
		if(.not.(allocated(rtReEmThr))) allocate(rtReEmThr(meshNumDoms))
		if(.not.(allocated(rtBeta))) allocate(rtBeta(meshNumDoms))
		if(.not.(allocated(rtAbsThr))) allocate(rtAbsThr(meshNumDoms))
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtKappa
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtSigma
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtRefrInd
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtReEmThr

		call skipReadLines(rtDatFileNum,3)
		read(rtDatFileNum,*) rtNumEmSfs
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumCTEmSfs
		if(rtNumCTEmSfs .gt. 0) then
			allocate(rtCTEmSfIds(rtNumCTEmSfs))
			allocate(rtCTEmSfVals(rtNumCTEmSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtCTEmSfIds
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtCTEmSfVals
		else
			call skipReadLines(rtDatFileNum,4)
		end if
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumCQEmSfs
		if(rtNumCQEmSfs .gt. 0) then
			allocate(rtCQEmSfIds(rtNumCQEmSfs))
			allocate(rtCQEmSfVals(rtNumCQEmSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtCQEmSfIds
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtCQEmSfVals
		else
			call skipReadLines(rtDatFileNum,4)
		end if

		call skipReadLines(rtDatFileNum,3)
		read(rtDatFileNum,*) rtNumBBSfs
		if(rtNumBBSfs .gt. 0) then
			allocate (rtBBSfIds(rtNumBBSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtBBSfIds
		else
			call skipReadLines(rtDatFileNum,2)
		end if
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumTrSfs
		if(rtNumTrSfs .gt. 0) then
			allocate (rtTrSfIds(rtNumTrSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtTrSfIds
		else
			call skipReadLines(rtDatFileNum,2)
		end if
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumNPSfs
		if(rtNumNpSfs .gt. 0) then
			allocate (rtNpSfIds(rtNumNpSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtNpSfIds
		else
			call skipReadLines(rtDatFileNum,2)
		end if
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumParSfs
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumParDiffSfs
		if(rtNumParDiffSfs .gt. 0) then
			allocate(rtParDiffSfIds(rtNumParDiffSfs))
			allocate(rtParDiffSfVals(rtNumParDiffSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtParDiffSfIds
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtParDiffSfVals
		else
			call skipReadLines(rtDatFileNum,4)
		end if
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumParSpecSfs
		if(rtNumParSpecSfs .gt. 0) then
			allocate(rtParSpecSfIds(rtNumParSpecSfs))
			allocate(rtParSpecSfVals(rtNumParSpecSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtParSpecSfIds
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtParSpecSfVals
		else
			call skipReadLines(rtDatFileNum,4)
		end if
		close(rtDatFileNum)

	end subroutine readRtData

	subroutine readRtSingleDomData(mFileName,mBinNum)
		integer,parameter :: rtDatFileNum=102
		integer :: mBinNum(3)
		character(*),parameter :: rtDatFile='../data/ledRtData.dat'		
		character(72) :: mFileName

		open(rtDatFileNum,file=rtDatFile)
		read(rtDatFileNum,*)
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) mFileName
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) mBinNum
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtKappa
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtSigma
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtRefrInd
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtReEmThr
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtMCNumRays

		call skipReadLines(rtDatFileNum,3)

		read(rtDatFileNum,*) rtNumEmSfs
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumCTEmSfs
		if(rtNumCTEmSfs .gt. 0) then
			allocate(rtCTEmSfIds(rtNumCTEmSfs))
			allocate(rtCTEmSfVals(rtNumCTEmSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtCTEmSfIds
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtCTEmSfVals
		else
			call skipReadLines(rtDatFileNum,4)
		end if
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumCQEmSfs
		if(rtNumCQEmSfs .gt. 0) then
			allocate(rtCQEmSfIds(rtNumCQEmSfs))
			allocate(rtCQEmSfVals(rtNumCQEmSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtCQEmSfIds
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtCQEmSfVals
		else
			call skipReadLines(rtDatFileNum,4)
		end if

		call skipReadLines(rtDatFileNum,3)

		read(rtDatFileNum,*) rtNumBBSfs
		if(rtNumBBSfs .gt. 0) then
			allocate (rtBBSfIds(rtNumBBSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtBBSfIds
		else
			call skipReadLines(rtDatFileNum,2)
		end if
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumTrSfs
		if(rtNumTrSfs .gt. 0) then
			allocate (rtTrSfIds(rtNumTrSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtTrSfIds
		else
			call skipReadLines(rtDatFileNum,2)
		end if
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtNumNPSfs
		if(rtNumNpSfs .gt. 0) then
			allocate (rtNpSfIds(rtNumNpSfs))
			read(rtDatFileNum,*)
			read(rtDatFileNum,*) rtNpSfIds
		else
			call skipReadLines(rtDatFileNum,2)
		end if

		close(rtDatFileNum)
	end subroutine readRtSingleDomData

	subroutine readFEMData(mFileName)
		integer,parameter :: femDatFileNum=103
		character(*),parameter :: femDatFile='../data/gridLedFemData.dat'		
		character(72) :: mFileName

		open(femDatFileNum,file=femDatFile)
		read(femDatFileNum,*)
		read(femDatFileNum,*) mFileName
		read(femDatFileNum,*)
		read(femDatFileNum,*) femByFile
		read(femDatFileNum,*)
		read(femDatFileNum,*) femResFile
		read(femDatFileNum,*)
		read(femDatFileNum,*) femTransient
		read(femDatFileNum,*)
		read(femDatFileNum,*) femSolverMaxIter
		read(femDatFileNum,*)
		read(femDatFileNum,*) femSolverType
		if(femTransient) then
			read(femDatFileNum,*)
			read(femDatFileNum,*) femTrScheme			
		end if

	end subroutine readFEMData

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
		read(rtDatFileNum,*) rtMCNumRays
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtKappa
		read(rtDatFileNum,*)
		read(rtDatFileNum,*) rtSigma
		close(rtDatFileNum)

		tempFileName = femResFile

		do mainCtr = 1,rtIterNum
			write(ctrString,'(i4.4)') mainCtr
			femResFile = trim(adjustl(tempFileName))//trim(adjustl(ctrString))
			if(mainCtr .eq. 1) then
!				call rtInitMesh(mFileName,mBinNum,mConstTSfIds,			&
!				mSfConstTs)
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
				rtRefRayPow = sigB*totEmPow*emSurfArea/rtMCNumRays
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

		fScrPts=commResDir//trim(adjustl(rtfResPre))//"scrPts.out"
		fScrBins=commResDir//trim(adjustl(rtfResPre))//"scrPolBins.out"
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
