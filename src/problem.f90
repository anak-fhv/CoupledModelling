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
		allocate(rtSfPowRatio(rtNumEmSfs))
		rtSfPowRatio = (/1.d0/)
		rtRefRayPow = rtSysRadPow/real(rtMCNumRays,8)
		call traceFromSurfLED()
		call binScreenPolar(20)
		call setMeshNodalValues(rtNodalSrc,"S")
		rtNodalSrc = 0.d0
		femBCsKnown = .true.
		call runFem()
	end subroutine runCoupledModel

	subroutine femSimple(fMesh,fFem)
		character(32),intent(in) :: fMesh,fFem
		character(32) :: fMeshDat,fFemDat

		fMeshDat = commDatDir//trim(adjustl(fMesh))//commDatExt
		fFemDat = commDatDir//trim(adjustl(fFem))//commDatExt
		
		call readMeshDataFile(fMeshDat)
		call readFemDataFile(fFemDat)
		femBCsKnown = .true.
		call runFem()
	end subroutine femSimple

	subroutine rtSimple(fMesh,fRt)
		character(32),intent(in) :: fMesh,fRt
		character(32) :: fMeshDat,fRtDat

		fMeshDat = commDatDir//trim(adjustl(fMesh))//commDatExt
		fRtDat = commDatDir//trim(adjustl(fRt))//commDatExt
		call readMeshDataFile(fMeshDat)
		call readRtDataFile(fRtDat)
		call rtInit()
!		allocate(rtSfPowRatio(rtNumEmSfs))
!		rtSfPowRatio = (/1.d0/)
		rtRefRayPow = rtSysRadPow/real(rtMCNumRays,8)
		call traceFromSurfLED()
		call binScreenPolar(20)
	end subroutine rtSimple

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
			allocate(rtCTEmSfTemps(rtNumCTEmSfs))
			read(nRtDat,*)
			read(nRtDat,*) rtCTEmSfIds
			read(nRtDat,*)
			read(nRtDat,*) rtCTEmSfTemps
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
