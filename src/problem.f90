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
		call binScreenPolar()
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
		character(72) :: fMeshDat,fRtDat

		fMeshDat = commDatDir//trim(adjustl(fMesh))//commDatExt
		fRtDat = commDatDir//trim(adjustl(fRt))//commDatExt
		call readMeshDataFile(fMeshDat)
		call readRtDataFile(fRtDat)
		call rtInit()
		call traceFromSurfLED()
		call binScreenPolar()
		call getScreenBinFluxes()
		call rtClose()
		rtRepeatMesh = .true.
		call readRtDataFile(fRtDat)
		call rtInit()
		call traceFromSurfLED()
		call binScreenPolar()
		call getScreenBinFluxes()
		call rtClose()
	end subroutine rtSimple

	subroutine simAnnealing(fMesh,fRt)
		integer,parameter :: nSARuns = 400,nBins = 26,nSpecDat = 161,	&
		nRunsDat = 162, nCFs = 163, nParams = 5
		integer :: i,j,k,saCt
		real(8),parameter :: T0=0.1d0,kapL=2500.d0,kapH=10000.d0,		&
		sigL=2000.d0,sigH=12000.d0,gL=-0.99d0,gH=0.99d0,nL=1.5d0,		&
		nH = 1.8d0
		real(8) :: fCost0,fCostn,T1,Tk,vSpaceRad,simRad,dev,rad,		&
		probFunc,selCheck,stoppingCrit,zSel
		real(8),allocatable :: rKnown(:),rRun(:),uniVec(:),newVals(:),	&
		lVals(:),hVals(:),x0(:),xBetter(:),z1(:),z2(:)
		logical :: move,inspace,templogic
		character(*),parameter :: fSpec="ybr",fRuns="simAnnRuns",		&
		fcfs="simAnnCFs"
		character(32),intent(in) :: fMesh,fRt
		character(72) :: fMeshDat,fRtDat,fSpecDat,fRunsDat,fCosts

!		Read the spectrometry results for YBR
		fSpecDat = commDatDir//trim(adjustl(fSpec))//commDatExt
		allocate(rKnown(nBins))
		allocate(rRun(nBins))
		open(nSpecDat,file=fSpecDat)
		do i=1,nBins
			read(nSpecDat,*) rKnown(i)
		end do
		close(nSpecDat)

		fMeshDat = commDatDir//trim(adjustl(fMesh))//commDatExt
		fRtDat = commDatDir//trim(adjustl(fRt))//commDatExt
		call readMeshDataFile(fMeshDat)
		allocate(uniVec(nParams))					! kapBlue,sigBlue,sigYell,g
		allocate(lVals(nParams))
		allocate(hVals(nParams))
		allocate(x0(nParams))
		allocate(xBetter(nParams))
		allocate(z1(nParams))
		allocate(z2(nParams))
		allocate(newVals(nParams))
		lVals = (/kapL,sigL,sigL,gL,nL/)	! Note: 2 sigmas used, one for each lambda
		hVals = (/kapH,sigH,sigH,gH,nH/)	! Only one kappa (assuming kappaY = 10.d0)
		fCosts = commResDir//trim(adjustl(fcfs))//commResExt
		open(nCFs,file=fCosts)
		fRunsDat = commResDir//trim(adjustl(fRuns))//commResExt
		open(nRunsDat,file=fRunsDat)
!-----------------------------------------------------------------------!
!	Random search in space
!-----------------------------------------------------------------------!
		do i=1,nSARuns
			if(i.gt.1) then
				rtRepeatMesh = .true.
			end if
			call readRtDataFile(fRtDat)
			call random_number(uniVec)
			rtKappa(2,1) = uniVec(1)*(hVals(1)-lVals(1)) + lVals(1)
			rtSigma(2,:) = uniVec(2:3)*(hVals(2:3)-lVals(2:3)) + lVals(2:3)
			rtAnisFac(2,:) = uniVec(4)*(hVals(4)-lVals(4)) + lVals(4)
			rtRefrInd(2,:) = uniVec(5)*(hVals(5)-lVals(5)) + lVals(5)
			x0 = (/rtKappa(2,1),rtSigma(2,:),rtAnisFac(2,1),rtRefrInd(2,1)/)
			call rtInit()
			call traceFromSurfLED()
			call binScreenPolar()
			call getScreenBinFluxes()
			rRun = rtYBR
			fCostn = norm2(rRun-rKnown)
			write(nCFs,'(6(e14.6,2x))') fCostn,x0
			call rtClose()
		end do

!-----------------------------------------------------------------------!
!	Simulated annealing algorithm
!-----------------------------------------------------------------------!

!!		Calculate maximum radius in variable space
!		vSpaceRad = norm2(hVals-lVals)
!		write(*,*)"vSpaceRad: ",vSpaceRad
!		dev = 0.2d0*vSpaceRad			! standard deviation, initial step
!		call readRtDataFile(fRtDat)
!		x0 = (/rtKappa(2,1),rtSigma(2,:),rtAnisFac(2,1),rtRefrInd(2,1)/)	! Hard coded for domain 2 right now.
!		xBetter = x0					! At the start, the better function loc is the same
!		call rtInit()
!		call traceFromSurfLED()
!		call binScreenPolar()
!		call getScreenBinFluxes()
!		rRun = rtYBR
!		fCost0 = norm2(rRun-rKnown)
!		call rtClose()
!		write(nRunsDat,'(5(e14.6,2x))') x0
!		write(nRunsDat,*) "Starting Cost Function: ", fCost0
!		write(nCFs,'(6(e14.6,2x))') fCost0,x0
!		stoppingCrit = fCost0*0.01d0
!		T1 = T0 - 0.01d0*T0

!!		Now to start the iterations
!		rtRepeatMesh = .true.
!		do i=1,nSARuns
!			Tk = T1
!			inspace = .false.
!			do while(.not.(inspace))
!				call getBoxMuellerNormals(z1,z2)
!				call random_number(zSel)
!				if(zSel .lt. 0.5d0) then
!					uniVec = z1/norm2(z1)
!					rad = dev*z2(1)
!				else
!					uniVec = z2/norm2(z2)
!					rad = dev*z1(1)
!				end if
!				newVals = x0 + rad*uniVec
!				call valsInSpaceCheck(lVals,hVals,newVals,inspace)
!				if(.not.(inspace)) then
!					write(*,*) "Not in space"
!					write(*,'(6(e14.6,2x))') rad, newVals
!					stop
!				end if
!			end do
!			write(*,'(a)')"New values in simulation: "
!			write(*,'(4(e14.6,2x))') newVals
!			write(nRunsDat,*) "Values in the run number: ", i
!			write(nRunsDat,'(5(e14.6,2x))') newVals
!			call readRtDataFile(fRtDat)
!			rtKappa(2,1) = newVals(1)
!			rtSigma(2,:) = newVals(2:3)
!			rtAnisFac(2,:) = newVals(4)
!			rtRefrInd(2,:) = newVals(5)
!			call rtInit()
!			call traceFromSurfLED()
!			call binScreenPolar()
!			call getScreenBinFluxes()
!			rRun = rtYBR
!			fCostn = norm2(rRun-rKnown)
!			call rtClose()
!			probFunc = exp((fCost0-fCostn)/Tk)
!			call random_number(selCheck)
!			write(nRunsDat,*) "New Cost Function: ", fCostn
!			write(nRunsDat,*) "New Probability Function: ", probFunc
!			write(nRunsDat,*) "Selection check: ",selCheck
!			if(fCostn .lt. fCost0) then
!				x0 = newVals
!				dev = dev*exp(-1.d0/(real(nSARuns/2,8)))
!				fCost0 = fCostn
!				Tk = ((T1/T0)**real(i,8))*Tk
!				write(nRunsDat,*)"New point chosen since fCostn < fCost0."
!				write(nCFs,'(6(e14.6,2x))') fCostn,x0
!			else
!				if(selCheck .lt. probFunc) then
!					xBetter = x0
!					x0 = newVals
!					Tk = ((T1/T0)**real(i,8))*Tk
!					write(nRunsDat,*)"Tk = ",Tk
!					write(nRunsDat,*)"New point chosen from selCheck."
!					write(nCFs,'(6(e14.6,2x))') fCostn,x0
!				else
!					cycle
!				end if
!			end if
!			if(fCostn .lt. stoppingCrit) then
!				write(*,'(a)') "Simulation reached stopping criterion."
!				write(*,'(a)')"New values in simulation: "
!				write(*,'(5(e14.6,2x))') newVals				
!			end if
!		end do
	end subroutine simAnnealing

	subroutine valsInSpaceCheck(lows,highs,inVals,check)
		integer :: i
		real(8),intent(in) :: lows(:),highs(:),inVals(:)
		logical,intent(out) :: check

		do i=1,size(lows,1)
			if((inVals(i).lt.lows(i)).or.(inVals(i).gt.highs(i))) then
				check = .false.
				return
			else
				check = .true.
			end if
		end do
	end subroutine valsInSpaceCheck

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
		read(nMeshDat,*) meshScalFac
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
		read(nRtDat,*) rtNumPolBins
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
		if(.not.(allocated(rtDomColCon))) then
			allocate(rtDomColCon(meshNumDoms,rtNumLambdas))
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
		read(nRtDat,*)
		do i=1,meshNumDoms
			read(nRtDat,*) rtDomColCon(i,:)
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

	subroutine binScreenPolar()
		integer,parameter :: nScrPts=151,nScrBins=152,lcYellow=6,		&
		lcBlue=3
		integer :: i,j,nLines,error,lambda,bTh,bPh
		real(8) :: th,ph,stepTh,stepPh,tempD1,tempD2,pt(3),dir(3)
		real(8),allocatable :: thetas(:),phis(:)
		character :: fScrPts*72,fScrBins*72

		fScrPts=commResDir//trim(adjustl(rtfResPre))//"_scrPts.out"
		fScrBins=commResDir//trim(adjustl(rtfResPre))//"_scrPolBins.res"
		stepTh = pi/(2.d0*real(rtNumPolBins,8))
		stepPh = 2.d0*pi/real(rtNumPolBins,8)
		allocate(thetas(rtNumPolBins+1))
		allocate(phis(rtNumPolBins+1))
		thetas(1) = 0.d0
		phis(1) = 0.d0
		do i=1,rtNumPolBins
			thetas(i+1) = i*stepTh
			phis(i+1) = i*stepPh
		end do
		if(.not.(allocated(rtPolarBins))) then
			allocate(rtPolarBins(rtNumPolBins,rtNumPolBins,rtNumLambdas))
		end if
		rtPolarBins = 0
		open(nScrPts,file=fScrPts)
		do i=1,rtNumScrPts
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
			do j=1,rtNumPolBins
				if((th.gt.thetas(j)).and.(th.lt.thetas(j+1)))then
					bTh = j
				end if
				if((ph.gt.phis(j)).and.(ph.lt.phis(j+1))) then
					bPh = j
				end if
			end do
			if(lambda .eq. lcYellow) then
				rtPolarBins(bTh,bPh,2)=rtPolarBins(bTh,bPh,2) + 1
			else
				rtPolarBins(bTh,bPh,1)=rtPolarBins(bTh,bPh,1) + 1
			end if
		end do
		close(nScrPts)

		open(nScrBins,file=fScrBins)
		do i=1,rtNumPolBins
			if(rtNumLambdas.eq. 2) then
				write(nScrBins,'(<rtNumPolBins>(i8,2x))') rtPolarBins(i,:,2)
				write(nScrBins,'(<rtNumPolBins>(i8,2x))') rtPolarBins(i,:,1)
			else
				write(nScrBins,'(<rtNumPolBins>(i8,2x))') rtPolarBins(i,:,1)
			end if
		end do
		close(nScrBins)

	end subroutine binScreenPolar

	subroutine getScreenAngleFluxes()
		integer,parameter :: nScrPts=155,nScrAngFlx=156,lcYellow=6,		&
		lcBlue=3,nAngles=26	! note that the number of angles is a parameter
		integer :: i,j,lambda,nLines
		real(8) :: th,dTh,stepTh,pt(3),dir(3),thetas(nAngles),			&
		areas(nAngles)
		real(8),allocatable :: ptFluxes(:,:)
		character :: fScrPts*72,fScrAngFluxes*72

		fScrPts=commResDir//trim(adjustl(rtfResPre))//"_scrPts.out"
		fScrAngFluxes=commResDir//trim(adjustl(rtfResPre))//"_scrPtFluxes.res"

		allocate(ptFluxes(nAngles,2))
		ptFluxes = 0
		stepTh = pi/(2.d0*(nAngles-1))
		dTh = pi/180.d0
		do j=1,nAngles
			thetas(j) = (j-1)*stepTh
			if((j.gt.1).and.(j.lt.nAngles)) then
				areas(j) = 2.d0*pi*(cos(thetas(j)-dTh) - 				&
				cos(thetas(j)+dTh))
			else
				if(j.eq.1) then
					areas(j) = 2.d0*pi*(cos(thetas(j)) - 				&
					cos(thetas(j)+dTh))
				elseif(j.eq.nAngles) then
					areas(j) = 2.d0*pi*(cos(thetas(j)-dTh) -			&
					cos(thetas(j)))
				end if
			end if
		end do
		open(nScrPts,file=fScrPts)
		do i=1,rtNumScrPts
			read(nScrPts,'(3(e16.9,2x),i2)') pt,lambda
			! Remember to change next line to use an arbitrary origin
			dir = pt/(norm2(pt))
			! Points going backward have to be ignored 
			! (although there shouldn't be any)
			if(dir(3).lt. 0.d0) cycle
			! otherwise, theta is just the angle subtended
			th = acos(dir(3))
			do j=1,nAngles
				if(abs(th-thetas(j)).lt.dTh) then
					if(lambda.eq.lcBlue) then
						ptFluxes(j,1) = ptFluxes(j,1) + rtRefRayPow
					elseif(lambda.eq.lcYellow) then
						ptFluxes(j,2) = ptFluxes(j,2) + rtRefRayPow
					else
						write(*,'(a)') "Point colour not recognised."
					end if
					exit
				end if
			end do
		end do
		close(nScrPts)

		ptFluxes(:,1) = ptFluxes(:,1)/areas
		ptFluxes(:,2) = ptFluxes(:,2)/areas
		open(nScrAngFlx,file=fScrAngFluxes)
		do i=1,nAngles
			write(nScrAngFlx,'(5(e16.9,2x))') ptFluxes(i,:),			&
			sum(ptFluxes(i,:)),ptFluxes(i,:)/sum(ptFluxes(i,:))
		end do
		close(nScrAngFlx)
	end subroutine getScreenAngleFluxes

	subroutine getScreenBinFluxes()
		integer,parameter :: nScrFluxes=153
		integer :: i,j,k
		real(8) :: th,ph,stepTh
		real(8),allocatable :: thetas(:),areas(:)
		character :: fScrFluxes*72

		fScrFluxes=commResDir//trim(adjustl(rtfResPre))//"_scrFlx.out"
		if(.not.(allocated(rtScrFluxes))) then
			allocate(rtScrFluxes(rtNumPolBins,rtNumLambdas))
		end if
		allocate(thetas(rtNumPolBins+1))
		allocate(areas(rtNumPolBins))
		thetas = 0.d0;
		stepTh = pi/(2.d0*real(rtNumPolBins,8))
		
		do i=1,rtNumPolBins
			thetas(i+1) = i*stepTh
			areas(i) = 2*pi*(cos(thetas(i))-cos(thetas(i+1)))
			if(rtNumLambdas.eq. 2) then
				rtScrFluxes(i,1) = rtRefRayPow*sum(rtPolarBins(i,:,1))/areas(i)
				rtScrFluxes(i,2) = rtLamPowRatio*rtRefRayPow*sum(rtPolarBins(i,:,2))/areas(i)
			else
				rtScrFluxes(i,1) = rtRefRayPow*sum(rtPolarBins(i,:,1))/areas(i)
			end if
		end do
		
		open(nScrFluxes,file=fScrFluxes)
		do i=1,rtNumPolBins
			write(nScrFluxes,'(5(e16.9,2x))') rtScrFluxes(i,:),			&
			sum(rtScrFluxes(i,:)),rtScrFluxes(i,:)/sum(rtScrFluxes(i,:))
		end do
		close(nScrFluxes)
		call getAngleFluxes()
	end subroutine getScreenBinFluxes

	subroutine getAngleFluxes()
		integer,parameter :: nMatchedFluxes=154,nMatchBins=26
		integer :: i,j,k
		real(8) :: adjFluxes(nMatchBins,2)
		character :: fAdjFluxes*72

		fAdjFluxes=commResDir//trim(adjustl(rtfResPre))//"_scrAdjFlx.res"
		j = 1
		do i=2,rtNumPolBins-2,2
			j = j+1
			adjFluxes(j,1) = (rtScrFluxes(i,1) + rtScrFluxes(i+1,1))/2.d0
			adjFluxes(j,2) = (rtScrFluxes(i,2) + rtScrFluxes(i+1,2))/2.d0
		end do
		adjFluxes(1,:) = rtScrFluxes(1,:)
		adjFluxes(nMatchBins,:) = rtScrFluxes(rtNumPolBins,:)
		open(nMatchedFluxes,file=fAdjFluxes)
		do i=1,nMatchBins
			write(nMatchedFluxes,'(6(e16.9,2x))') adjFluxes(i,:),		&
			sum(adjFluxes(i,:)),adjFluxes(i,:)/sum(adjFluxes(i,:)),		&
			adjFluxes(i,2)/adjFluxes(i,1)
		end do
		close(nMatchedFluxes)
		if(.not.(allocated(rtYBR))) then
			allocate(rtYBR(nMatchBins))
		end if
		if(rtNumLambdas.eq. 2) then
			rtYBR = adjFluxes(:,2)/adjFluxes(:,1)
		end if
	end subroutine getAngleFluxes

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

!		open(nScrIncs,file=fSfEmPts)
!		call getFileNumLines(nScrIncs,nLines)
!		write(*,*) "nLines: ", nLines
		write(*,'(a,2x,i8)')"Number of screen points: ", rtNumScrPts
		write(*,'(a,2x,i8)')"Number of bottom points: ", rtNumBotPts
		edges = scrMax - scrMin
		allocate(bY(nBins,nBins))
		allocate(bB(nBins,nBins))
		bY = 0
		bB = 0
		open(nScrIncs,file=fSfEmPts)
		do i=1,rtNumScrPts
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
