module rt

	use utilities
	use mesh

	implicit none

!	Declare all relevant types
	type emissionSurface
		real(8) :: totEmPow
		real(8),allocatable :: cuSumFcEmPow(:)
		type(surface) :: emSurf
	end type emissionSurface

!	Module variables
	integer,parameter :: rtNumPFTable=10000,rtElemMinRays=5,			&
	rtSctLrgSphThr = 1, rtSctMieThr = 2
	integer :: rtMCNumRays,rtNumEmSfs,rtNumCTEmSfs,rtNumCQEmSfs,		&
	rtNumBBSfs,rtNumTrSfs,rtNumNpSfs,rtNumParSfs,rtLimOutPts,			&
	rtLimReEmDrops,rtNtraces,rtBBAbsNum,rtTrExitNum,rtNumParDiffSfs,	&
	rtNumParSpecSfs,rtParSfAbsNum,rtAbsNoReEmCt,rtInVolAbsCt,			&
	rtTransBetnDomsCt,rtLoggerMode,rtNumActBys,rtCurrDom,rtNumLambdas,	&
	rtCurrLambda,rtSctThr
	integer,allocatable :: rtEmSfIds(:),rtCTEmSfIds(:),rtCQEmSfIds(:),	&
	rtBBSfIds(:),rtTrSfIds(:),rtNpSfIds(:),rtElemAbs(:),rtElemSto(:),	&
	rtWallInf(:),rtParSfIds(:),rtParDiffSfIds(:),rtParSpecSfIds(:),		&
	rtActBys(:)
	real(8) :: rtRefRayPow,rtSysRadPow
	real(8),allocatable :: rtWallSrc(:),rtNodalSrc(:),rtPFTable(:),		&
	rtCTEmSfVals(:),rtCQEmSfVals(:),rtParDiffSfVals(:),					&
	rtParSpecSfVals(:),rtEmSpectr(:,:),rtReEmSpectr(:,:),rtBeta(:,:),	&
	rtKappa(:,:),rtSigma(:,:),rtRefrInd(:,:),rtReEmThr(:,:),			&
	rtAbsThr(:,:),rtAnisFac(:,:),rtSfPowRatio(:)
	character :: rtfResPre*32,rtEmSpFile*72,rtReEmSpFile*72
	type(emissionSurface),allocatable :: rtEmSurfs(:)

	contains

	subroutine rtInit()
		integer :: i
		call rtInitMesh()
		if(.not.(allocated(rtEmSfIds))) then
			allocate(rtEmSfIds(rtNumEmSfs))
		end if
		if(rtNumParSfs .gt. 0) then
			if(.not.(allocated(rtParSfIds))) then
				allocate(rtParSfIds(rtNumParSfs))
			end if
		end if
		if(.not.(allocated(rtActBys))) then
			allocate(rtActBys(rtNumActBys))
		end if
		rtEmSfIds = (/rtCTEmSfIds,rtCQEmSfIds/)
		rtParSfIds = (/rtParDiffSfIds,rtParSpecSfIds/)
		rtActBys = (/rtBBSfIds,rtTrSfIds,rtNpSfIds,rtParSfIds/)
		call popLargeSphDiffMuTable()
		if(rtNumCTEmSfs .gt. 0) then
			do i=1,rtNumCTEmSfs
				call setSurfaceConstTemperature(rtCTEmSfIds(i),			&
				rtCTEmSfVals(i))
			end do
		end if
		if(.not.(allocated(rtElemAbs))) then
			allocate(rtElemAbs(meshNumElems))
		end if
		if(.not.(allocated(rtElemSto))) then
			allocate(rtElemSto(meshNumElems))
		end if
		if(.not.(allocated(rtWallInf))) then
			allocate(rtWallInf(meshNumElems))
		end if
		if(.not.(allocated(rtWallSrc))) then
			allocate(rtWallSrc(meshNumNodes))
		end if
		if(.not.(allocated(rtNodalSrc))) then
			allocate(rtNodalSrc(meshNumNodes))
		end if
		rtElemAbs = 0
		rtElemSto = 0
		rtWallInf = 0
		rtWallSrc = 0.0d0
		rtNodalSrc = 0.0d0
		rtParSfAbsNum = 0
		rtNtraces = 0
		rtTrExitNum = 0
		rtBBAbsNum = 0
		rtAbsNoReEmCt = 0
		rtInVolAbsCt = 0
		rtTransBetnDomsCt = 0
		write(*,*) "rtLimOutPts: ",rtLimOutPts
		write(*,*) "rtLimReEmDrops: ",rtLimReEmDrops

		call CreateUniformEmissionSurfaces()
	end subroutine rtInit

	subroutine rtInitMesh()
		call readMesh()
		call getElementNeighbours()
		call populateSurfaceFaceAreas()
		call populateElementVolumes()
	end subroutine rtInitMesh

	subroutine handleExit(outPt,outPtCt,endPt,dirOut,lambda,nExitPts,	&
	nScrPts,nBotPts)
		integer,intent(in) :: lambda,nExitPts,nScrPts,nBotPts
		integer,intent(inout) :: outPtCt
		real(8),parameter :: scrDist=1.d0
		real(8) :: bottomSurf,ptScr(3),ptBot(3)
		real(8),intent(in) :: endPt(3),dirOut(3)
		logical,intent(in) :: outPt

		bottomSurf=minval(meshVerts(:,3),1)
		if(outPt) then
			outPtCt = outPtCt + 1
			if(outPtCt .gt. rtLimOutPts) then
				write(*,*)"Count of out points reached limit."
				return
			end if
		else
			write(nExitPts,'(6(e16.9,2x),i2)') endPt,dirOut,lambda
			call getScreenPoint(endPt,dirOut,scrDist,ptScr)
			if(ptScr(3).gt. bottomSurf) then
				write(nScrPts,'(3(e16.9,2x),i2)') ptScr,lambda
			else
				write(nBotPts,'(3(e16.9,2x),i2)') ptScr,lambda
			end if
		end if
	end subroutine handleExit

	subroutine getScreenPoint(p,d,rScr,pScr)
		integer :: i,j
		real(8) :: a,b,c,x(2)
		real(8),intent(in) :: rScr,p(3),d(3)
		real(8),intent(out) :: pScr(3)

!		See derivation: a = s^2(dx^2+dy^2+dz^2), 
!		but dx,dy,dz are direction coordinates
		a = 1.d0
		b = 2.d0*(p(1)*d(1)+p(2)*d(2)+p(3)*d(3))
		c = (p(1)**2.d0+p(2)**2.d0+p(3)**2.d0 - rScr**2.d0)
		call solveQuadratic(a,b,c,x)
		if(x(1) .le. 0.d0) then
			pScr = p + x(2)*d
		else
			pScr = p + x(1)*d
		end if

!		Important to check the point
		if(abs(sum(pScr**2.d0)-rScr**2.d0) .gt. NANO) then
			write(*,'(a)') "Quadratic solved wrong"
		end if

	end subroutine getScreenPoint

	subroutine getNextFace(ec,pt,dir,newFc,lToFc)
		integer :: fc,remNode,fcNodes(3)
		integer,intent(out) :: newFc
		real(8) :: ml,ptToFc,remVert(3),fcNorm(3),cent(3),fcVerts(3,3)
		real(8),intent(in) :: pt(3),dir(3),ec(4,3)
		real(8),intent(out) :: lToFc
		logical :: inward

		inward = .true.
		ml = LARGE
		do fc =1,4
			fcNodes = getFaceNodes(fc)
			remNode = 10-sum(fcNodes)
			fcVerts = ec(fcNodes,:)
			remVert = ec(remNode,:)
			fcNorm = getFaceNorm(fcVerts,remVert,inward)
			ptToFc = dot_product(fcNorm,(pt-fcVerts(1,:)))
			lToFc = ptToFc/dot_product(dir,-fcNorm)
			if(lToFc < 0) cycle
			if (lToFc < ml) then
    			ml = lToFc
    			newFc = fc
			end if
		end do
		lToFc = ml
	end subroutine getNextFace

	subroutine readSpectrumData(specType)
		integer,parameter :: nEmSpecFile=135,nReEmSpecFile=157
		integer :: i,nLines
		character(*),intent(in) :: specType

		if(specType .eq. "emission") then
			open(nEmSpecFile,file=rtEmSpFile)
			call getFileNumLines(nEmSpecFile,nLines)
			if(.not.(allocated(rtEmSpectr))) then
				allocate(rtEmSpectr(nLines,2))
			else
				deallocate(rtEmSpectr)
				allocate(rtEmSpectr(nLines,2))
			end if
			open(nEmSpecFile,file=rtEmSpFile)
			do i=1,nLines
				read(nEmSpecFile,*) rtEmSpectr(i,:)
			end do
		elseif(specType .eq. "reemission") then
			open(nReEmSpecFile,file=rtReEmSpFile)
			call getFileNumLines(nReEmSpecFile,nLines)
			if(.not.(allocated(rtReEmSpectr))) then
				allocate(rtReEmSpectr(nLines,2))
			else
				deallocate(rtReEmSpectr)
				allocate(rtReEmSpectr(nLines,2))
			end if
			open(nReEmSpecFile,file=rtReEmSpFile)
			do i=1,nLines
				read(nReEmSpecFile,*) rtReEmSpectr(i,:)
			end do
		else
			write(*,'(a)')"Wrong selector for spectrum type."
			write(*,'(a)')"Default values will be used."
		end if
	end subroutine readSpectrumData

	subroutine getRayWavelength(newRay,lEm)
		integer :: ind1,ind2
		real(8) :: randLam
		real(8),intent(out) :: lEm
		logical,intent(in) :: newRay

		call random_number(randLam)
		if(newRay) then
			call bisLocReal(rtEmSpectr(:,2),randLam,ind1,ind2)
			lEm = rtEmSpectr(ind2,1)
		else
			call bisLocReal(rtReEmSpectr(:,2),randLam,ind1,ind2)
			lEm = rtReEmSpectr(ind2,1)
		end if
	end subroutine getRayWavelength

	function checkNewPt(pt,dir,ec,newFc) result(inFc)
		integer :: fcNodes(3)
		integer,intent(in) :: newFc
		real(8) :: fcVerts(3,3)
		real(8),intent(in):: dir(3),pt(3),ec(4,3)
		logical :: inFc

		fcNodes = getFaceNodes(newFc)
		fcVerts = ec(fcNodes,:)
		inFc = insideFaceCheck(fcVerts,pt)
	end function checkNewPt

	subroutine startRayInVolume(emEl,pInt,pt,dir)
		integer :: elNodes(4)
		integer,intent(in) :: emEl
		real(8) :: ec(4,3)
		real(8),intent(out) :: pInt,pt(3),dir(3)

		elNodes = meshElems(emEl)%nodes
		ec = meshVerts(elNodes,:)
		pt = selTetraPoint(ec)
		dir = getRaySphDir()
		pInt = getRayPathIntegral()
	end subroutine startRayInVolume

	subroutine startRayFromSurf(emEl,emFc,pInt,pt,dir)
		integer :: remNo,fcNodes(3),elNodes(4)
		integer,intent(out) :: emEl,emFc
		real(8) :: remVert(3),fcNorm(3),ec(4,3),fcVerts(3,3)
!		real(8),intent(in) :: pRatio(:)
		real(8),intent(out) :: pInt,pt(3),dir(3)


		call chooseEmittingFace(emEl,emFc)
! 	Added to the code to make a near point source simulation
! 		For 1040 mu hemisphere
!		emEl = 67
!		emFc = 4

! 		For 1040 mu cylinder
!		emEl = 272
!		emFc = 2

!		For hem.msh
!		emEl = 50865
!		emFc = 1
! 	Added to the code to make a near point source simulation

		elNodes = meshElems(emEl)%nodes
		ec = meshVerts(elNodes,:)
		fcNodes = getFaceNodes(emFc)
		fcVerts = ec(fcNodes,:)
		remNo = 10-sum(fcNodes)
		remVert = ec(remNo,:)
		fcNorm = getFaceNorm(fcVerts,remVert,.true.)
		pt = selTriPoint(fcVerts)
		dir = getFaceDiffRayDir(fcVerts,fcNorm)
		pInt = getRayPathIntegral()
	end subroutine startRayFromSurf

	subroutine chooseEmittingFace(emEl,emFc)
		integer :: i,nSf
		integer,intent(out) :: emEl,emFc
		real(8) :: rs
		type(emissionSurface) :: emSf

		nSf = size(rtSfPowRatio,1)
		call random_number(rs)
		do i=1,nSf
			if(rs .lt. rtSfPowRatio(i)) then
				emSf = rtEmSurfs(i)
				call getEmFace(emSf,emEl,emFc)
				exit
			end if
		end do
	end subroutine chooseEmittingFace

	subroutine getEmFace(emSf,emEl,emFc)
		integer :: i,ind1,ind2,sz
		integer,intent(out) :: emEl,emFc
		real(8):: psi
		type(emissionSurface) :: emSf

		call random_number(psi)
		call bisLocReal(emSf%cuSumFcEmPow,psi,ind1,ind2)
		emEl = emSf%emSurf%elNum(ind2)
		emFc = emSf%emSurf%fcNum(ind2)
	end subroutine getEmFace

	subroutine createUniformEmissionSurfaces()
		integer :: i,j,nEmSf,currSurf,nEmFcs,elNum,fcNum,fcNodes(3),	&
		elNodes(4)
		real(8) :: fcEmPow,fcArea,fcCentT,emSfPow,fcNoTs(3)

		nEmSf = size(rtEmSfIds,1)
		if(.not.(allocated(rtEmSurfs))) then
			allocate(rtEmSurfs(nEmSf))
		end if
		do i=1,nEmSf
			currSurf = rtEmSfIds(i)
			rtEmSurfs(i)%emSurf = meshSurfs(currSurf)
			nEmFcs = meshSurfs(currSurf)%numFcs
			if(.not.(allocated(rtEmSurfs(i)%cuSumFcEmPow))) then
				allocate(rtEmSurfs(i)%cuSumFcEmPow(nEmFcs))
			end if
			emSfPow = 0.d0
			do j=1,nEmFcs
				elNum = rtEmSurfs(i)%emSurf%elNum(j)
				elNodes = meshElems(elNum)%nodes
				fcNum = rtEmSurfs(i)%emSurf%fcNum(j)
				fcArea = rtEmSurfs(i)%emSurf%fcArea(j)
				fcEmPow = fcArea*1.d0					! Temporary to ease notation
				emSfPow = emSfPow+fcEmPow
				rtEmSurfs(i)%cuSumFcEmPow(j) = emSfPow
			end do
			rtEmSurfs(i)%totEmPow = emSfPow
			rtEmSurfs(i)%cuSumFcEmPow=rtEmSurfs(i)%cuSumFcEmPow/emSfPow
		end do		
	end subroutine createUniformEmissionSurfaces

	function scatterRayIsotropic() result(dir)
		real(8) :: rth,rph,th,ph,dir(3)

		call random_number(rth)
		th = acos(1.d0 - 2.d0*rth)
		call random_number(rph)
		ph = 2.d0*pi*rph
		dir = getDirectionCoords(th,ph)
	end function scatterRayIsotropic

	function scatterRayLargeSphereCloud() result(dir)
		real(8) :: rMu,rPhi,mu,th,rph,ph,dir(3)

		call random_number(rMu)
		call getLargeDiffSphMu(rMu,mu)
		th = acos(mu)
		call random_number(rph)
		ph = 2.d0*pi*rph
		dir = getDirectionCoords(th,ph)
	end function scatterRayLargeSphereCloud

	function scatterRayHeyneyGreenStein() result(dir)
		real(8) :: rMu,rPhi,mu,th,rph,ph,dir(3)

		call random_number(rMu)
		call getHeyneyGreensteinMu(rMu,mu)
		th = acos(mu)
		call random_number(rph)
		ph = 2.d0*pi*rph
		dir = getDirectionCoords(th,ph)				
	end function scatterRayHeyneyGreenStein

	function specularReflection(ec,fcNum,dirIn) result(dirOut)
		integer :: remNo,fcNodes(3)
		integer,intent(in) :: fcNum
		real(8) :: cosInc,fcVerts(3,3),remVert(3),fcNorm(3),dirOut(3)
		real(8),intent(in) :: ec(4,3),dirIn(3)

		fcNodes = getFaceNodes(fcNum)
		remNo = 10-sum(fcNodes)
		fcVerts = ec(fcNodes,:)
		remVert = ec(remNo,:)
		fcNorm = getFaceNorm(fcVerts,remVert,.true.)
		cosInc = -dot_product(dirIn,fcNorm)
		if (cosInc .lt. 0.0d0) then
            write(*,*) "change direction of normal"
            cosInc = -cosInc
            fcNorm = -fcNorm
        end if
		dirOut = dirIn + 2.0d0*cosInc*fcNorm
	end function specularReflection

	function diffuseReflection(ec,fcNum) result(dirOut)
		integer,intent(in) :: fcNum
		integer :: remNo,fcNodes(3)
		real(8),intent(in) :: ec(4,3)
		real(8) :: fcVerts(3,3),remVert(3),fcNorm(3),dirOut(3)

		fcNodes = getFaceNodes(fcNum)
		remNo = 10-sum(fcNodes)
		fcVerts = ec(fcNodes,:)
		remVert = ec(remNo,:)
		fcNorm = getFaceNorm(fcVerts,remVert,.true.)
		dirOut = getFaceDiffRayDir(fcVerts,fcNorm)
	end function diffuseReflection

	subroutine transBetweenDoms(ec,fcNum,nhbrDom,dirIn,trans,dirOut)
		integer,intent(in) :: fcNum,nhbrDom
		integer :: remNo,fcNodes(3)
		real(8) :: rIndRatio,th1,th2,cosInc,sinTht,rhoFresnel,frRefChk,	&
		n1,n2,temp1,temp2,remVert(3),fcNorm(3),fcVerts(3,3)
		real(8),intent(in) :: ec(4,3),dirIn(3)
		real(8),intent(out):: dirOut(3)
		logical,intent(out) :: trans

		trans = .false.
		fcNodes = getFaceNodes(fcNum)
		remNo = 10-sum(fcNodes)
		fcVerts = ec(fcNodes,:)
		remVert = ec(remNo,:)
		fcNorm = getFaceNorm(fcVerts,remVert,.true.)
		cosInc = -dot_product(dirIn,fcNorm)
		if (cosInc .lt. 0.0d0) then
            write(*,*) "change direction of normal"
            cosInc = -cosInc
            fcNorm = -fcNorm
        end if
		th1 = acos(cosInc)
		n1 = rtRefrInd(rtCurrDom,rtCurrLambda)
		n2 = rtRefrInd(nhbrDom,rtCurrLambda)
		rIndRatio = n1/n2
		sinTht = rIndRatio*sqrt(1-cosInc**2.d0)
		th2 = asin(sinTht)
		if(sinTht .gt. 1.d0) then
			trans = .false.
			dirOut = specularReflection(ec,fcNum,dirIn)
			return
		else
			call getFresnelReflectivity(th1,th2,rhoFresnel)
			call random_number(frRefChk)
			if(frRefChk .lt. rhoFresnel) then
				trans = .false.
				dirOut = specularReflection(ec,fcNum,dirIn)
				return
			else
				trans = .true.
				dirOut = (rIndRatio*cosInc-sqrt(1.d0-sinTht**2.d0))*fcNorm
				dirOut = dirOut + rIndRatio*dirIn
				! Check for correct direction of refraction
				temp1 = acos(dot_product(dirOut,-fcNorm))
				temp2 = sin(temp1)
				if(abs(temp2-sinThT) .gt. NANO) then
					write(*,'(a)') "trouble refracting"
				end if
				return
			end if
		end if
	end subroutine transBetweenDoms

	subroutine transExit(ec,fcNum,dirIn,trans,dirOut)
		integer,intent(in) :: fcNum
		integer :: remNo,fcNodes(3)
		real(8),parameter :: n2 = 1.d0	! Air is the second medium
		real(8) :: rIndRatio,th1,th2,cosInc,sinTht,rhoFresnel,frRefChk,&
		temp1,temp2,remVert(3),fcNorm(3),fcVerts(3,3)
		real(8),intent(in) :: ec(4,3),dirIn(3)
		real(8),intent(out):: dirOut(3)
		logical,intent(out) :: trans

		fcNodes = getFaceNodes(fcNum)
		remNo = 10-sum(fcNodes)
		fcVerts = ec(fcNodes,:)
		remVert = ec(remNo,:)
		fcNorm = getFaceNorm(fcVerts,remVert,.true.)
		cosInc = -dot_product(dirIn,fcNorm)
		if (cosInc .lt. 0.0d0) then
            write(*,*) "change direction of normal"
            cosInc = -cosInc
            fcNorm = -fcNorm
        end if
		th1 = acos(cosInc)
		rIndRatio = rtRefrInd(rtCurrDom,rtCurrLambda)/n2
		sinTht = rIndRatio*sqrt(1-cosInc**2.d0)
		th2 = asin(sinTht)
		if(sinTht .gt. 1.d0) then
			trans = .false.
			dirOut = specularReflection(ec,fcNum,dirIn)
			return
		else
			call getFresnelReflectivity(th1,th2,rhoFresnel)
			call random_number(frRefChk)
			if(frRefChk .lt. rhoFresnel) then
				trans = .false.
				dirOut = specularReflection(ec,fcNum,dirIn)
				return
			else
				trans = .true.
				dirOut = (rIndRatio*cosInc-sqrt(1.d0-sinTht**2.d0))*fcNorm
				dirOut = dirOut + rIndRatio*dirIn
				! Check for correct direction of refraction
				temp1 = acos(dot_product(dirOut,-fcNorm))
				temp2 = sin(temp1)
				if(abs(temp2-sinThT) .gt. NANO) then
					write(*,'(a)') "trouble refracting"
				end if
				return
			end if
		end if
	end subroutine transExit

	subroutine getFresnelReflectivity(th1,th2,rhoFresnel)
		real(8) :: thDiff,thSum
		real(8),intent(in) :: th1,th2
		real(8),intent(out) :: rhoFresnel

		thDiff = th1-th2
		thSum = th1+th2
		rhoFresnel = 0.5d0*((tan(thDiff)**2.d0/tan(thSum)**2.d0) + 		&
		(sin(thDiff)**2.d0/sin(thSum)**2.d0))
	end subroutine getFresnelReflectivity

	subroutine setMediumValues(valIn,valNameFlag)
		real(8),intent(in) :: valIn
		character(*),intent(in) :: valNameFlag

		if(valNameFlag == "k") then
			rtKappa = valIn
		elseif(valNameFlag == "s") then
			rtSigma = valIn
		else
			write(*,'(a)') "Property to be set not recognised."
			stop
		end if
	end subroutine setMediumValues

	subroutine getElementNumRays(elNum,elNumRays,elRayPow)
		integer :: elDom,elNodes(4)
		integer,intent(in) :: elNum
		integer,intent(out) :: elNumRays
		real(8) :: elVol,elCentTemp,elEmPow,elKappa,elNoTemps(4)
		real(8),intent(out) :: elRayPow

		elDom = meshElems(elNum)%domain
		elNodes = meshElems(elNum)%nodes
		elNoTemps = meshTemperatures(elNodes)
		elVol = getElementVolume(elNum)
		elCentTemp = sum(elNoTemps)/4.0d0
		elKappa = rtKappa(elDom,rtCurrLambda)
		elEmPow = 4.0d0*elKappa*sigB*(elCentTemp**4.0d0)*elVol
		elNumRays = nint(elEmPow/rtRefRayPow)
		if(elNumRays .lt. rtElemMinRays) then
			elNumRays = rtElemMinRays
		end if
		elRayPow = elEmPow/real(elNumRays,8)
	end subroutine getElementNumRays

	subroutine getLargeDiffSphMu(rMu,mu)
		integer :: interLoc
		real(8),parameter :: muMin = 0.0d0, muMax = 1.0d0
		real(8) :: randLv,randhV,lV,hV,ratio,intLocReal,numReal
		real(8),intent(in) :: rMu
		real(8),intent(out) :: mu

		interLoc = floor(rMu*real(rtNumPFTable,8))
		if(interLoc == 0) then
			interLoc = 1
		end if
		lV = rtPFTable(interLoc)
		hV = rtPFTable(interLoc+1)
		intLocReal = real(interLoc,8)
		numReal = real(rtNumPFTable,8)
		randLv = (intLocReal-1)*1.0d0/numReal
		randhV = (intLocReal)*1.0d0/numReal
		mu = lV + ((rMu-randLv)/(randHv-randLv))*(hV-lV)
		if(mu .lt. -1.d0) mu = -1.d0
		if(mu .gt. 1.d0) mu = 1.d0
	end subroutine getLargeDiffSphMu

	subroutine getHeyneyGreensteinMu(rMu,mu)
		real(8),intent(in) :: rMu
		real(8),intent(out) :: mu
		real(8) :: s,rMuNew,g

		g = rtAnisFac(rtCurrDom,rtCurrLambda)
		if(abs(g) .lt. MICRO) then
			mu = 1.d0 - 2.d0*rMu
			return
		end if

		s = 2.d0*rMu - 1.d0
		mu = (1+g**2.d0-((1-g**2.d0)/(1+g*s))**2.d0)
		mu = (1/(2*g))*mu
		do while(abs(mu) .gt. 1.d0)
			call random_number(rMuNew)
			s = 2.d0*rMuNew - 1.d0
			mu = (1+g**2.d0-((1-g**2.d0)/(1+g*s))**2.d0)
			mu = (1/(2*g))*mu
		end do
	end subroutine getHeyneyGreensteinMu

	subroutine popLargeSphDiffMuTable()
		integer :: i,nSplSteps,lV
		real(8) :: j,h,m,a,b,c,d,r,stepSize
		real(8),allocatable :: x(:),y(:),yy(:)

		nSplSteps = 100
		h = 2.0d0/nSplSteps
		allocate(x(nSplSteps+1))
		allocate(y(nSplSteps+1))
		do i=0,nSplSteps
			x(i+1) = -1.0d0 + h*(i)
			y(i+1) = (2.d0/(3.d0*pi))*((x(i+1)**2.d0)*acos(x(i+1)) 		&
			-asin(x(i+1))/2.d0-(3.d0/2.d0)*x(i+1)*sqrt(1-x(i+1)**2.d0)	&
			+ pi/4.d0)
		end do

		call cubicSplineNaturalFit(y,x,yy)

		if(.not.(allocated(rtPFTable))) then
			allocate(rtPFTable(rtNumPFTable))
		end if
		stepSize = 1.0d0/real(rtNumPFTable,8)
		i = 0
		do m=0.0d0,0.999999999d0,stepSize
			i = i+1
			do lV=1,nSplSteps-1
				r = (m-y(lV))/(y(lV+1)-y(lV))
				if((r .gt. 0.0d0) .and. (r .lt. 1.0d0)) then
					exit
				end if
			end do
			a = (y(lV+1) - m)/(y(lV+1)-y(lV))
			b = 1.d0-a
			c = (1.d0/6.d0)*(a**3.d0-a)*(y(lV+1)-y(lV))**2.d0
			d = (1.d0/6.d0)*(b**3.d0-b)*(y(lV+1)-y(lV))**2.d0
			rtPFTable(i) = a*x(lV) + b*x(lV+1) + c*yy(lV) + d*yy(lV+1)
		end do
		rtPFTable(rtNumPFTable) = -1.0d0
		open(1864,file="../results/rtPfTable.out")
		write(1864,'(f12.9)') rtPfTable
		close(1864)
	end subroutine popLargeSphDiffMuTable

	subroutine getUnifCloudCoeffs(alpha,rho,dP,Nt)
		integer,intent(in) :: Nt
		real(8),intent(in) :: alpha,rho,dp

		rtKappa = pi*((dp/2.0d0)**2.0d0)*Nt*alpha
		rtSigma = pi*((dp/2.0d0)**2.0d0)*Nt*rho
		rtBeta = rtKappa + rtSigma
	end subroutine getUnifCloudCoeffs

	function getRaySphDir() result(dir)
		real(8) :: r1,r2,th,ph,dir(3)

		call random_number(r1)
		call random_number(r2)
		ph = 2*pi*r1
		th = acos(1.0d0 - 2.0d0*r2)
		dir = getDirectionCoords(th,ph)
	end function getRaySphDir

	function getFaceDiffRayDir(fcVerts,fcNorm) result(dir)
		real(8),intent(in) :: fcVerts(3),fcNorm(3)
		real(8):: th,ph,dir(3)

		call random_number(th)
		th = asin(sqrt(th))
		call random_number(ph)
    	ph = 2.0d0*pi*ph
		dir = getDirectionCoords(th,ph)

		if(dot_product(fcNorm,dir) .lt. 0.d0) then
			do while(dot_product(fcNorm,dir).lt. 0.d0)
				call random_number(th)
				th = asin(sqrt(th))
				call random_number(ph)
				ph = 2.0d0*pi*ph
				dir = getDirectionCoords(th,ph)				
			end do

		end if

	end function getFaceDiffRayDir

	function getDirectionCoords(th,ph) result(dir)
		real(8) :: th,ph,dir(3)

		dir(1) = sin(th)*cos(ph)
		dir(2) = sin(th)*sin(ph)
		dir(3) = cos(th)
		dir = dir/norm2(dir)
	end function getDirectionCoords

	function getRayPathIntegral() result(pInt)
		real(8) :: rBeta,pInt

		call random_number(rBeta)
		pInt = log(1.d0/rBeta)
	end function getRayPathIntegral


	subroutine createNonUniformEmissionSurface()
!	Placeholder for future use in case of non-uniform surface 
!	temperatures
	end subroutine createNonUniformEmissionSurface

!------------------------------------------------------------------------
!	LED algorithms for multidomain tracing
!------------------------------------------------------------------------

	subroutine traceFromSurfLED()
		integer,parameter :: nSfEmFil=191,nSfAbFil=192,nScrPts=193,		&
		nSctDirs=194,nReEmDrp=195,nExitPts=196,nOutPts=197,nBotPts=198,	&
		lcYellow=6,lcBlue=3
		integer :: i,j,stEl,emEl,emFc,endEl,outPtCt,reEmCt,reEmDropCt,	&
		lambda,elNodes(4),plcHlder
		real(8) :: pInt,rAbs,rReEm,tcos,pt(3),dir(3),endPt(3),dirOut(3),&
		endDir(3),ptScr(3),spFnVals(4)
		logical :: outPt,vExit,byAbs,scatter,reEmission
		character(*),parameter :: wFmt = '(a,2x,i8)'
		character :: fSfEmPts*72,fSfAbPts*72,fSctDirs*72,fReEmDrp*72,	&
					 fExitPts*72,fScrPts*72,fOutPts*72,fBotPts*72


		fSfEmPts=commResDir//trim(adjustl(rtfResPre))//"_surfems.out"
		fSfAbPts=commResDir//trim(adjustl(rtfResPre))//"_surfpts.out"
		fSctDirs=commResDir//trim(adjustl(rtfResPre))//"_sctDirs.out"
		fReEmDrp=commResDir//trim(adjustl(rtfResPre))//"_reEmDrp.out"
		fExitPts=commResDir//trim(adjustl(rtfResPre))//"_exitPts.out"
		fScrPts=commResDir//trim(adjustl(rtfResPre))//"_scrPts.out"
		fOutPts=commResDir//trim(adjustl(rtfResPre))//"_outPts.out"
		fBotPts=commResDir//trim(adjustl(rtfResPre))//"_botPts.out"
		open(nSfEmFil,file=fSfEmPts)
		open(nSfAbFil,file=fSfAbPts)
		open(nScrPts,file=fScrPts)
		open(nSctDirs,file=fSctDirs)
		open(nReEmDrp,file=fReEmDrp)
		open(nExitPts,file=fExitPts)
		open(nBotPts,file=fBotPts)
		open(nOutPts,file=fOutPts,position='append')
		outPtCt = 0
		reEmDropCt = 0

		do i=1,rtMCNumRays
			if(mod(i,rtMCNumRays/10).eq.0) write(*,*) "Ray: ", i
			call startRayFromSurf(emEl,emFc,pInt,pt,dir)
			write(nSfEmFil,'(6(f15.12,2x))') pt,dir
			lambda = lcBlue
			reEmCt = 0
			rtCurrLambda = 1		! Hard coded, needs to change
			rtCurrDom = meshElems(emEl)%domain

!			Tracing line
			stEl = emEl
100			rtNtraces = rtNtraces + 1
			if(rtNtraces .eq. 61) then
				plcHlder = rtNtraces
			end if
			call traceTrans(i,stEl,pt,dir,pInt,nOutPts,outPt,vExit,		&
			byAbs,endEl,endPt,endDir)
!			Handle the exit of a ray from the volume
			if(vExit) then
				call handleExit(outPt,outPtCt,endPt,endDir,lambda,		&
				nExitPts,nScrPts,nBotPts)
				if(outPtCt .gt. rtLimOutPts) then
					write(*,*)"Leaving tracing routine now."
					return
				end if
				cycle
			end if

			if(byAbs) then
				rtElemAbs(endEl) = rtElemAbs(endEl) + 1
				write(nSfAbFil,'(3(f15.12,2x))') endPt
				call shapeFunctionsAtPoint(endEl,endPt,spFnVals)
				elNodes = meshElems(endEl)%nodes
				rtNodalSrc(elNodes) = rtNodalSrc(elNodes) + &
				rtRefRayPow*spFnVals
				rtAbsNoReEmCt = rtAbsNoReEmCt + 1
				cycle
			end if

			call checkScatter(scatter)
			if(scatter) then
				reEmCt = reEmCt+1
				if(reEmCt .gt. MEGA) then
					write(*,*) "Too many scatters for one ray: ", i
					write(nReEmDrp,'(3(e16.9,2x),i8)') pt,i
					reEmDropCt = reEmDropCt + 1
					cycle
					if(reEmDropCt .gt. rtLimReEmDrops) then
						write(*,*)"Re-emission error count reached limit."
						stop
					end if
				end if
				pInt = getRayPathIntegral()
				pt = endPt
				stEl = endEl
				if(rtSctThr .eq. rtSctLrgSphThr) then
					dir = scatterRayLargeSphereCloud()
				elseif(rtSctThr .eq. rtSctMieThr) then
					dir = scatterRayHeyneyGreenStein()
				else
					dir = scatterRayIsotropic()
				end if
				tcos = dot_product(endDir,dir)
				write(nSctDirs,'(4(e16.9,2x))') dir,tcos
				goto 100
			else
				call checkReEmission(reEmission)
				if(reEmission) then
					reEmCt = reEmCt+1
					if(reEmCt .gt. MEGA) then
						write(*,*) "Point dropped in ReEm loop: ", i
						write(nReEmDrp,'(3(e16.9,2x),i8)') pt,i
						reEmDropCt = reEmDropCt + 1
						cycle
						if(reEmDropCt .gt. rtLimReEmDrops) then
							write(*,*)"Too many points dropped in reEmission."
							stop
						end if
					end if
					lambda = lcYellow
					rtCurrLambda = 2		! Hard coded, needs to change
					pInt = getRayPathIntegral()
					pt = endPt
					stEl = endEl
					dir = getRaySphDir()
					goto 100
				else
					rtElemAbs(endEl) = rtElemAbs(endEl) + 1
					write(nSfAbFil,'(3(f15.12,2x))') endPt
					call shapeFunctionsAtPoint(endEl,endPt,spFnVals)
					elNodes = meshElems(endEl)%nodes
					rtNodalSrc(elNodes) = rtNodalSrc(elNodes) + &
					rtRefRayPow*spFnVals
					rtInVolAbsCt = rtInVolAbsCt + 1
					rtAbsNoReEmCt = rtAbsNoReEmCt + 1
				end if
			end if
		end do
		if(reEmDropCt .eq. 0) then
			write(nReEmDrp,'(a)') "No points dropped in reEmission loop."
		end if
		write(*,wFmt)"Out of face points: ", outPtCt
		write(*,wFmt)"Total number of absorptions: ", rtAbsNoReEmCt
		write(*,wFmt)"Absorptions inside volume: ", rtInVolAbsCt
		write(*,wFmt)"Absorptions at black surfaces: ",rtBBAbsNum
		write(*,wFmt)"Absorptions at non-black surfaces: ",rtParSfAbsNum
		write(*,wFmt)"Exits at transparent surfaces: ",rtTrExitNum
		write(*,wFmt)"Number of transmissions across neighbouring domains: ", rtTransBetnDomsCt
		write(*,wFmt)"Number of actual traces: ", rtNtraces
		close(nSfEmFil)
		close(nSfAbFil)
		close(nScrPts)
		close(nSctDirs)
		close(nReEmDrp)
		close(nExitPts)
		close(nOutPts)
		close(nBotPts)
	end subroutine traceFromSurfLED

	subroutine checkScatter(sctr)
		real(8) :: rAbs,thrVal
		logical,intent(out) :: sctr

		call random_number(rAbs)
		thrVal = rtAbsThr(rtCurrDom,rtCurrLambda)
		if(rAbs .gt. thrVal) then
			sctr = .true.
		else
			sctr = .false.
		end if
	end subroutine checkScatter

	subroutine checkReEmission(reEm)
		real(8) :: rReEm,thrVal
		logical,intent(out) :: reEm

		call random_number(rReEm)
		thrVal = rtReEmThr(rtCurrDom,rtCurrLambda)
		if(rReEm .lt. thrVal) then
			reEm = .true.
		else
			reEm = .false.
		end if
	end subroutine checkReEmission

!------------------------------------------------------------------------
!	END algorithms for multidomain tracing
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!	The following routine is a general multi-domain tracing routine
!	with refraction across domain boundaries and strategies for exit
!	from system from the previous code improved upon.
!------------------------------------------------------------------------

	subroutine traceTrans(rN,stEl,stPt,stDir,pInt,nFOutPts,outPt,vExit,	&
	byAbs,endEl,endPt,endDir)
		integer :: i,noNum,rayIterCt,nhbrFc,cEl,newFc,cDom,nhbrEl,		&
		nhbrDom,nhbrSf,parById,elNodes(4)
		integer,intent(in) :: stEl,rN,nFOutPts
		integer,intent(out) :: endEl
		real(8) :: lTrav,lToFc,optPath,cBeta,pt(3),dir(3),newDir(3),	&
		ec(4,3)
		real(8),intent(in) :: pInt,stPt(3),stDir(3)
		real(8),intent(out) :: endPt(3),endDir(3)
		logical :: inFc,bbSfChk,trSfChk,npSfChk,trans,parSfChk,			&
		parDiffSfChk,parSpecSfChk,byAbsChk,actByChk
		logical,intent(out) :: outPt,vExit,byAbs

		cEl = stEl
		pt = stPt
		dir = stDir
		rayIterCt = 0
		endEl = 0
		endPt = 0.0d0
		optPath = 0.0d0
		endDir = 0.0d0
		outPt = .false.
		vExit = .false.
		byAbs = .false.
		do while(optPath.lt.pInt)
			rtCurrDom = meshElems(cEl)%domain
			if(rayIterCt .gt. MEGA) then
				outPt = .true.
				vExit = .true.
				write(*,*)"Ray entered interminable loop"
				exit					
			end if
			bbSfChk = .false.
			trSfChk = .false.
			npSfChk = .false.
!			Initialise the face-check values
			elNodes = meshElems(cEl)%nodes
!			cDom = meshElems(cEl)%domain
!			cBeta = rtBeta(cDom)
			cBeta = rtBeta(rtCurrDom,rtCurrLambda)
			ec = meshVerts(elNodes,:)
			pt = pt + DPICO*dir
			optPath = optPath + DPICO*cBeta
			call getNextFace(ec,pt,dir,newFc,lToFc)
			optPath = optPath + lToFc*cBeta
			if(optPath .ge. pInt) then
				pt = pt + (optPath/cBeta-pInt/cBeta)*dir
				endEl = cEl
				endPt = pt
				endDir = dir
				exit
			end if
			pt = pt + lToFc*dir
			inFc = checkNewPt(pt,dir,ec,newFc)
			if(.not. inFc) then
				write(*,*) "Point traced not within face."
				write(*,*) "Elnum: ",cEl, "Raynum: ",rN
				write(nFOutPts,'(a)') "Element, face and length: "
				write(nFOutPts,'(i8,2x,i2,2x,f16.9)') cEl,newFc,lToFc
				write(nFOutPts,'(a)') "Point and direction: "
				write(nFOutPts,'(3(f16.9))') pt
				write(nFOutPts,'(3(f16.9))') dir
				write(nFOutPts,'(a)') "Element coords: "
				do noNum=1,4
					write(nFOutPts,'(3(f16.9))') ec(noNum,:)
				end do
				outPt = .true.
				vExit = .true.
				optPath = MEGA
				exit
			end if
!			Check if the neighbouring face is a system boundary or iface
			nhbrFc = meshElems(cEl)%neighbours(newFc,2)
			nhbrSf = meshElems(cEl)%neighbours(newFc,3)
			actByChk = (count(nhbrSf==rtActBys) .gt. 0)
			if(actByChk) then
				bbSfChk = (count(nhbrSf==rtBBSfIds) .gt. 0)
				trSfChk = (count(nhbrSf==rtTrSfIds) .gt. 0)
				npSfChk = (count(nhbrSf==rtNPSfIds) .gt. 0)
				parSfChk = (count(nhbrSf==rtParSfIds) .gt. 0)
				if(bbSfChk) then
					endEl = cEl
					endPt = pt
					endDir = dir
					byAbs = .true.
					optPath = MEGA
					rtBBAbsNum = rtBBAbsNum + 1
					exit
				end if
				if(trSfChk) then
					call transExit(ec,newFc,dir,trans,newDir)
					if(trans) then
						endEl = cEl
						endPt = pt
						endDir = newDir
						optPath = MEGA
						vExit = .true.
						rtTrExitNum = rtTrExitNum + 1
						exit
					else
						dir = newDir
					end if
				end if
				if(npSfChk) then
					newDir = specularReflection(ec,newFc,dir)
					dir = newDir
				end if
				if(parSfChk) then
					parDiffSfChk = (count(nhbrSf==rtParDiffSfIds).gt. 0)
					parSpecSfChk = (count(nhbrSf==rtParSpecSfIds).gt. 0)
					if(parDiffSfChk) then
						parById = minloc(rtParDiffSfIds,1,rtParDiffSfIds==nhbrSf)
						call diffParBoundary(ec,newFc,parById,newDir,	&
						byAbsChk)
						if(byAbsChk) then
							endEl = cEl
							endPt = pt
							endDir = dir
							optPath = MEGA
							rtParSfAbsNum = rtParSfAbsNum + 1
							byAbs = .true.
							exit
						else
							dir = newDir
						end if
					end if
					if(parSpecSfChk) then
						parById = minloc(rtParSpecSfIds,1,rtParSpecSfIds==nhbrSf)
						call specParBoundary(ec,newFc,parById,dir,		&
						newDir,byAbsChk)
						if(byAbsChk) then
							endEl = cEl
							endPt = pt
							endDir = dir
							optPath = MEGA
							rtParSfAbsNum = rtParSfAbsNum + 1
							byAbs = .true.
							exit
						else
							dir = newDir
						end if
					end if
				end if
			else
!			Check for neighbouring faces
				nhbrEl = meshElems(cEl)%neighbours(newFc,1)
				nhbrDom = meshElems(nhbrEl)%domain
				if(nhbrDom .ne. rtCurrDom) then
					call transBetweenDoms(ec,newFc,nhbrDom,dir,trans,	&
					newDir)
					dir = newDir
					if(trans) then
						rtTransBetnDomsCt = rtTransBetnDomsCt+1
						cEl = nhbrEl
					else
! 						Statement added for Total Internal Reflection
						cEl = cEl
					end if
				else
					cEl = nhbrEl
				end if
				rayIterCt = rayIterCt + 1
			end if
		end do
	end subroutine traceTrans

	subroutine diffParBoundary(ec,fcNum,parById,dirOut,byAbsChk)
		integer,intent(in) :: fcNum,parById
		real(8) :: rByAbs
		real(8),intent(in) :: ec(4,3)
		real(8),intent(out) :: dirOut(3)
		logical,intent(out) :: byAbsChk

		call random_number(rByAbs)
		if(rByAbs .lt. rtParDiffSfVals(parById)) then
			dirOut = diffuseReflection(ec,fcNum)
			byAbsChk = .false.
		else
			dirOut = 0.d0
			byAbsChk = .true.
		end if
	end subroutine diffParBoundary

	subroutine specParBoundary(ec,fcNum,parById,dirIn,dirOut,byAbsChk)
		integer,intent(in) :: fcNum,parById
		real(8) :: rByAbs
		real(8),intent(in) :: dirIn(3),ec(4,3)
		real(8),intent(out) :: dirOut(3)
		logical,intent(out) :: byAbsChk

		call random_number(rByAbs)
		if(rByAbs .lt. rtParSpecSfVals(parById)) then
			dirOut = specularReflection(ec,fcNum,dirIn)
			byAbsChk = .false.
		else
			dirOut = 0.d0
			byAbsChk = .true.
		end if
	end subroutine specParBoundary

!------------------------------------------------------------------------
!	END multidomain tracing routine
!------------------------------------------------------------------------


!------------------------------------------------------------------------
!	Ray tracing routines for 1-D gray medium without transmission lie
!	hereon. These are not used for the LED case
!------------------------------------------------------------------------

	subroutine traceFromVol()
		integer,parameter :: limoutpt=5
		integer,parameter :: nVlAbPtsFil=286,nVlEmPtsFil=208
		integer :: i,j,k,cEl,elNr,endEl,outPtCt,elNodes(4)
		integer,allocatable :: emSfIds(:)
		real(8) :: pL,pt(3),dir(3),endPt(3),spFnVals(4)
		logical :: outPt
		character(*),parameter :: vlEmPtsFil=commResDir//"volems.out",	&
								  vlAbPtsFil=commResDir//"volpts.out"

		open(nVlEmPtsFil,file=vlEmPtsFil)
		open(nVlAbPtsFil,file=vlAbPtsFil)
		outPtCt = 0
		do i=1,meshNumElems
			elNr = rtElemSto(i)
			if(elNr .ne. 0) then
				do j=1,elNr
					call startRayInVolume(i,pL,pt,dir)
					write(nVlEmPtsFil,'(6(f15.12,2x))') pt,dir
					cEl = i
					call traceSingleRay(pt,dir,pL,cEl,outPt,endEl,endPt)
					if(outPt) then
						outPtCt = outPtCt + 1
						if(outptct .ge. limoutpt) then
							write(*,*)"Count of dropped points reached&
							& limit."
							stop
						end if
					end if
					if(endEl .ne. 0) then
						rtElemAbs(endEl) = rtElemAbs(endEl) + 1
						write(nVlAbPtsFil,'(3(f15.12,2x))') endPt
					end if
				end do
			end if
		end do
		close(nVlAbPtsFil)
		close(nVlEmPtsFil)
	end subroutine traceFromVol

	subroutine traceVolPowerBased()
		integer,parameter :: limoutpt=5
		integer,parameter :: nVlAbPtsFil=286,nVlEmPtsFil=208
		integer :: i,j,k,cEl,elNumRays,endEl,outPtCt,elNodes(4)
		integer,allocatable :: emSfIds(:)
		real(8) :: pL,elRayPow,pt(3),dir(3),endPt(3),spFnVals(4)
		logical :: outPt
		character(*),parameter :: vlEmPtsFil=commResDir//"volems.out",	&
								  vlAbPtsFil=commResDir//"volpts.out"

		open(nVlEmPtsFil,file=vlEmPtsFil)
		open(nVlAbPtsFil,file=vlAbPtsFil)
		rtNodalSrc = 0.0d0
		outPtCt = 0
		do i=1,meshNumElems
			call getElementNumRays(i,elNumRays,elRayPow)
				do j=1,elNumRays
					call startRayInVolume(i,pL,pt,dir)
					write(nVlEmPtsFil,'(6(f15.12,2x))') pt,dir
					call shapeFunctionsAtPoint(i,pt,spFnVals)
					elNodes = meshElems(i)%nodes
					rtNodalSrc(elNodes) = rtNodalSrc(elNodes) - 		&
					elRayPow*spFnVals
					cEl = i
					call traceSingleRay(pt,dir,pL,cEl,outPt,endEl,endPt)
					if(outPt) then
						outPtCt = outPtCt + 1
						if(outptct .ge. limoutpt) then
							write(*,*)"Count of dropped points reached&
							& limit."
							stop
						end if
					end if
					if(endEl .ne. 0) then
						rtElemAbs(endEl) = rtElemAbs(endEl) + 1
						write(nVlAbPtsFil,'(3(f15.12,2x))') endPt
						call shapeFunctionsAtPoint(endEl,endPt,spFnVals)
						elNodes = meshElems(endEl)%nodes
						rtNodalSrc(elNodes) = rtNodalSrc(elNodes) + 	&
						elRayPow*spFnVals
					end if
				end do
		end do
		close(nVlAbPtsFil)
		close(nVlEmPtsFil)
		rtNodalSrc = rtNodalSrc + rtWallSrc
		open(975,file="../obj/tempRadSrc.out")
			write(975,'(f20.13)') rtNodalSrc
		close(975)
	end subroutine traceVolPowerBased

!	subroutine traceFromSurf(pRatio)
!		integer,parameter :: limoutpt=5
!		integer,parameter :: nSfAbPtsFil=175,nSfEmPtsFil=197
!		integer :: i,j,cEl,emEl,emFc,endEl,outPtCt,elNodes(4)
!		real(8) :: pL,pt(3),dir(3),endPt(3),spFnVals(4)
!		real(8),intent(in) :: pRatio(:)
!		logical :: outPt
!		character(*),parameter :: sfEmPtsFil=commResDir//"surfems.out",	&
!								  sfAbPtsFil=commResDir//"surfpts.out"

!		open(nSfAbPtsFil,file=sfAbPtsFil)
!		open(nSfEmPtsFil,file=sfEmPtsFil)
!		rtNodalSrc = 0.0d0
!		outPtCt = 0
!		do i=1,rtMCNumRays
!			call startRayFromSurf(pRatio,emEl,emFc,pL,pt,dir)
!			write(nSfEmPtsFil,'(6(f15.12,2x))') pt,dir
!			cEl = emEl
!			call traceSingleRay(pt,dir,pL,cEl,outPt,endEl,endPt)
!			if(outPt) then
!				outPtCt = outPtCt + 1
!				if(outptct .ge. limoutpt) then
!					write(*,*)"Count of dropped points reached limit."
!					stop
!				end if
!			end if
!			if(endEl .ne. 0) then
!				rtElemAbs(endEl) = rtElemAbs(endEl) + 1
!				write(nSfAbPtsFil,'(3(f15.12,2x))') endPt
!				call shapeFunctionsAtPoint(endEl,endPt,spFnVals)
!				elNodes = meshElems(endEl)%nodes
!				rtNodalSrc(elNodes) = rtNodalSrc(elNodes) + 			&
!				rtRefRayPow*spFnVals
!			end if
!		end do
!		close(nSfEmPtsFil)
!		close(nSfAbPtsFil)
!		open(975,file="../obj/tempRadSrc.out")
!			write(975,'(f20.13)') rtNodalSrc
!		close(975)
!		rtWallSrc = rtNodalSrc
!	end subroutine traceFromSurf

	subroutine traceSingleRay(pt,dir,pL,cEl,outPt,endEl,endPt)
		integer :: i,rayIterCt,chCt,newFc,nEmSfs,nhbrFc,elNodes(4)
		integer,intent(inout) :: cEl
		integer,intent(out) :: endEl
		integer,allocatable :: emSfIds(:)
		real(8) :: lTrav,lToFc,newDir(3),ec(4,3)
		real(8),intent(in) :: pL
		real(8),intent(out) :: endPt(3)
		real(8),intent(inout) :: pt(3),dir(3)
		logical :: inFc
		logical,intent(out) :: outPt

		rayIterCt = 0
		endEl = 0
		endPt = 0.0d0
		lTrav = 0.0d0
		outPt = .false.
		do while(lTrav.lt.pL)
			elNodes = meshElems(cEl)%nodes
			ec = meshVerts(elNodes,:)
			pt = pt + DPICO*dir
			lTrav = lTrav+DPICO
			chCt = 0
			if(rayIterCt .gt. MEGA) then
				outPt = .true.
				write(*,*)"Ray entered interminable loop"
				exit					
			end if
			call getNextFace(ec,pt,dir,newFc,lToFc)
			lTrav = lTrav + lToFc
			if(lTrav.ge.pL) then
				pt = pt + (pL-(lTrav-lToFc))*dir
				endEl = cEl
				endPt = pt
				exit
			end if
			pt = pt + lToFc*dir
			inFc = checkNewPt(pt,dir,ec,newFc)
			if(.not. inFc) then
				write(*,*) "Point traced not within face."
				write(*,*) "Elnum: ",cEl, "FaceNum: ",newFc
				outPt = .true.
				exit
			end if
			nhbrFc = meshElems(cEl)%neighbours(newFc,2)
			if(nhbrFc .lt. 0) then
				chCt = count(-nhbrFc==rtEmSfIds)
				if(chCt == 0) then
					newDir = specularReflection(ec,newFc,dir)
					dir = newDir
				else
					lTrav = MEGA
					exit
				end if
			end if
			cEl = meshElems(cEl)%neighbours(newFc,1)
			rayIterCt = rayIterCt + 1
		end do
	end subroutine traceSingleRay

	subroutine createEmissionSurfaces()
		integer :: i,j,nEmSf,currSurf,nEmFcs,elNum,fcNum,fcNodes(3),	&
		elNodes(4)
		real(8) :: fcEmPow,fcArea,fcCentT,emSfPow,fcNoTs(3)

		nEmSf = size(rtEmSfIds,1)
		if(.not.(allocated(rtEmSurfs))) then
			allocate(rtEmSurfs(nEmSf))
		end if
		do i=1,nEmSf
			currSurf = rtEmSfIds(i)
			rtEmSurfs(i)%emSurf = meshSurfs(currSurf)
			nEmFcs = meshSurfs(currSurf)%numFcs
			if(.not.(allocated(rtEmSurfs(i)%cuSumFcEmPow))) then
				allocate(rtEmSurfs(i)%cuSumFcEmPow(nEmFcs))
			end if
			emSfPow = 0.d0
			do j=1,nEmFcs
				elNum = rtEmSurfs(i)%emSurf%elNum(j)
				elNodes = meshElems(elNum)%nodes
				fcNum = rtEmSurfs(i)%emSurf%fcNum(j)
				fcArea = rtEmSurfs(i)%emSurf%fcArea(j)
				fcNodes = getFaceNodes(fcNum)
				fcNoTs = meshTemperatures(elNodes(fcNodes))
				fcCentT = sum(fcNoTs)/3.0d0
				fcEmPow = sigB*fcArea*(fcCentT**4.0d0)
				emSfPow = emSfPow+fcEmPow
				rtEmSurfs(i)%cuSumFcEmPow(j) = emSfPow
			end do
			rtEmSurfs(i)%totEmPow = emSfPow
			rtEmSurfs(i)%cuSumFcEmPow=rtEmSurfs(i)%cuSumFcEmPow/emSfPow
		end do
	end subroutine createEmissionSurfaces
!------------------------------------------------------------------------
!	END of 1-D gray medium subroutines
!------------------------------------------------------------------------

end module rt
