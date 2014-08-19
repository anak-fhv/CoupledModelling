module rtHelper

	use commonRoutines
	use mesh

	implicit none

!	Declare all relevant types
	type emissionSurface
		real(8) :: totEmPow
		real(8),allocatable :: cuSumFcEmPow(:)
		type(surface) :: emSurf
	end type emissionSurface

!	Module variables
	integer :: rtNumRays
	integer,allocatable :: 
	real(8) :: kappa,sigma
	real(8),allocatable ::
	type(emissionSurface),allocatable :: rtEmSurfs(:)
!	Inputs
	character(*),parameter ::

	contains

	subroutine traceFromSurf(pRatio,elemAbsNums)
		integer,parameter :: limoutpt=5
		integer :: i,j,cEl,emEl,emFc,endEl,outPtCt
		integer,intent(inout) :: elemAbsNums(:)
		real(8) :: pL,pt(3),dir(3),endPt(3)
		real(8),intent(in) :: pRatio(:)
		logical :: outPt

		outPtCt = 0
		do i=1,nRays
			call startRayFromSurf(pRatio,emEl,emFc,pL,pt,dir)
			cEl = emEl
			call traceSingleRay(pt,dir,pL,cEl,outPt,endEl,endPt)
			if(outPt) then
				outPtCt = outPtCt + 1
				if(outptct .ge. limoutpt) then
					write(*,*)"Count of dropped points reached limit."
					stop
				end if
			end if
			if(endEl .ne. 0) then
				elemAbsNums(endEl) = elemAbsNums(endEl) + 1
			end if

		end do
	end subroutine traceFromSurf

	subroutine traceFromVol(elemStoNums,elemAbsNums)
		integer,parameter :: limoutpt=5
		integer :: i,j,k,cEl,elNr,endEl,outPtCt
		integer,intent(in) :: elemStoNums(:)
		integer,intent(inout) :: elemAbsNums(:)
		real(8) :: pL,pt(3),dir(3),endPt(3)
		logical :: outPt

		outPtCt = 0
		do i=1,meshNumElems
			elNr = elemStoNums(i)
			if(elNr .ne. 0) then
				do j=1,elNr
					call startRayInVolume(i,pL,pt,dir)
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
						elemAbsNums(endEl) = elemAbsNums(endEl) + 1
					end if
				end do
			end if
		end do
	end subroutine traceFromVol

	subroutine traceSingleRay(pt,dir,pL,cEl,outPt,endEl,endPt)
		integer :: i,rayIterCt,chCt,cEl,newFc,nEmSfs,nhbrFc,elNodes(4)
		integer,intent(in) :: cEl
		integer,intent(out) :: endEl
		real(8) :: lTrav,lToFc,newDir(3),ec(4,3)
		real,intent(in) :: pL,pt(3),dir(3)
		real,intent(out) :: endPt(3)
		logical :: inFc
		logical,intent(out) :: outPt

		nEmSfs = size(rtEmSurfs,1)
		rayIterCt = 0
		endEl = 0
		endPt = 0.0d0
		lTrav = 0.0d0
		outPt = .false.
		do while(lTrav.lt.pL)
			elNodes = meshElems(cEl)%nodes
			ec = vertices(elNodes,:)
			pt = pt + PICO*dir
			pL = pL-PICO
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
				write(*,*) "Elnum: ",cEl, "Raynum: ",i
				outPt = .true.
				exit
			end if
			nhbrFc = meshElems(cEl)%neighbours(newFc,2)
			if(nhbrFc .lt. 0) then
				chCt = count(-nhbrFc.eq.rtEmSurfs%emSurf%sfId)
				if(chCt .eq. 0) then
					call specRef(ec,newFc,dir,newDir)
					dir = newDir
				else
					pL = 0.0d0
					exit
				end if
			end if
			cEl = meshElems(cEl)%neighbours(newFc,1)
			rayIterCt = rayIterCt + 1
		end do
	end subroutine traceSingleRay

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

	subroutine startRayInVolume(emEl,pL,pt,dir)
		integer :: elNodes(4)
		integer,intent(in) :: emEl
		real(8),intent(out) :: pL,pt(3),dir(3)

		elNodes = meshElems(emEl)%nodes
		ec = meshVerts(elNodes,:)
		pt = selTetraPoint(ec)
		dir = getRaySphDir()
		pL = getRayPathLength()
	end subroutine startRayInVolume

	subroutine startRayFromSurf(pRatio,emEl,emFc,pL,pt,dir)
		integer :: remNo,fcNodes(3),elNodes(4)
		integer,intent(out) :: emEl,emFc
		real(8) :: remVert(3),fcNorm(3)ec(4,3),fcVerts(3,3)
		real(8),intent(in) :: pRatio(:)
		real(8),intent(out) :: pL,pt(3),dir(3)

		call chooseEmittingFace(pRatio,emEl,emFc)
		elNodes = meshElems(emEl)%nodes
		ec = meshVerts(elNodes,:)
		fcNodes = getFaceNodes(emFc)
		fcVerts = ec(fcNodes,:)
		remNo = 10-sum(fcNodes)
		remVert = ec(remNo,:)
		fcNorm = getFaceNorm(fcVerts,remVert,.true.)
		pt = selTriPoint(fcVerts)
		dir = getFaceRayDir(fcVerts,fcNorm)
		pL = getRayPathLength()
	end subroutine startRayFromSurf

	subroutine chooseEmittingFace(pRatio,emEl,emFc)
		integer :: i,nSf
		integer,intent(out) :: emEl,emFc
		real(8) :: rs
		real(8),intent(in) :: pRatio(:)
		type(emissionSurface) :: emSf

		nSf = size(pRatio,1) + 1
		call random_number(rs)
		do i=1,nSf-1
			if(rs .lt. pRatio(i)) then
				emSf = rtEmSurfs(i)
				call getEmFace(emSf,emEl,emFc)
				exit
			elseif((rs.gt.pRatio(i)).and.(rs.lt.pRatio(i+1))) then
				emSf = rtEmSurfs(i+1)
				call getEmFace(emSf,emEl,emFc)
				exit
			end if
		end do
	end subroutine chooseEmittingFace

	subroutine getEmFace(emSf,emEl,emFc)
		integer :: i
		integer,intent(out) :: emEl,emFc
		real(8):: psi
		type(emissionSurface) :: emSf

		call random_number(psi)
	    if (psi > 0.5d0) then
	        do i=size(emSf%cuSumFcEmPow,1)-1,1,-1
	            if (emSf%cuSumFcEmPow(i) .lt. psi) exit
	        end do
	        emEl = emSf%emSurf%elNum(i+1)
	        emFc = emSf%emSurf%fcNum(i+1)
	    else           
	        do i=1,size(emSf%cuSumFcEmPow,1)
	            if (emSf%cuSumFcEmPow(i) .gt. psi) exit
	        end do
	        emEl = emSf%emSurf%elNum(i)
	        emFc = emSf%emSurf%fcNum(i)
	    end if
	end subroutine getEmFace

	subroutine createEmissionSurfaces(emSurfOriginalIds)
		integer :: i,j,nEmSf,currSurf,nEmFcs,fcNodes(3),elNodes(4)
		integer,intent(in) :: emSurfOriginalIds(:)
		real(8) :: fcEmPow,fcArea,fcCentT,fcNoTs(3)

		nEmSf = size(emSurfOriginalIds)
		if(.not.(allocated(rtEmSurfs))) then
			allocate(rtEmSurfs(nEmSf))
		end if
		do i=1,nEmSf
			currSurf = emSurfOriginalIds(i)
			rtEmSurfs(i)%emSurf = meshSurfs(currSurf)
			nEmFcs = meshSurfs(currSurf)%numFcs
			if(.not.(allocated(rtEmSurfs(i)%fcEmPower))) then
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
				fcEmPow = kappa*sigB*fcrArea*(fcCentT**4.0d0)
				emSfPow = emSfPow+fcEmPow
				rtEmSurfs(i)%cuSumFcEmPow(j) = emSfPow
			end do
			rtEmSurfs(i)%totEmPow = emSfPow
			rtEmSurfs(i)%cuSumFcEmPow=rtEmSurfs(i)%cuSumFcEmPow/emSfPow
		end do
	end subroutine createEmissionSurface

	function specularReflection(ec,fcNum,dirIn) result(dirOut)
		integer :: remNo,fcNodes(3)
		integer,intent(in) :: fcNum
		real(8) :: cosInc,fcVerts(3,3),remVert(3),fcNorm(3)
		real(8),intent(in) :: ec(4,3),dirIn(3)
		real(8),intent(out) :: dirOut

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

	function getRayPathLength() result(pathLength)
		real(8) :: randL,pathLength

		call random_number(randL)
		pathLength = (1.0d0/(kappa+sigma))*log(1.0d0/randL)
	end function getRayPathLength

	subroutine setMediumValues(valIn,valNameFlag)
		real(8),intent(in) :: valIn
		character(*),intent(in) :: valNameFlag

		if(valNameFlag .eq. "k") then
			kappa = valIn
		elseif(valNameFlag .eq. "s") then
			sigma = valIn
		else
			write(*,'(a)') "Property to be set not recognised."
			stop
		end if
	end subroutine setMediumValues

end module rtHelper
