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
	integer,parameter :: rtNumPFTable=10000
	integer :: rtNumRays,rtNumEmSfs,rtNumCTSfs,rtNumCQSfs,rtNumTrSfs,	&
	rtNumNpSfs,rtElemMinRays
	integer,allocatable :: rtElemAbs(:),rtElemSto(:),rtEmSfIds(:),		&
	rtConstTSfIds(:),rtConstQSfIds(:),rtTrSfIds(:),rtNpSfIds(:),		&
	rtWallInf(:)
	real(8) :: rtBeta,rtKappa,rtSigma,rtRefrInd,rtRefRayPow,rtAbsThr,	&
	rtReEmThr
	real(8),allocatable :: rtWallSrc(:),rtNodalSrc(:),rtPFTable(:)
	type(emissionSurface),allocatable :: rtEmSurfs(:)

	contains

	subroutine rtInit(mBinNum,mSfConstTs,mSfConstQs,mFileName)
		integer :: i
		integer,intent(in) :: mBinNum(3)
		character(*),intent(in) :: mFileName
		real(8),intent(in) :: mSfConstTs(:),mSfConstQs(:)

		call rtInitMesh(mFileName,mBinNum)
		if(rtNumCTSfs .gt. 0) then
			do i=1,size(rtEmSfIds,1)
				call setSurfaceConstTemperature(rtConstTSfIds(i),		&
				mSfConstTs(i))
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
		if(.not.(allocated(rtEmSfIds))) then
			allocate(rtEmSfIds(rtNumEmSfs))
		end if
		rtEmSfIds = (/rtConstTSfIds,rtConstQSfIds/)
		call CreateUniformEmissionSurfaces()
		call popLargeSphDiffMuTable()
	end subroutine rtInit

	subroutine rtInitMesh(mFileName,mBinNum)
		integer,intent(in) :: mBinNum(3)
		character(*),intent(in) :: mFileName
		meshFile = trim(adjustl(mFileName))
		meshNBins = mBinNum
		call readMesh()
		call getElementNeighbours()
		call populateSurfaceFaceAreas()
	end subroutine rtInitMesh

	subroutine traceFromSurf(pRatio)
		integer,parameter :: limoutpt=5
		integer,parameter :: nSfAbPtsFil=175,nSfEmPtsFil=197
		integer :: i,j,cEl,emEl,emFc,endEl,outPtCt,elNodes(4)
		real(8) :: pL,pt(3),dir(3),endPt(3),spFnVals(4)
		real(8),intent(in) :: pRatio(:)
		logical :: outPt
		character(*),parameter :: sfEmPtsFil=commResDir//"surfems.out",	&
								  sfAbPtsFil=commResDir//"surfpts.out"

		open(nSfAbPtsFil,file=sfAbPtsFil)
		open(nSfEmPtsFil,file=sfEmPtsFil)
		rtNodalSrc = 0.0d0
		outPtCt = 0
		do i=1,rtNumRays
			call startRayFromSurf(pRatio,emEl,emFc,pL,pt,dir)
			write(nSfEmPtsFil,'(6(f15.12,2x))') pt,dir
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
				rtElemAbs(endEl) = rtElemAbs(endEl) + 1
				write(nSfAbPtsFil,'(3(f15.12,2x))') endPt
				call shapeFunctionsAtPoint(endEl,endPt,spFnVals)
				elNodes = meshElems(endEl)%nodes
				rtNodalSrc(elNodes) = rtNodalSrc(elNodes) + 			&
				rtRefRayPow*spFnVals
			end if
		end do
		close(nSfEmPtsFil)
		close(nSfAbPtsFil)
		open(975,file="../obj/tempRadSrc.out")
			write(975,'(f20.13)') rtNodalSrc
		close(975)
		rtWallSrc = rtNodalSrc
	end subroutine traceFromSurf

	subroutine traceFromSurfLED(pRatio)
		integer,parameter :: limoutpt=5
		integer,parameter :: nSfAbPtsFil=175,nSfEmPtsFil=197
		integer :: i,j,cEl,emEl,emFc,endEl,outPtCt,elNodes(4),cb,cy,reEmCt,reEmDropCt
		real(8) :: pL,rAbs,rReEm,pt(3),dir(3),endPt(3),dirOut(3),ptScr(3),spFnVals(4)
		real(8),intent(in) :: pRatio(:)
		logical :: outPt,blue,trans
		character(*),parameter :: sfEmPtsFil=commResDir//"surfems.out",	&
								  sfAbPtsFil=commResDir//"surfpts.out"

		open(nSfAbPtsFil,file=sfAbPtsFil)
		open(nSfEmPtsFil,file=sfEmPtsFil)
		rtNodalSrc = 0.0d0
		outPtCt = 0
		cb = 1
		cy = 2
		open(1975,file="../results/tempOutPts.out")
		open(1995,file="../results/screenpts_yellow.out")
		open(1996,file="../results/screenpts_blue.out")
		open(1985,file="../results/scatterDirs.out")
		open(1965,file="../results/reEmDrop.out")
		open(1955,file="../results/ptsBackplate.out")
		reEmDropCt = 0
		do i=1,rtNumRays
!			write(*,*) "ray: ", i
			if(mod(i,rtNumRays/10) .eq. 0) write(*,*) "Nrays: ",i
			call startRayFromSurf(pRatio,emEl,emFc,pL,pt,dir)
			write(nSfEmPtsFil,'(6(f15.12,2x))') pt,dir
			blue = .true.
			cEl = emEl
			if(i .eq. 13848) then
				write(*,*) "First dropped: ", i
			end if
			do while(blue)
				call traceSingleRayWithTrans(pt,dir,pL,cEl,outPt,endEl,	&
				endPt,trans,dirOut)
				if(outPt) then
					outPtCt = outPtCt + 1
					if(outptct .ge. limoutpt) then
						write(*,*)"Count of dropped points reached limit."
						stop
					end if
					blue = .false.
				end if
				if(endEl .ne. 0) then
!					write(*,*) "Endel not zero"
					call random_number(rAbs)
					if(rAbs .lt. rtAbsThr) then
!						write(*,*) "Absorbed"
						call random_number(rReEm)
						if(rReEm .lt. rtReEmThr) then
!							write(*,*) "Reemitted"
							reEmCt = 0
							do while(endEl .ne. 0)
								reEmCt = reEmCt+1
								if(reEmCt .gt. MEGA) then
									write(*,*) "Ray dropped in reemission loop"
									write(1965,'(2(i8,2x),2(3f12.9,2x))') i,endEl,pt,dir
									reEmDropCt = reEmDropCt + 1
									if(reEmDropCt .gt. 10) then
										write(*,*) "Too many rays dropped in reEm loop"
										stop
									end if
									exit
								end if
								dir = getRaySphDir()
								pL = getRayPathLength()
								pt = endPt
								cEl = endEl
								call traceSingleRayWithTrans(pt,dir,pL,	&
								cEl,outPt,endEl,endPt,trans,dirOut)
!								write(*,*) "Stuck in re-emission loop"
							end do
							if(trans) then
								write(1975,'(6(f12.9,2x),i2)') pt,dirOut,cy
								if(dir(3) .gt. 0.d0) then
									ptScr = ((3.d0-pt(3))/dirOut(3))*dirOut + pt
									write(1995,'(3(f16.9,2x),i2)') ptScr, cy
								end if
							else
								write(1955,'(2(i8,2x),6(f12.9,2x))') i,cEL,pt,dir
							end if
						else
							rtElemAbs(endEl) = rtElemAbs(endEl) + 1
							write(nSfAbPtsFil,'(3(f15.12,2x))') endPt
							call shapeFunctionsAtPoint(endEl,endPt,		&
							spFnVals)
							elNodes = meshElems(endEl)%nodes
							rtNodalSrc(elNodes) = rtNodalSrc(elNodes) + &
							rtRefRayPow*spFnVals
						end if
						blue = .false.
					else
						pL = getRayPathLength()
						pt = endPt
						cEl = endEl
						dir = scatterRayLargeSphereCloud()
						write(1985,'(3(f15.12,2x))') dir
					end if
				else
					if(trans) then
						write(1975,'(6(f12.9,2x),i2)') pt,dirOut,cb
						if(dir(3) .gt. 0.d0) then
							ptScr = ((3.d0-pt(3))/dirOut(3))*dirOut + pt
							write(1996,'(3(f12.9,2x),i2)') ptScr,cb
						end if
					end if
					blue = .false.
				end if
			end do
		end do
		close(1995)
		close(1996)
		close(1985)
		close(1965)
		close(1955)
		close(1975)
		close(nSfEmPtsFil)
		close(nSfAbPtsFil)
		open(975,file="../data/tempRadSrc.out")
			write(975,'(f20.13)') rtNodalSrc
		close(975)
		write(*,*) "Rays dropped in re-emission loop: ", reEmDropCt
		write(*,*) "Out points in run: ", outPtCt
		rtWallSrc = rtNodalSrc
	end subroutine traceFromSurfLED

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

	subroutine traceSingleRayWithTrans(pt,dir,pL,cEl,outPt,endEl,endPt,	&
	trans,dirOut)
		integer :: i,rayIterCt,chCt,trCt,newFc,nEmSfs,nhbrFc,elNodes(4)
		integer,intent(inout) :: cEl
		integer,intent(out) :: endEl
		integer,allocatable :: emSfIds(:)
		real(8) :: lTrav,lToFc,newDir(3),ec(4,3)
		real(8),intent(in) :: pL
		real(8),intent(out) :: endPt(3),dirOut(3)
		real(8),intent(inout) :: pt(3),dir(3)
		logical :: inFc
		logical,intent(out) :: outPt,trans

		rayIterCt = 0
		endEl = 0
		endPt = 0.0d0
		lTrav = 0.0d0
		dirOut = 0.0d0
		outPt = .false.
		trans = .false.
		do while(lTrav.lt.pL)
			elNodes = meshElems(cEl)%nodes
			ec = meshVerts(elNodes,:)
			pt = pt + PICO*dir
			lTrav = lTrav+PICO
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
				chCt = count(-nhbrFc==rtEmSfIds)
				trCt = count(-nhbrFc==rtTrSfIds)
				if((chCt == 0) .and. (trCt == 0)) then
					newDir = specularReflection(ec,newFc,dir)
					dir = newDir
				end if
				if(trCt .gt. 0) then
					call transSurface(ec,newFc,dir,trans,newDir)
					if(trans) then
						lTrav = MEGA
						dirOut = newDir
						exit
					end if
					dir = newDir
				end if
				if(chCt .gt. 0) then
!					endPt = pt
!					dirOut = dir
!					endEl = cEl
					lTrav = MEGA
					exit
				end if
			end if
			cEl = meshElems(cEl)%neighbours(newFc,1)
			rayIterCt = rayIterCt + 1
		end do
	end subroutine traceSingleRayWithTrans

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
			pt = pt + PICO*dir
			lTrav = lTrav+PICO
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

	subroutine traceOut(pRatio)
		integer,parameter :: limoutpt=5
		integer,parameter :: nAbsPtsFil=124,nEmPtsFil=248
		integer :: i,j,cEl,emEl,emFc,endEl,outPtCt,elNodes(4)
		real(8) :: pL,pt(3),dir(3),endPt(3),spFnVals(4)
		real(8),intent(in) :: pRatio(:)
		real(8),allocatable :: nodalSources(:)
		logical :: outPt
		character(*),parameter :: emPtsFil=commResDir//"empts.out",	&
								  absPtsFil=commResDir//"abspts.out"

		open(nAbsPtsFil,file=absPtsFil)
		open(nEmPtsFil,file=emPtsFil)
		allocate(nodalSources(meshNumNodes))
		nodalSources = 0.0d0
		outPtCt = 0
		do i=1,rtNumRays
			call startRayFromSurf(pRatio,emEl,emFc,pL,pt,dir)
			write(nEmPtsFil,'(6(f15.12,2x))') pt,dir
			cEl = emEl
			call traceSingleRay(pt,dir,pL,cEl,outPt,endEl,endPt)
			if(outPt) then
				outPtCt = outPtCt + 1
				if(outptct .ge. limoutpt) then
					write(*,*)"Count of dropped points reached limit."
					stop
				end if
			end if
			do while(endEl .ne. 0)
				rtElemAbs(endEl) = rtElemAbs(endEl) + 1
				write(nAbsPtsFil,'(3(f15.12,2x))') endPt
				call shapeFunctionsAtPoint(endEl,endPt,spFnVals)
				elNodes = meshElems(endEl)%nodes
				nodalSources(elNodes) = nodalSources(elNodes) + 		&
				rtRefRayPow*spFnVals
				dir = getRaySphDir()
				pL = getRayPathLength()
				pt = endPt
				cEl = endEl
				call traceSingleRay(pt,dir,pL,cEl,outPt,endEl,endPt)
			end do
		end do
		close(nEmPtsFil)
		close(nAbsPtsFil)
		open(975,file="../obj/tempRadSrc.out")
			write(975,'(f20.13)') nodalSources
		close(975)
		call setMeshNodalValues(nodalSources,"S")
		rtWallSrc = nodalSources
	end subroutine traceOut

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
		real(8) :: ec(4,3)
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
		real(8) :: remVert(3),fcNorm(3),ec(4,3),fcVerts(3,3)
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

	subroutine createNonUniformEmissionSurface()
!	Placeholder for future use in case of non-uniform surface 
!	temperatures
	end subroutine createNonUniformEmissionSurface

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

	subroutine transSurface(ec,fcNum,dirIn,trans,dirOut)
		integer,intent(in) :: fcNum
		integer :: remNo,fcNodes(3)
		real(8),parameter :: n2 = 1.d0	! Air is the second medium
		real(8) :: rIndRatio,thi,trRatio,cosInc,sinTht,remVert(3),	&
		fcNorm(3),fcVerts(3,3)
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
		rIndRatio = (rtRefrInd/n2)
		trRatio = rIndRatio*sqrt(1-cosInc**2.d0)
		if(trRatio .gt. 1.d0) then
			dirOut = specularReflection(ec,fcNum,dirIn)
		else
			trans = .true.
			sinTht = rIndRatio*sin(acos(cosInc))
			dirOut = (rIndRatio*cosInc - sqrt(1.d0-sinTht**2.d0))*fcNorm
			dirOut = dirOut + rIndRatio*dirIn
		end if
	end subroutine transSurface

	function getRayPathLength() result(pathLength)
		real(8) :: randL,pathLength

		call random_number(randL)
		pathLength = (1.0d0/(rtKappa+rtSigma))*log(1.0d0/randL)
	end function getRayPathLength

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
		integer :: elNodes(4)
		integer,intent(in) :: elNum
		integer,intent(out) :: elNumRays
		real(8) :: elVol,elCentTemp,elEmPow,elNoTemps(4)
		real(8),intent(out) :: elRayPow

		elNodes = meshElems(elNum)%nodes
		elNoTemps = meshTemperatures(elNodes)
		elVol = getElementVolume(elNum)
		elCentTemp = sum(elNoTemps)/4.0d0
		elEmPow = 4.0d0*rtKappa*sigB*(elCentTemp**4.0d0)*elVol
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

		ratio = (rMu-muMin)/(muMax-muMin)
		interLoc = floor(ratio*rtNumPFTable)
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

	subroutine popLargeSphDiffMuTable()
		integer :: i,nSplSteps,lV
		real(8) :: j,h,m,a,b,c,d,r,stepSize
		real(8),allocatable :: x(:),y(:),yy(:)

		nSplSteps = 10
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
		stepSize = 1.0d0/rtNumPFTable
		i = 0
		do m=0.0d0,0.99999d0,stepSize
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

	function getFaceRayDir(fcVerts,fcNorm) result(dir)
		real(8),intent(in) :: fcVerts(3),fcNorm(3)
		real(8):: th,ph,dir(3)

		call random_number(th)
		th = asin(sqrt(th))
		call random_number(ph)
    	ph = 2.0d0*pi*ph
		dir = getDirectionCoords(th,ph)

		if(dot_product(fcNorm,dir) .lt. 0) then
			th = pi-th
			dir = getDirectionCoords(th,ph)
		end if

	end function getFaceRayDir

	function getDirectionCoords(th,ph) result(dir)
		real(8) :: th,ph,dir(3)

		dir(1) = sin(th)*cos(ph)
		dir(2) = sin(th)*sin(ph)
		dir(3) = cos(th)
		dir = dir/norm2(dir)
	end function getDirectionCoords

end module rt
