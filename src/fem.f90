include 'mkl_pardiso.f90'	! Added for ParDiSo

module fem

	use utilities
	use mesh

	implicit none

!	Declare all relevant types
	type nodeRow
		integer,allocatable :: col(:)
		real(8),allocatable :: val(:)
	end type nodeRow

!	Module variables
	integer :: femSolverMaxIter
	integer,allocatable :: femBCs(:),femStRowPtr(:),femStColPtr(:),		&
	femCpRowPtr(:),femCpColPtr(:)
	real(8) :: femAmbT
	real(8),allocatable :: femKs(:),femBVals(:),femRhos(:),femCaps(:),	&
	femSrc(:),femTvals(:),femStVals(:),femCpVals(:)
	logical :: femBCsKnown,femSysAssembled,femTransient,femUnifGen
	character(72) :: femByFile="boundaries",femResFile="res",			&
	femSolverType="BiCGStab",femTrScheme="rk"
	type(nodeRow),allocatable :: femSt(:),femCp(:)

	contains

	subroutine runFem(mFileName)
		integer :: pbDatFileNum
		character(*),intent(in) :: mFileName

!		if(femSysAssembled) then
!			if(norm2(meshSources-femSrc) .lt. NANO) then
!				write(*,*)"System already solved, no changes detected."
!				stop
!			end if
!			femSrc = meshSources
!			call setupFinalEquations()
!			call solveFemSystem()
!			call writeVtkResults()
!		else
			if(.not.(allocated(meshElems))) then
				call femInitMesh(mFileName)
			end if
			if(.not.(femBCsKnown)) then
				call getBoundaryConditions()
			end if
			call assembleFemSystem()
!			femSysAssembled = .true.
			call solveFemSystem()
			call setMeshNodalValues(femTvals,"T")
			call writeVtkResults()
			call writeNodalResults()
!			deallocate(femSt)
			deallocate(femSrc)
			deallocate(femTvals)
			deallocate(femStVals)
			deallocate(femStRowPtr)
			deallocate(femStColPtr)
!			call writeElementNeighbourhoodData()
!		end if
	end subroutine runFem

	subroutine femInitMesh(mFileName)
		character(*),intent(in) :: mFileName

		meshFile = trim(adjustl(mFileName))
		call readMesh()
	end subroutine femInitMesh

	subroutine assembleFemSystem()
		integer :: i,j,k,elDom,elNodes(4),elBys(4)
		real(8) :: elK,elC,elRho,elVol,elByTs(4),elBySrc(4),elGenSrc(4),&
		elSt(4,4),elBySt(4,4),elCp(4,4)
		type(nodeRow) :: stElem(4),cpElem(4)

		allocate(femSt(meshNumNodes))
		allocate(femSrc(meshNumNodes))
		allocate(femTvals(meshNumNodes))
		femSrc = 0.0d0
		femTvals = 0.0d0
		if(femTransient) then
			allocate(femCp(meshNumNodes))
		end if
		do i=1,meshNumElems
			elNodes = meshElems(i)%nodes
			elDom = meshElems(i)%domain
			elK = femKs(elDom)
			call getElementUnitStiffness(i,elK,elVol,elSt)
!			elBys = meshElems(i)%neighbours(:,2)
!			if(any(elBys<0)) then
				call boundaryConditions(i,elBys,elBySt,elByTs,elBySrc)
				elSt = elSt + elBySt
				call addToGlobalTemperature(elNodes,elByTs)
				call addToGlobalSource(elNodes,elBySrc)
!				femSrc(elNodes) = elBySrc
!				femTvals(elNodes) = elByTs
!			end if
			
			stElem = femSt(elNodes)
			call assembleNoderows(elNodes,elSt,stElem)
			femSt(elNodes) = stElem
			if(femTransient) then
				elC = femCaps(elDom)
				elRho = femRhos(elDom)
				call getElementCapacitance(elC,elRho,elVol,elCp)
				cpElem = femCp(elNodes)
				call assembleNoderows(elNodes,elCp,cpElem)
				femCp(elNodes) = cpElem
			end if
			if(femUnifGen) then
				call getUniformSource(100.0d0,elVol,elGenSrc)
				call addToGlobalSource(elNodes,elGenSrc)
			end if
		end do

		if(femTransient) then
			call collapseNodeRows("st")
			deallocate(femSt)
			call collapseNodeRows("cp")
			deallocate(femCp)
			call getInitialCondition()
		end if
		if(allocated(meshSources)) then
			if(norm2(meshSources) .gt. MICRO) then
				femSrc = femSrc + meshSources
			end if
		end if
		call setupFinalEquations()
		call collapseNodeRows("st")
!		open(111,file="../obj/tempSrc.out")
!		write(111,'(a)')"femStRowPtr:"
!		write(111,'(i8)')femStRowPtr
!		write(111,*)
!		write(111,'(a)')"femStColPtr:"
!		write(111,'(i8)')femStColPtr
!		write(111,*)
!		write(111,'(a)')"femStVals:"
!		write(111,'(f9.5)')femStVals
!		write(111,'(a)')"femTemperatures:"
!		write(111,'(f9.5)')femTvals
!		close(111)
		deallocate(femSt)
	end subroutine assembleFemSystem

	subroutine setupFinalEquations()
		integer :: i

		do i=1,meshNumNodes
			if(femTvals(i) .ne. 0.0d0) then
				if(allocated(femSt(i)%col)) then
					deallocate(femSt(i)%val)
					deallocate(femSt(i)%col)
				end if
				allocate(femSt(i)%col(1))
				allocate(femSt(i)%val(1))
				femSt(i)%col(1) = i
				femSt(i)%val(1) = 1.0d0
				femSrc(i) = femTvals(i)
			end if
		end do
	end subroutine setupFinalEquations

!-----------------------------------------------------------------------!
!	Routines for boundary handling
!-----------------------------------------------------------------------!
    subroutine getBoundaryConditions()
        integer,parameter :: byFilNum=401
		integer :: i
		character(72) :: femByFileName

		allocate(femKs(meshNumDoms))
		allocate(femRhos(meshNumDoms))
		allocate(femCaps(meshNumDoms))
		allocate(femBCs(meshNumSurfs))
		allocate(femBVals(meshNumSurfs))
        femBVals = 0.0d0
        femBCs = 4
		femByFileName = commDatDir//trim(adjustl(femByFile))//commDatExt
        open(byFilNum,file=femByFileName)
        read(byFilNum,*)
        read(byFilNum,*)
        do i=1,meshNumDoms
			read(byFilNum,*)
            read(byFilNum,*) femKs(i)
        end do
		do i=1,meshNumDoms
			read(byFilNum,*)
			read(byFilNum,*) femRhos(i)
			read(byFilNum,*)
			read(byFilNum,*) femCaps(i)
		end do
		read(byFilNum,*)
        read(byFilNum,*) femBCs
        read(byFilNum,*)
        read(byFilNum,*) femBVals
        read(byFilNum,*)
        read(byFilNum,*) femAmbT
		close(byFilNum)
		femBCsKnown = .true.
    end subroutine getBoundaryConditions

	subroutine boundaryConditions(elNum,elBys,elBySt,elByTs,elBySrc)
		integer :: i,j,bSfId,bSfType,fcNodes(3),boSfs(4),elNodes(4) 
		integer,intent(in) :: elNum,elBys(4)
		real(8) :: bSfVal,ec(4,3),elQ(4),htA(4)
		real(8) :: elByTs(4),elBySrc(4),elBySt(4,4)

		boSfs = meshElems(elNum)%neighbours(:,2)
		elNodes = meshElems(elNum)%nodes
		ec = meshVerts(elNodes,:)
		elByTs = 0.d0
		elBySrc = 0.d0
		elBySt = 0.d0
		elQ = 0.d0
		htA = 0.d0
		do i=1,4
			if(boSfs(i) < 0) then
				bSfId = -boSfs(i)
				bSfType = femBCs(bSfId)
				bSfVal = femBVals(bSfId)
				fcNodes = getFaceNodes(i)
				if(abs(bSfVal).ge.PICO) then
					if(bSfType == 1) then
						elByTs(fcNodes) = bSfVal
					elseif(bSfType == 2) then
						elQ = getFluxSource(ec,fcNodes,bSfVal)
						elBySrc = elBySrc - elQ
					elseif(bSfType == 3) then
						call getConvectiveSource(ec,fcNodes,bSfVal,htA,	&
						elBySt)
						elBySrc = elBySrc + htA
					elseif(bSfType == 4) then
						continue
					else
						write(*,*)"Unrecognised boundary type."
					end if
				else
					continue
				end if
			else
				continue
			end if
		end do
	end subroutine boundaryConditions

	function getFluxSource(ec,fcNodes,bSfVal) result(elQ)
		integer,intent(in) :: fcNodes(3)
		real(8) :: fcArea,fc(3,3),elQ(4)
		real(8),intent(in) :: bSfVal,ec(4,3)

		elQ = 0.0d0
		fc = ec(fcNodes,:)
		fcArea = triangleArea(fc)
		elQ(fcnodes) = (1.0d0/6.0d0)*(2.0d0*fcarea)*bSfVal
	end function getFluxSource

	subroutine getConvectiveSource(ec,fcNodes,bSfVal,htA,elBySt)
		integer,intent(in) :: fcNodes(3)
		real(8) :: fcArea,fc(3,3),spFnAreaInt(4,4)
		real(8),intent(in) :: bSfVal,ec(4,3)
		real(8),intent(out) :: htA(4),elBySt(4,4)

		elBySt = 0.0d0
		fc = ec(fcNodes,:)
		fcArea = triangleArea(fc)
		htA(fcnodes) = (1.0d0/6.0d0)*(2.0d0*fcArea)*bSfVal*femAmbT
		spFnAreaInt = areaIntShapeFuncs(fcnodes)
		elBySt = (2.0d0*fcArea/24.0d0)*bSfVal*spFnAreaInt
	end subroutine getConvectiveSource

	subroutine getUniformSource(uniGenVal,elVol,elGenSrc)
		real(8),intent(in) :: uniGenVal,elVol
		real(8),intent(out) :: elGenSrc(4)

		elGenSrc = 1.0d0
		elGenSrc = uniGenVal*(elvol/4.0d0)*elGenSrc
	end subroutine getUniformSource

	subroutine addToGlobalTemperature(elNodes,elByTs)
		integer,intent(in) :: elNodes(4)
		real(8),intent(in) :: elByTs(4)

		if(any(abs(elByTs).ge.PICO)) then
			femTvals(elNodes) = elByTs
		end if
	end subroutine addToGlobalTemperature

	subroutine addToGlobalSource(elNodes,elBySrc)
		integer :: i
		integer,intent(in) :: elNodes(4)
		real(8),intent(in) :: elBySrc(4)

		do i=1,4
			femSrc(elnodes(i)) = femSrc(elnodes(i))+elBySrc(i)
		end do
	end subroutine addToGlobalSource
!-----------------------------------------------------------------------!
!	End of boundary handling routines
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
!	Routines required to assemble the sparse global stiffness matrix
!-----------------------------------------------------------------------!

	subroutine assembleNoderows(elNodes,elSt,nodeRowElem)
		integer :: i,sz
		integer,intent(in) :: elnodes(4)
		integer,allocatable :: tempCol(:)
		real(8) :: elSt(4,4)
		real(8),allocatable :: tempVal(:)
		type(nodeRow),intent(inout) :: nodeRowElem(4)

		do i=1,4
			if(.not.(allocated(nodeRowElem(i)%col))) then
				allocate(nodeRowElem(i)%col(4))
				allocate(nodeRowElem(i)%val(4))
				nodeRowElem(i)%col = elnodes
				nodeRowElem(i)%val = elSt(i,:)
			else
				sz = size(nodeRowElem(i)%col,1)
				allocate(tempcol(sz+4))
				allocate(tempval(sz+4))
				tempcol = (/nodeRowElem(i)%col,elnodes/)
				tempval = (/nodeRowElem(i)%val,elSt(i,:)/)
				call move_alloc(tempCol,nodeRowElem(i)%col)
				call move_alloc(tempVal,nodeRowElem(i)%val)
				if(allocated(tempCol)) deallocate(tempCol)
				if(allocated(tempVal)) deallocate(tempVal)
			end if
		end do

	end subroutine assembleNoderows

	subroutine collapseNodeRows(whichMat)
		integer :: i,indStart,indEnd,numUnique,numVal,sz
		integer,allocatable :: tempcol(:),ind(:),rowPtr(:),colPtr(:)
		real(8),allocatable :: tempval(:),vals(:)
		character(*),intent(in) :: whichMat
		type(nodeRow),allocatable :: sysNodeRows(:)

		allocate(sysNodeRows(meshNumNodes))
		if(whichMat == "st") then
			sysNodeRows = femSt
		elseif(whichMat == "cp") then
			sysNodeRows = femCp
		else
			write(*,*)"Unrecognised call to collapseNodeRows."
			stop
		end if
		allocate(rowPtr(meshNumNodes+1))
		rowPtr(1) = 1
		numVal = 0
		do i=1,meshNumNodes
			sz = size(sysNodeRows(i)%col,1)
			allocate(ind(sz))
			call indexedSortInteger(sysNodeRows(i)%col,ind)
			sysNodeRows(i)%val = sysNodeRows(i)%val(ind)
			deallocate(ind)
			call addNodeEntries(sysNodeRows(i)%col,sysNodeRows(i)%val,	&
			numUnique)
			numVal = numVal+numUnique
			rowPtr(i+1) = rowPtr(i)+numUnique
		end do

		allocate(vals(numVal))
		allocate(colPtr(numVal))
		do i=1,meshNumNodes
			indStart = rowPtr(i)
			indEnd = (rowPtr(i+1)-1)
			vals(indStart:indEnd) = sysNodeRows(i)%val
			colPtr(indStart:indEnd) = sysNodeRows(i)%col
			deallocate(sysNodeRows(i)%col)
			deallocate(sysNodeRows(i)%val)
		end do

		if(whichMat == "st") then
			allocate(femStRowPtr(meshNumNodes+1))
			allocate(femStColPtr(numVal))
			allocate(femStVals(numVal))
			femStRowPtr = rowPtr
			femStColPtr = colPtr
			femStVals = vals
		elseif(whichMat == "cp") then
			allocate(femCpRowPtr(meshNumNodes+1))
			allocate(femCpColPtr(numVal))
			allocate(femCpVals(numVal))
			femCpRowPtr = rowPtr
			femCpColPtr = colPtr
			femCpVals = vals
		end if

		deallocate(rowPtr)
		deallocate(colPtr)
		deallocate(vals)
		deallocate(sysNodeRows)
	end subroutine collapseNodeRows

	subroutine addNodeEntries(nodeCols,nodeVals,numUnique)
		integer :: i,j,k
		integer,intent(out) :: numUnique
		integer,allocatable :: tempCols(:),starts(:)
		integer,allocatable,intent(inout) :: nodeCols(:)
		real(8),allocatable :: tempVals(:)
		real(8),allocatable,intent(inout) :: nodeVals(:)

		if(size(nodeCols,1)==1) then
			numUnique = 1
			return
		end if

		call getNodeRepeatEntries(nodeCols,starts,numUnique)
		allocate(tempCols(numUnique))
		allocate(tempVals(numUnique))
		tempCols = nodeCols(starts)

		do i=1,numUnique
			tempVals(i) = sum(nodeVals(starts(i):starts(i+1)-1))
		end do
		call move_alloc(tempCols,nodeCols)
		call move_alloc(tempVals,nodeVals)
	end subroutine addNodeEntries

	subroutine getNodeRepeatEntries(nodeCols,starts,numUnique)
		integer :: i,j,ind
		integer,intent(in) :: nodeCols(:)
		integer,intent(out) :: numUnique
		integer,allocatable :: temp(:)
		integer,allocatable,intent(out) :: starts(:)
		integer,dimension(1) :: ml

		allocate(starts(size(nodeCols,1)))
		starts = 0
		ind=1
		starts(1) = 1
		i = 2
		do while(nodeCols(ind).lt.maxval(nodeCols))
			ml = minloc(nodeCols,nodeCols.gt.nodeCols(ind))
			starts(i) = ml(1)
			ind = ml(1)
			i = i+1
		end do
		allocate(temp(i))
		temp = (/starts(1:i-1),size(nodeCols,1)+1/)
		call move_alloc(temp,starts)
		numUnique = i-1
	end subroutine getNodeRepeatEntries
!-----------------------------------------------------------------------!
!	End of the assembly routines
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
!	Routines to solve transient heat transfer problems
!-----------------------------------------------------------------------!

	subroutine getElementCapacitance(elC,elRho,elVol,elCp)
		integer :: i
		real(8),intent(in) :: elC,elVol,elRho
		real(8),intent(out) :: elCp(4,4)

		elCp = 1.d0;
		do i=1,4
			elCp(i,i) = 2.d0
		end do
		elCp = elC*elRho*(elVol/20.d0)*elCp
	end subroutine getElementCapacitance

	subroutine getInitialCondition()

		femTvals = 100.0d0			! Only a placeholder, to be replaced
	end subroutine getInitialCondition

	subroutine solveTransientSystem()
!		integer :: nDoms,iaSt(:),jaSt(:),iaCp(:),jaCp(:),doElems(:),&
!		connTab(:,:)
!		real(8) :: sySt(:),syCp(:),sySrc(:),syInit(:),noverts(:,:)
!		logical :: useRK
!		character(*) :: resFilePre

		if(femTrScheme == "rk") then
			call transientRK()
		elseif(femTrScheme == "fd") then
			call transientFD()
		else
			write(*,*)"Transient solver type unrecognised."
			write(*,*)"Using default Finite Difference solver."
			call transientFD()
		end if
	end subroutine solveTransientSystem

	subroutine transientFD()
!		integer,parameter :: trfileno=888
!		integer :: i,j,fno,n,nv,ntstep,nDoms,iter,iaSt(:),jaSt(:),	&
!		iaCp(:),jaCp(:),doElems(:),connTab(:,:)
!		real(8) :: theta,tstep,tfinal,sySt(:),syCp(:),sySrc(:),		&
!		syInit(:),noVerts(:,:)
!		real(8),allocatable :: CKLhs(:),CKRhs(:),FRhs(:),Tnew(:),	&
!		kTRhs(:),initGuess(:)
!		character(*),parameter ::fres="res",fext=".vtk"
!		character(len=100) :: fName,sysCall
!		character(*) :: resFilePre

!		n = size(sySrc,1)
!		allocate(FRhs(n))
!		allocate(kTRhs(n))
!		allocate(initGuess(n))
!		nv = size(sySt,1)
!		allocate(CKLhs(nv))
!		allocate(CKRhs(nv))

!		tstep = 0.001d0
!		tfinal= 2.d0
!		ntstep = tfinal/tstep + 1
!		theta = 1.d0
!		CKLhs = theta*tstep*sySt + syCp
!		CKRhs = syCp - (1.d0-theta)*tstep*sySt
!		initGuess = 0.d0

!		do i=1,ntstep
!			FRhs = theta*sySrc + (1.d0-theta)*sySrc
!			call mkl_dcsrgemv("N",n,CKRhs,iaSt,jaSt,syInit,kTRhs)
!			FRhs = tstep*FRhs + kTRhs
!			call bicgstab(CKLhs,iaCp,jaCp,FRhs,100000,initGuess,	&
!			Tnew,iter)
!			syInit = Tnew
!			deallocate(Tnew)
!			if(mod(i,10).eq.0) then
!				fno = (i+1)/10
!				write(fname,*) fno
!				fName = trim(adjustl(resFilePre))//trim(adjustl(fname))
!				fName = trim(adjustl(fName))//fext
!				call writeresultsvtk(noVerts,connTab,nDoms,doElems,		&
!				syInit,trim(adjustl(fname)))
!			end if
!		end do
		
	end subroutine transientFD

	subroutine transientRK()
!		integer,parameter :: trfileno=888
!		integer :: i,j,k,n,ntstep,nDoms,fno,iaSt(:),jaSt(:),iaCp(:),&
!		jaCp(:),doElems(:),connTab(:,:)
!		real(8) :: tstep,tfinal,sySt(:),syCp(:),sySrc(:),syInit(:),	&
!		noVerts(:,:)
!		real(8),allocatable :: kTprod(:),rhs(:),diffT(:),newT(:),	&
!		k1(:),k2(:),k3(:),k4(:)
!		character(*),parameter :: fres="res",fext=".vtk"
!		character(len=100) :: fName,sysCall
!		character(*) :: resFilePre

!		n = size(sySrc,1)
!		allocate(kTprod(n))
!		allocate(rhs(n))
!		allocate(newT(n))

!		tstep = 0.0001
!		tfinal = 2
!		ntstep = tfinal/tstep + 1

!		do i = 1,ntstep
!			call mkl_dcsrgemv("N",n,sySt,iaSt,jaSt,syInit,kTprod)
!			rhs = sySrc - kTprod
!			call timegradient(syCp,iaCp,jaCp,rhs,k1)
!			newT = syInit + (tstep/2.d0)*k1
!			call mkl_dcsrgemv("N",n,sySt,iaSt,jaSt,newT,kTprod)
!			rhs = sySrc - kTprod
!			call timegradient(syCp,iaCp,jaCp,rhs,k2)
!			newT = syInit + (tstep/2.d0)*k2
!			call mkl_dcsrgemv("N",n,sySt,iaSt,jaSt,newT,kTprod)
!			rhs = sySrc - kTprod
!			call timegradient(syCp,iaCp,jaCp,rhs,k3)
!			newT = syInit + tstep*k3
!			call mkl_dcsrgemv("N",n,sySt,iaSt,jaSt,newT,kTprod)
!			rhs = sySrc - kTprod
!			call timegradient(syCp,iaCp,jaCp,rhs,k4)
!			syInit = syInit + (tstep/6)*(k1 + 2*k2 + 2*k3 + k4)
!			if(allocated(k1)) deallocate(k1)
!			if(allocated(k2)) deallocate(k2)
!			if(allocated(k3)) deallocate(k3)
!			if(allocated(k4)) deallocate(k4)
!			fno = (i+1)/10
!			write(fname,*) fno
!			fName = trim(adjustl(resFilePre))//trim(adjustl(fname))
!			fName = trim(adjustl(fName))//fext
!			call writeresultsvtk(noVerts,connTab,nDoms,doElems,		&
!			syInit,trim(adjustl(fname)))

!		end do

	end subroutine transientRK

!-------------------------------------------------------------------
!	Subroutine to get time gradient at each initial step of RK
!-------------------------------------------------------------------

	subroutine timegradient(syCp,iaCp,jaCp,rhs,diffT)
		integer :: n,iter,iaCp(:),jaCp(:)
		real(8) :: syCp(:),rhs(:)
		real(8),allocatable :: initGuess(:),diffT(:)

		n = size(rhs,1)
		allocate(initGuess(n))
		initGuess = 0.d0

!		call bicgstab(syCp,iaCp,jaCp,rhs,100000,initGuess,diffT,iter)
		diffT = 0.d0	! Placeholder

	end subroutine timegradient

!-----------------------------------------------------------------------!
!	End of transient solver routines
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
!	Different solver routines
!-----------------------------------------------------------------------!

	subroutine solveFemSystem()
		integer :: iter
		integer,allocatable :: ia(:),ja(:)
		real(8),allocatable :: b(:),x(:),initGuess(:),acsr(:)

		allocate(acsr(size(femStVals,1)))
		acsr = femStVals
		allocate(ia(size(femStRowPtr,1)))
		ia = femStRowPtr
		allocate(ja(size(femStColPtr)))
		ja = femStColPtr
		allocate(b(meshNumNodes))
		b = femSrc
		allocate(initGuess(meshNumNodes))
		initGuess = 0.0d0
		if(femSolverType == "BiCGStab") then
			call solveBiCGStab(acsr,ia,ja,b,initGuess,x,iter)
		elseif(femSolverType == "MinRes") then
			call solveMinRes(acsr,ia,ja,b,initGuess,x,iter)
		elseif(femSolverType == "ParDiSo") then
			call solveParDiSo(acsr,ia,ja,b,x)
		else
			write(*,*)"Solver type not recognised."
			write(*,*)"Using default BiCGStab solver."
			call solveBiCGStab(acsr,ia,ja,b,initGuess,x,iter)
		end if
		femTvals = x
	end subroutine solveFemSystem

	subroutine solveMinRes(acsr,ia,ja,b,initGuess,x,iter)
		integer :: i,j,k
		integer,intent(in) :: ia(:),ja(:)
		integer,intent(out) :: iter
		real(8),parameter :: cc=1e-7
		real(8) :: alpha
		real(8),intent(in) :: acsr(:),b(:),initGuess(:)
		real(8),allocatable :: r(:),p(:),temp(:)
		real(8),allocatable,intent(out) :: x(:)

		allocate(x(meshNumNodes))
		allocate(r(meshNumNodes))
		allocate(p(meshNumNodes))
		allocate(temp(meshNumNodes))
		x = initGuess
		call mkl_dcsrgemv("N",meshNumNodes,acsr,ia,ja,x,p)
		r = b-p
		call mkl_dcsrgemv("N",meshNumNodes,acsr,ia,ja,r,p)
		do i=1,femSolverMaxIter
			if(mod(i,1000).eq.0) then
				print *, "iteration: ", i
				print *, "residual ratio: ", norm2(r)/cc
			end if
			alpha = dot_product(p,r)/dot_product(p,p)
			x = x + alpha*r
			r = r-alpha*p
			if(norm2(r) .lt. cc) then
				iter = i
				return
			end if
			if(i.eq.femSolverMaxIter) then
				print *, "Maximum iterations reached."
				print *, "Convergence not achieved."
				print *, "Norm of residual: ", norm2(r)
				print *, "Convergence criterion: ", cc
				if((norm2(r)/cc) .lt. 2.d0) then
					write(*,*)"The residual is within a small range of &
					&the convergence criterion."
					write(*,*)"Increasing the iteration count may help."
				end if
			end if
			call mkl_dcsrgemv("N",meshNumNodes,acsr,ia,ja,r,p)
		end do
	end subroutine solveMinRes

	subroutine solveBiCGStab(acsr,ia,ja,b,initGuess,x,iter)
		integer :: i,j,k
		integer,intent(in) :: ia(:),ja(:)
		integer,intent(out) :: iter
		real(8),parameter :: cc=1e-9
		real(8) :: alpha,beta,delta0,delta,delta_old,omega
		real(8),intent(in) :: acsr(:),b(:),initGuess(:)
		real(8),allocatable :: r(:),p(:),s(:),rst(:),temp1(:),temp2(:)
		real(8),allocatable,intent(out) :: x(:)

		allocate(x(meshNumNodes))
		allocate(r(meshNumNodes))
		allocate(p(meshNumNodes))
		allocate(s(meshNumNodes))
		allocate(rst(meshNumNodes))
		allocate(temp1(meshNumNodes))
		allocate(temp2(meshNumNodes))

		x = initGuess

		call mkl_dcsrgemv("N",meshNumNodes,acsr,ia,ja,x,temp1)
		r = b-temp1
		call random_number(rst)
		p = r
		delta = dot_product(rst,r)
		write(*,'(a,1x,f15.3)') "Starting delta: ", delta
		delta0 = delta

		do i=1,femSolverMaxIter
			if(norm2(r) /= norm2(r)) then
				write(*,'(a)') "Error in solver: residual NaN"
				write(*,'(a)') "Check problem definition and mesh &
				&refinement, for error source."
				exit
			end if
			if(mod(i,1000).eq.0) then
				write(*,'(a,1x,i6)') 'Iteration number: ',i
				write(*,'(a,1x,f15.3)') "Residual ratio: ", norm2(r)/cc
			end if
			call mkl_dcsrgemv("N",meshNumNodes,acsr,ia,ja,p,temp1)	! temp1=A*p
			alpha = delta/dot_product(rst,temp1)
			s = r - alpha*temp1
			call mkl_dcsrgemv("N",meshNumNodes,acsr,ia,ja,s,temp2)	! temp2=A*s
			omega = dot_product(s,temp2)/dot_product(temp2,temp2)
			x = x + alpha*p + omega*s
			r = s - omega*temp2
			delta_old = delta
			delta = dot_product(rst,r)
			beta = (delta/delta_old)*(alpha/omega)
			p = r + beta*(p - omega*temp1)
			if(norm2(r) .lt. cc) then
				iter = i
				return
			end if
			if(i.eq.femSolverMaxIter) then
				write(*,'(a)') "Maximum iterations reached."
				write(*,'(a)') "Convergence not achieved."
				write(*,'(a,1x,f15.3)') "Norm of residual: ", norm2(r)
				write(*,'(a,1x,f15.3)') "Convergence criterion: ", cc
				if((norm2(r)/cc) .lt. 2.d0) then
					write(*,*)"The residual is within a small range of &
					&the convergence criterion."
					write(*,*)"Increasing the iteration count may help."
				end if
			end if
		end do
	end subroutine solveBiCGStab

	subroutine solveParDiSo(acsr,ia,ja,b,x)
		use mkl_pardiso
		integer,parameter :: dp = kind(1.0d0)
		type(mkl_pardiso_handle),allocatable  :: pt(:)
		integer :: i,maxfct,mnum,mtype,phase,n,nrhs,error,msglvl,nnz,	&
		error1
		integer,allocatable :: iparm(:)
		integer,intent(in) :: ia(:),ja(:)
		real(8) :: acsr(:),b(:)
		real(8),allocatable,intent(out) :: x(:)
		integer,dimension(1) :: idum
		real(8),dimension(1) :: ddum

		nnz = size(acsr,1)
		nrhs = 1
		maxfct = 1
		mnum = 1

		if(not(allocated(x))) allocate(x(meshNumNodes))
		allocate(iparm(64))		!set up pardiso control parameter

		do i=1,64
		   iparm(i) = 0
		end do
		iparm(1) = 1 ! no solver default
		iparm(2) = 2 ! fill-in reordering from metis
		iparm(4) = 0 ! no iterative-direct algorithm
		iparm(5) = 0 ! no user fill-in reducing permutation
		iparm(6) = 0 ! =0 solution on the first n compoments of x
		iparm(8) = 2 ! numbers of iterative refinement steps
		iparm(10) = 13 ! perturbe the pivot elements with 1e-13
		iparm(11) = 1 ! use nonsymmetric permutation and scaling mps
		iparm(13) = 0 ! maximum weighted matching algorithm is
					  !switched-off (default for symmetric).
					  ! try iparm(13) = 1 in case of inaccuracy
		iparm(14) = 0 ! output: number of perturbed pivots
		iparm(18) = -1 ! output: number of nonzeros in the factor lu
		iparm(19) = -1 ! output: mflops for lu factorization
		iparm(20) = 0 ! output: numbers of cg iterations

		error  = 0 ! initialize error flag
		msglvl = 0 ! 0=no output, 1=print statistical information
		mtype  = 11 ! real and unsymmetric matrix

		! Initiliaze the internal solver memory pointer.
		! This is only necessary for the first call of the solver.
		allocate (pt(64))
		do i=1,64
		   pt(i)%dummy =  0
		end do

		phase = 11 ! Only reordering and symbolic factorization
		call pardiso (pt,maxfct,mnum,mtype,phase,meshNumNodes,acsr,ia,	&
		ja,idum,nrhs,iparm,msglvl,ddum,ddum,error)
		if (error /= 0) then
		   write(*,*) 'the following error was detected: ', error
		   goto 1000
		end if

		phase = 22 ! only factorization
		call pardiso (pt,maxfct,mnum,mtype,phase,meshNumNodes,acsr,ia,	&
		ja,idum,nrhs,iparm,msglvl,ddum,ddum,error)
		if (error /= 0) then
		   write(*,*) 'the following error was detected: ', error
		   goto 1000
		endif
		! back substitution and iterative refinement
		iparm(8) = 2 ! max numbers of iterative refinement steps

		phase = 33 ! only solving
		call pardiso (pt,maxfct,mnum,mtype,phase,meshNumNodes,acsr,ia,	&
		ja,idum,nrhs,iparm,msglvl,b,x,error)
		write(*,*) 'solve completed ... '
		if (error /= 0) then
		   write(*,*) 'the following error was detected: ', error
		   goto 1000
		endif
1000 	continue
		! termination and release of memory

		phase = -1 ! release internal memory
		call pardiso (pt,maxfct,mnum,mtype,phase,meshNumNodes,ddum,		&
		idum,idum,idum,nrhs,iparm,msglvl,ddum,ddum,error1)
		if (error1 /= 0) then
		   write(*,*) 'the following release error was detected: ',	&
		   error1
		   stop 1
		endif
		if ( error /= 0 ) stop 1
	end subroutine solveParDiSo
!-----------------------------------------------------------------------!
!	End of solver routines
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
!	Routines for postprocessing of results
!-----------------------------------------------------------------------!

	subroutine writeVtkResults()
		integer,parameter :: resFid=246,nCorners=4,tetType=10
		integer :: i

		open(resFid,file=commResDir//trim(adjustl(femResFile))//commVTKExt)
		write(resFid,'(a)')"# vtk DataFile Version 1.0"
		write(resFid,'(a)')"3D Unstructured Grid of Linear Tetrahedrons"
		write(resFid,'(a)')"ASCII"
		write(resFid,'(a)')""
		write(resFid,'(a)')"DATASET UNSTRUCTURED_GRID"
		write(resFid,'(a,2x,i8,2x,a)')"POINTS ",meshNumNodes," double"
		do i=1,meshNumNodes
			write(resFid,*) meshVerts(i,:)
		end do
		write(resFid,'(a)')""
		write(resFid,'(a,2x,i8,2x,i8)')"CELLS ",meshNumElems,5*meshNumElems
		do i=1,meshNumElems
			write(resFid,'(i2,2x,3(i8,2x),i8)')nCorners,meshElems(i)%nodes-1
		end do
		write(resFid,'(a)')""
		write(resFid,'(a,2x,i8)')"CELL_TYPES ",meshNumElems
		do i=1,meshNumElems
			write(resFid,'(i4)')tetType
		end do
		if(meshNumDoms.gt.1) then
			write(resFid,'(a,2x,i8)')"CELL_DATA",meshNumElems
			write(resFid,'(a,2x,i2)')"FIELD FieldData", 1
			write(resFid,'(a,2x,i2,2x,i8,2x,a)')"Material",1,meshNumElems,"int"
			write(resFid,'(5(i4,2x))')meshElems%domain-1
		end if
		write(resFid,'(a)')""
		write(resFid,'(a,2x,i8)')"POINT_DATA ",meshNumNodes
		write(resFid,'(a)')"SCALARS temperature double"
		write(resFid,'(a)')"LOOKUP_TABLE default"
		do i=1,meshNumNodes
			write(resFid,*) femTvals(i)
		end do
		close(resFid)
	end subroutine writeVtkResults

	subroutine writeNodalResults()
		integer,parameter :: resFid=135
		integer :: i

		open(resFid,file=commResDir//trim(adjustl(femResFile))//commOutExt)
		do i=1,meshNumNodes
			write(resFid,'(4(e20.12,2x))') meshVerts(i,:),femTvals(i)
		end do
		close(resFid)
	end subroutine writeNodalResults

!-----------------------------------------------------------------------!
!	Routines for postprocessing of results
!-----------------------------------------------------------------------!

end module fem
