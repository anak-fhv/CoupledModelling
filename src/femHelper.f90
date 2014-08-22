include 'mkl_pardiso.f90'	! Added for ParDiSo

module femHelper

	use commonRoutines, only: pi,sigb,LARGE,MEGA,PICO,commResDir,		&
						commOutExt,selTetraPoint,selTriPoint,			&
						getRaySphDir,getFaceRayDir,getFaceNorm,			&
						getFaceNodes,insideFaceCheck
	use mesh, only: meshFile,meshNBins,meshNumNodes,meshNumElems,		&
			  meshNumDoms,meshNumSurfs,meshElems,meshVerts,meshSurfs,	&
			  meshTemperatures,readMesh,getElementNeighbours,			&
			  populateSurfaceFaceAreas,setSurfaceConstTemperature

	implicit none

!	Declare all relevant types
	type nodeRow
		integer,allocatable :: col(:)
		real(8),allocatable :: val(:)
	end type nodeRow

	contains

!	Module variables
	integer :: rtNumRays
	integer,allocatable :: rtElemAbs(:),rtElemSto(:),rtEmSfIds(:),		&
	rtWallInf(:)
	real(8) :: rtKappa,rtSigma
	type(nodeRow),allocatable :: femSysSt(:),femSysCp(:)

	subroutine runFem()
		integer :: pbDatFileNum

		if(femSysAssembled) then
			if(norm2(meshSources-sySrc) .lt. NANO) then
				write(*,*)"System already solved, no changes detected."
				stop
			end if
			call setupFinalEquations()
			call solveFemSystem()
			call writeResults()
		else
			call femInitMesh()
			call assembleFemSystem()
			femSysAssembled = .true.
			call solveFemSystem()
			call writeResults()
		end if
	end subroutine runFem

	subroutine femInitMesh()
		character(*),intent(in) :: mFileName

		meshFile = trim(adjustl(mFileName))
		call readMesh()
	end subroutine femInitMesh

	subroutine assembleFemSystem(transient)
		integer :: i,j,k,elNodes(4)
		logical,intent(in) :: transient
		type(nodeRow) :: stElem(4),cpElem(4)
		type(nodeRow),allocatable :: sySt(:),cpGlo(:)

		allocate(sySt(meshNumNodes))
		if(transient) then
			allocate(cpGlo(meshNumNodes))
		end if
		do i=1,meshNumElems
			elNodes = meshElems(i)%nodes
			elDom = meshElems(i)%domain
			elK = systemKs(:,:,elDom)
			call getElementStiffness(i,elK,elVol,elSt)
			elBys = meshElems(i)%neighbours(:,2)
			if(all(elBys>0)) then
				continue
			else
				call boundaryConditions(i,elBys,elBySt,elByTs,elBySrc)
				elSt = elSt + elBySt
				sySrc(elNodes) = elBySrc
				syTvals(elNodes) = elByTs
			end if
			stElem = sySt(elNodes)
			call assembleNoderows(elNodes,elSt,stElem)
			sySt(elNodes) = stElem
			if(transient) then
				elC = systemCs(elDom)
				elRho = systemRhos(elDom)
				call getElementCapacitance(elC,elRho,elVol,elCp)
				cpElem = cpGlo(elNodes)
				call assembleNoderows(elNodes,elCp,cpElem)
				cpGlo(elNodes) = cpElem
			end if
		end do
		call setupFinalEquations()
		call collapseNodeRows()
	end subroutine assembleFemSystem

	subroutine setupFinalEquations(sySt,sySrc,syTvals)
		integer :: i
		real(8) :: sySrc(:),syTvals(:)
		type(nodeRow),intent(inout):: sySt(:)

		do i=1,meshNumNodes
			if(syTvals(i) .ne. 0.0d0) then
				if(allocated(stNo(i)%col)) then
					deallocate(stNo(i)%val)
					deallocate(stNo(i)%col)
				end if
				allocate(stNo(i)%col(1))
				allocate(stNo(i)%val(1))
				stNo(i)%col(1) = i
				stNo(i)%val(1) = 1.0d0
				sySrc(i) = syTvals(i)
			end if
		end do
	end subroutine setupFinalEquations

!-----------------------------------------------------------------------!
!	Routines for boundary handling
!-----------------------------------------------------------------------!
	subroutine boundaryConditions(elNum,elBys,elBySt,elByTs,elBySrc)
		integer :: i,j,bloc,fbtype,n1,n2,fcNodes(3),boSfs(4),elNodes(4) 
		integer,intent(in) :: elNum,elBys(4)
		real(8) :: ec(4,3),elQ(4),htA(4)
		real(8),intent(out) :: elByTs(4),elBySrc(4),elBySt(4,4)
		real(8),dimension(:),intent(in) :: bvals
		integer,dimension(4) :: bofcs
		integer,dimension(3) :: fcnodes

		bstiff = 0.0d0
		bforce = 0.0d0
		btemp = 0.0d0
		elq = 0.0d0
		hta = 0.0d0
		elht = 0.0d0
		do i=1,4
			if(bofcs(i) /= 0) then
				bloc = bofcs(i)
				fbtype = bcs(bloc)
				bv = bvals(bloc)
				call bfacenodes(i,fcnodes)
				if(bv.ne.0.0d0) then
					if(fbtype == 1) then
						btemp(fcnodes) = bv
					elseif(fbtype == 2) then
						call fluxboundary(ec,fcnodes,bv,elq)
						bforce = bforce - elq
					elseif(fbtype == 3) then
						call convectiveboundary(ec,fcnodes,bv,Tamb,hta,elht)
						bforce = bforce + hta
						bstiff = bstiff + elht
					else
						continue
					end if
				end if
			else
				continue
			end if
		end do
	end subroutine boundaryConditions

	subroutine fluxboundary(ec,fcnodes,bv,elq)
		real(8),dimension(4,3) :: ec
		real(8),dimension(3,3) :: fc
		real(8),dimension(4) :: elq
		real(8) :: fcarea,bv
		integer,dimension(3) :: fcnodes
		integer :: i,j,k,el

		elq = 0.0d0
		fc = ec(fcnodes,:)
		fcarea = facearea(fc)
		elq(fcnodes) = (1.0d0/6.0d0)*(2.0d0*fcarea)*bv
	end subroutine fluxboundary

	subroutine convectiveboundary(ec,fcnodes,bv,Tamb,hta,elht)
		real(8),dimension(4,3) :: ec
		real(8),dimension(4,4) :: sp,elht,surfint
		real(8),dimension(3,3) :: fc
		real(8),dimension(4) :: hta
		real(8) :: Tamb,fcarea,bv
		integer,dimension(3) :: fcnodes
		integer :: fnum,i,bl

		hta = 0.0d0
		elht= 0.0d0
		fc = ec(fcnodes,:)
		do i=1,3
		end do
		fcarea = facearea(fc)
		hta(fcnodes) = (1.0d0/6.0d0)*(2.0d0*fcarea)*bv*Tamb
		surfint = shapefuncsquaresurfint(fcnodes)
		elht = (2.0d0*fcarea/24.0d0)*bv*surfint
	end subroutine convectiveboundary
!-----------------------------------------------------------------------!
!	End of boundary handling routines
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
!	Routines required to assemble the sparse global stiffness matrix
!-----------------------------------------------------------------------!

	subroutine assembleNoderows(elNodes,elSt,nodeRowElem)
		integer :: i,sz
		integer,intent(in) :: elnodes(4)
		integer,allocatable :: tempcol(:)
		real(8) :: elSt(4,4)
		real(8),allocatable :: tempval(:)
		type(nodeRow),intent(inout) :: nodeRowElem(4)

		do i=1,4
			if(.not.(allocated(nodeRowElem(i)%col))) then
				allocate(nodeRowElem(i)%col(4))
				allocate(nodeRowElem(i)%val(4))
				nodeRowElem(i)%col = elnodes
				nodeRowElem(i)%val = btdb(i,:)
			else
				sz = size(nodeRowElem(i)%col,1)
				allocate(tempcol(sz+4))
				allocate(tempval(sz+4))
				tempcol = (/nodeRowElem(i)%col,elnodes/)
				tempval = (/nodeRowElem(i)%val,btdb(i,:)/)
				call move_alloc(tempcol,nodeRowElem(i)%col)
				call move_alloc(tempval,nodeRowElem(i)%val)
				if(allocated(tempcol)) deallocate(tempcol)
				if(allocated(tempval)) deallocate(tempval)
			end if
		end do

	end subroutine assemble_noderows

	subroutine collapseNodeRows(node,val,col,row_ptr)
		integer :: i,numUnique,numNo,numVal,sz
		integer,allocatable :: tempcol(:),ind(:),
		integer,allocatable,intent(out) :: sysRowPtr(:),sysColPtr(:)
		real(8),allocatable :: tempval(:),sysVals(:)
		type(nodeRow),intent(in) :: sysNodeRows(:)

		allocate(sysRowPtr(meshNumNodes+1))
		sysRowPtr(1) = 1
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
			sysRowPtr(i+1) = sysRowPtr(i)+numUnique
		end do

		allocate(sysVals(numVal))
		allocate(sysColPtr(numVal))
		do i=1,meshNumNodes
			indStart = sysRowPtr(i)
			indEnd = (sysRowPtr(i+1)-1)
			sysVals(indStart:indEnd) = sysNodeRows(i)%val
			sysColPtr(indStart:indEnd) = sysNodeRows(i)%col
			deallocate(sysNodeRows(i)%col)
			deallocate(sysNodeRows(i)%val)
		end do

	end subroutine collapseNodeRows

	subroutine addNodeEntries(nodeCols,nodeVals,numUnique)
		integer :: i,j,k
		integer,intent(out) :: numUnique
		integer,allocatable :: tempCols(:),starts(:)
		integer,allocatable,intent(inout) :: nodeCols(:)
		real(8),allocatable :: tempVals(:)
		real,allocatable,intent(inout) :: nodeVals(:)

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
	end subroutine add_duplicates

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
!	Different solver routines
!-----------------------------------------------------------------------!

	subroutine solveFemSystem()
		if(femSolverType == "BiCGStab") then
			call solveBiCGStab(acsr,ia,ja,b,maxNumIter,initGuess,x,iter)
		elseif(femSolverType == "MinRes") then
			call solveMinRes(acsr,ia,ja,b,maxNumIter,initGuess,x,iter)
		elseif(femSolverType == "ParDiSo") then
			call solveParDiSo(acsr,ia,ja,b,x)
		else
			write(*,*)"Solver type not recognised."
			write(*,*)"Using default BiCGStab solver."
			call solveBiCGStab(acsr,ia,ja,b,maxiter,initGuess,x,iter)
		end if
	end subroutine solveFemSystem

	subroutine solveMinRes(acsr,ia,ja,b,maxNumIter,initGuess,x,iter)
		integer :: i,j,k
		integer,intent(in) :: maxNumIter,ia(:),ja(:)
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
		do i=1,maxNumIter
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
			if(i.eq.maxiter) then
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

	subroutine solveBiCGStab(acsr,ia,ja,b,maxNumIter,initGuess,x,iter)
		integer :: i,j,k
		integer,intent(in) :: maxNumIter,ia(:),ja(:)
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

		do i=1,maxiter
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
			if(i.eq.maxiter) then
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
		real(8),intent(in) :: acsr(:),b(:)
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
		call pardiso(pt,maxfct,mnum,mtype,phase,meshNumNodes,acsr,ia,	&
		ja,idum,nrhs,iparm,msglvl,ddum,ddum,error)
		if (error /= 0) then
		   write(*,*) 'the following error was detected: ', error
		   goto 1000
		end if

		phase = 22 ! only factorization
		call pardiso(pt,maxfct,mnum,mtype,phase,meshNumNodes,acsr,ia,	&
		ja,idum,nrhs,iparm,msglvl,ddum,ddum,error)
		if (error /= 0) then
		   write(*,*) 'the following error was detected: ', error
		   goto 1000
		endif
		! back substitution and iterative refinement
		iparm(8) = 2 ! max numbers of iterative refinement steps

		phase = 33 ! only solving
		call pardiso(pt,maxfct,mnum,mtype,phase,meshNumNodes,Kcsr,ia,	&
		ja,idum,nrhs,iparm,msglvl,b,x,error)
		write(*,*) 'solve completed ... '
		if (error /= 0) then
		   write(*,*) 'the following error was detected: ', error
		   goto 1000
		endif
1000 	continue
		! termination and release of memory

		phase = -1 ! release internal memory
		call pardiso(pt,maxfct,mnum,mtype,phase,meshNumNodes,ddum,		&
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

end module femHelper
