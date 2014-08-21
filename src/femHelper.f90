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

	subroutine runFem()
		integer :: pbDatFileNum

		if(femSysAssembled) then
			if(norm2(meshSources-sySrc) .lt. NANO) then
				write(*,*) "System already solved, no changes detected."
				stop
			end if
			call setupFinalEquations()
		else
			call femInitMesh()
			call assembleFemSystem()
			femSysAssembled = .true.
			call setupFinalEquations()
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
			call assembleNoderows(elNodes,stElem,elSt)
			sySt(elNodes) = stElem
			if(transient) then
				elC = systemCs(elDom)
				elRho = systemRhos(elDom)
				call getElementCapacitance(elC,elRho,elVol,elCp)
				cpElem = cpGlo(elNodes)
				call assembleNoderows(elNodes,cpElem,elCp)
				cpGlo(elNodes) = cpElem
			end if
		end do
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

end module femHelper
