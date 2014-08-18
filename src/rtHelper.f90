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

	subroutine getNextFace(ec,pt,dir,newFc,lToFc)
		integer :: fc,remNode,fcNodes(3)
		integer,intent(out) :: newFc
		real(8) :: ml,ptToFc,remVert(3),fcNorm(3),cent(3),fcVerts(3,3)
		real(8),intent(in) :: pt(3),dir(3),ec(4,3)
		real(8),intent(out) :: lToFc
		logical :: inward

		inward = .true.
		ml = 1000.0d0
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
		end do
	end subroutine createEmissionSurface

	function getRayPathLength() result(pathLength)
		real(8) :: randL,pathLength

		call random_number(randL)
		pathLength = (1.0d0/(kappa+sigma))*log(1.0d0/randL)
	end function getRayPathLength

end module rtHelper
