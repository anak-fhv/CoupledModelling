! Mesh module: contains everything grid-related

module mesh

	implicit none

!	Declare all relevant types
	type tetraElement
		integer :: domain,nodes(4),neighbours(4,2)
		real(8) :: volume,centroid(3)
	end type tetraElement

	type surface
		character(len=16) :: sfName
		integer :: sfId,numFcs
		integer,allocatable :: elNum(:),fcNum(:)
	end type surface

!	Module variables
	integer :: meshNumNodes,meshNumElems,meshNumDoms,meshNumSurfs
	real(8),allocatable :: meshTemperatures(:),meshForces(:),			&
	meshVerts(:,:)
	type(tetraElement),allocatable :: meshElems(:)
	type(surface),allocatable :: meshSurfs(:)
!	Inputs read from problem data file
	integer :: nBins = 8
	character(32) :: meshFile = "a.msh"

	contains

!-----------------------------------------------------------------------!
!	The following subroutine and the subroutines it calls depend
!	on only using one OPEN and one CLOSE statement on the mesh
!	file and serially reading the mesh data into corresponding
!	module variables
!-----------------------------------------------------------------------!

	subroutine readMesh()
		integer,parameter :: fNum=100
		integer :: openStat,mDets(7)

        call openmeshfile(fNum,meshFile)
        call readmeshdetails(fNum,mDets)
		meshNumNodes = mDets(1)
		meshNumElems = mDets(2)
		meshNumDoms = mDets(6)
		meshNumSurfs = mDets(7)
        call readmeshvertices(fNum)
        call readMeshElements(fNum)
        call readMeshElementDomains(fNum)
        call readMeshSurfaces(fNum)
        call closeMeshFile(fNum)
	end subroutine readMesh

    subroutine openMeshFile(fNum)
        integer,intent(in) :: fNum
        character(*),intent(in) :: meshFile

        open(unit=fNum,file=fname,form='formatted',status='old',		&
		action='read',iostat=openStat)
		call check_io_error(openStat)
    end subroutine openMeshFile

    subroutine readMeshDetails(fNum,mDets)
        integer,intent(in) :: fNum
		integer,intent(out) :: mDets(7)

        read(fNum, *)
        read(fNum, *)
        read(fNum, *) mDets
    end subroutine readMeshDetails

    subroutine readMeshVertices(fNum)
		integer :: i
        integer,intent(in) :: fNum

		if(.not.(allocated(meshVerts))) then
	        allocate(meshVerts(meshNumNodes,3))
		end if
        do i=1,meshNumNodes
            read(fNum, *) meshVerts(i,:)
        end do
    end subroutine readMeshVertices

    subroutine readMeshElements(fNum)
		integer :: i
        integer,intent(in) :: fNum

		if(.not.(allocated(meshElems))) then
			allocate(meshElems(meshNumElems))
		end if
		do i=1,meshNumElems
            read(fNum, *) meshElems(i)%nodes
        end do
    end subroutine readMeshElements

	subroutine readMeshElementDomains(fNum)
		integer :: i,j,doSize,cSize,temp(10)
        integer,intent(in) :: fNum

        do i=1,meshNumDoms
            read(fNum,'(i8)') doSize
            if(mod(doSize,10) == 0) then
                cSize = doSize/10
            else
                cSize = 1 + (doSize/10)
            end if
            do j=1,cSize
                read(fNum,'(10(1x,i8))') temp
                do k=1,10
                    if(temp(k) .ne. 0) then
                        meshElems(temp(k))%domain = i
                    end if
                end do
            end do
        end do
	end subroutine readMeshElementDomains

	subroutine readMeshSurfaces(fNum)
		integer :: i,j,k,cs,ct,numSfFaces,temp(10)
        integer,intent(in) :: fNum
		character(16) :: surfName
		logical :: l
		
		if(.not.(allocated(meshSurfs))) then
	        allocate(meshSurfs(meshNumSurfs))
		end if
		do i=1,meshNumSurfs
			meshSurfs%sfId = i
            read(fNum,'(i8,1x,a)') numSfFaces,surfName
			meshSurfs(i)%numFcs = numSfFaces
			meshSurfs(i)%sfName = surfName
			l = (verify('iface',surfName) .ne. 0)
			l = l .and. (verify('bic',surfName) .ne. 0)
			allocate(meshSurfs(i)%elNum(numFc))
			allocate(meshSurfs(i)%fcNum(numFc))
			if(mod(numFc,5) == 0) then
                cs = numFc/5
            else
                cs = 1 + numFc/5
            end if
			ct = 0
            do j=1,cs
                read(fno,'(5(2x,i8,1x,i1))') temp
                do k=1,5
                    if(temp(2*k-1) .ne. 0) then
						ct = ct+1
                        meshSurfs(i)%elNum(ct) = temp(2*k-1)
						meshSurfs(i)%fcNum(ct) = temp(2*k)
						if(l) then
							meshElems(temp(2*k-1))%neighbours(2*k,:) = 	&
							(/temp(2*k-1),-i/)
						end if
                    end if
                end do
            end do
		end do
	end subroutine readMeshSurfaces

	subroutine closeMeshFile(fNum)
		integer,intent(in) :: fNum

		close(fNum)
	end subroutine closeMeshFile

!-----------------------------------------------------------------------!
!	End of the serial mesh file reading routines
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
!	The subroutines starting here are all auxiliary to the main mesh
!	module and may be called only if required by the implementation,
!	for example, the routines to populate the neighbours of tetrahedra
!-----------------------------------------------------------------------!

	subroutine binElements()
		integer :: i,j,bPerEdge,centRatios(3),elNodes(4)
		real(8) :: dmin(3),dmax(3),edges(3),elCent(3),elVerts(4,3)

		bPerEdge = (nBins+1)**(1.0d0/3.0d0)
		dmin = minval(meshVerts,1)
		dmax = maxval(meshVerts,1)
		edges = dmax-dmin
		do i=1,meshNumElems
			elNodes = meshElems(i)%nodes
			elVerts = meshVerts(elNodes,:)
			elCent = sum(elVerts,1)/3.0d0
			meshElems(i)%centroid = elCent
			centRatios = ceiling(((elCent-dmin)/edges)*bPerEdge)
			where(centRatios == 0)
				centRatios = 1
			end where
		end do
	end subroutine binElements

	subroutine getElementNeighbours()
	end subroutine getElementNeighbours

!-----------------------------------------------------------------------!
!	End of the serial mesh file reading routines
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
!	The subroutines that follow are to be used to set values of the 
!	mesh properties such as temperature and FEM source (forcing) term
!-----------------------------------------------------------------------!

	subroutine setMeshNodalValues(nodalVals,propNameFlag)
		real(8),intent(in) :: nodalVals(:)
		character(*),intent(in) :: propNameFlag

		if(size(nodalVals,1) .ne. meshNumNodes) then
			write(*,'(a)') "The temperature field returned does not &
			& match the mesh. Please verify your allocations."
			stop
		end if

		if(propNameFlag .eq. "T") then
			if(.not.(allocated(meshTemperatures))) then
				allocate(meshTemperatures(meshNumNodes))
			end if
			meshTemperatures = nodalVals
		elseif(propNameFlag .eq. "F") then
			if(.not.(allocated(meshForces))) then
				allocate(meshForces(meshNumNodes))
			end if
			meshForces = nodalVals
		else
			write(*,'(a)') "Property to be set not recognised."
			stop
		end if

	end subroutine setMeshNodalValues

!-----------------------------------------------------------------------!
!	End of the setter subroutines
!-----------------------------------------------------------------------!

end module mesh
