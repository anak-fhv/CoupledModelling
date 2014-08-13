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

	type elementBins
		integer,allocatable :: bin(:),adjacency(:)
	end type elementBins

!	Module variables
	integer :: meshNumNodes,meshNumElems,meshNumDoms,meshNumSurfs
	integer,allocatable :: meshStCols(:),meshStRowPtr(:)
	real(8),allocatable :: meshTemperatures(:),meshForces(:),			&
	meshVerts(:,:),meshSt(:)
	type(tetraElement),allocatable :: meshElems(:)
	type(surface),allocatable :: meshSurfs(:)
!	Inputs read from problem data file
	integer :: nBins = 8
	character(72) :: meshFile = "a.msh"

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

        open(unit=fNum,file=meshFile,form='formatted',status='old',		&
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

	subroutine getElementNeighbours()
		integer,parameter :: adjFilNum=200
		integer :: i,j
		integer,allocatable :: binAdjacency(:,:)
		type(elementBins),allocatable :: elBins(:)

		call binElements(elBins)
		call getBinAdjacency(elBins)
		do i=1,nBins
			call findNeighboursWithinBin()
			call findNeighboursAcrossBins()
		end do
	end subroutine getElementNeighbours

	subroutine binElements(elBins)
		integer :: i,j,nBPE,binNum,sz,cRs(3),elNodes(4)
		integer,allocatable :: tempBin(:)
		real(8) :: dmin(3),dmax(3),edges(3),elCent(3),elVerts(4,3)
		type(elementBins),allocatable,intent(out) :: elBins(:)

		allocate(elBins(nBins))
		nBPE = (nBins+1)**(1.0d0/3.0d0)
		dmin = minval(meshVerts,1)
		dmax = maxval(meshVerts,1)
		edges = dmax-dmin
		do i=1,meshNumElems
			elNodes = meshElems(i)%nodes
			elVerts = meshVerts(elNodes,:)
			elCent = sum(elVerts,1)/3.0d0
			meshElems(i)%centroid = elCent
			cRs = ceiling(((elCent-dmin)/edges)*nBPE)
			where(cRs == 0)
				cRs = 1
			end where
			binNum = cRs(1) + (cRs(2)-1)*nBPE + (cRs(3)-1)*nBPE^2
			sz = size(elBins(binNum)%bin,1)
			if(sz.gt.0) then
				allocate(tempBin(sz+1))
				tempBin(1:sz) = elBins(binNum)%bin
				tempBin(sz+1) = i
				call move_alloc(tempBin,elbins(binNum)%bin)
			else
				allocate(elBins(binNum)%bin(1))
				elBins(binNum)%bin = i
			end if
		end do
	end subroutine binElements

	subroutine getBinAdjacency(elBins)
		type(elementBins),allocatable,intent(inout) :: elBins(:)
	end subroutine getBinAdjacency

	subroutine findNeighboursWithinBin()
		integer :: i,j,p,q,k,binSize,el1,el2,fc1,fc2,elNo1(4),elNo2(4)
		logical :: shared
		type(elementBin) :: singleBin

		binSize = size(singleBin%bin,1)
		do i=1,binSize-1
			el1 = singleBin%bin(i)
			ct = count(meshElems(el1)%neighbours(:,1) > 0)
			if(ct == 4) cycle
			elNo1 = meshElems(el1)%nodes
			do j=i+1,binSize
				el2 = singleBin%bin(j)
				elNo2 = meshElems(el2)%nodes
				call checkForSharedFace(elNo1,elNo2,shared,fc1,fc2)
				if(shared) then
					call addNeighbour(el1,el2,fc1,fc2)
				end if
			end do
		end do
	end subroutine findNeighboursWithinBin

	subroutine findNeighboursAcrossBins()
		
	end subroutine findNeighboursAcrossBins

	subroutine addNeighbour(el1,el2,fc1,fc2)
		integer,intent(in) :: el1,el2,fc1,fc2

		if((any(/el1,el2,fc1,fc2/) == 0)) then
			write(*,'(a)') "Wrong arguments passed to add neighbour."
			stop
		end if
		meshElems(el1)%neighbours(fc1,1) = el2
		meshElems(el1)%neighbours(fc1,2) = fc2
		meshElems(el2)%neighbours(fc2,1) = el1
		meshElems(el2)%neighbours(fc2,2) = fc1
	end subroutine addNeighbour

	subroutine checkForSharedFace(elNo1,elNo2,shared,fc1,fc2)
		integer :: p,q
		integer,intent(in) :: elNo1(4),elNo2(4)
		integer,intent(out) :: fc1,fc2
		logical,intent(out) :: shared

		same1 = 0
		same2 = 0
		sameCt = 0
		shared = .false.
		fc1 = 0
		fc2 = 0
		do p=1,4
			do q=1,4
				if(elNo1(p) == elNo2(q)) then
					sameCt = sameCt+1
					same1(sameCt) = p
					same2(sameCt) = q
				end if
			end do
		end do
		if(sameCt == 3) then
			shared = .true.
			fc1 = getFaceFromNodes(same1)
			fc2 = getFaceFromNodes(same2)
		end if
	end subroutine compareNodes

	function getFaceFromNodes(fcNodes) result(fcNum)
		integer :: indices(3)
		integer,intent(in) :: fcNodes(3)
		integer,intent(out) :: fcNum

		call indexedSort(fcNodes,indices)
		if(fcNodes == (/1,2,3/)) then
			fcNum = 1
		elseif(fcNodes == (/1,2,4/)) then
			fcNum = 2
		elseif(fcnodes == (/2,3,4/)) then
			fcNum = 3
		elseif(fcnodes == (/1,3,4/)) then
			fcNum = 4
		else
			write(*,'(a)') "Face node indices not recognised."
			stop
		end if
	end function getFaceFromNodes

	subroutine indexedSort(a,b)
		integer,intent(in) :: a
		integer,intent(out) :: b(size(a,1))
		integer :: i,j,n,temp

		n = size(a,1)
		b = (/1:n/)
		do i=1,n-1
			do j=i+1,n
				if(a(i) > a(j)) then
					temp = a(j)
					a(j) = a(i)
					a(i) = temp
					temp = b(i)
					b(i) = b(j)
					b(j) = temp
				end if
			end do
		end do
	end subroutine indexedsort
!-----------------------------------------------------------------------!
!	End of the auxiliary routines
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

	subroutine saveMeshNodalStiffness(noSt,stCols,stRowPtr)
		integer,intent(in) :: stCols(:),stRowPtr
		real(8),intent(in) :: noSt(:)

		if(.not.(allocated(meshSt))) then
			allocate(meshSt(size(noSt,1)))
		end if
		if(.not.(allocated(meshStCols))) then
			allocate(meshStCols(size(stCols,1)))
		end if
		if(.not.(allocated(meshStRowPtr))) then
			allocate(meshStRowPtr(size(stRowPtr,1)))
		end if

		meshSt = noSt
		meshStCols = stCols
		meshStRowPtr = stRowPtr
	end subroutine saveMeshNodalStiffness

!-----------------------------------------------------------------------!
!	End of the setter subroutines
!-----------------------------------------------------------------------!

end module mesh
