! Mesh module: contains everything grid-related

module mesh

	use commonRoutines, only: indexedSort,checkIoError,invertReal4by4,&
						determinantReal4by4

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

	type elementBin
		integer,allocatable :: bin(:)
	end type elementBin

!	Module variables
	integer :: meshNumNodes,meshNumElems,meshNumDoms,meshNumSurfs
	integer,allocatable :: meshStCols(:),meshStRowPtr(:)
	real(8),allocatable :: meshTemperatures(:),meshForces(:),			&
	meshVerts(:,:),meshSt(:)
	type(tetraElement),allocatable :: meshElems(:)
	type(surface),allocatable :: meshSurfs(:)
!	Inputs
	integer :: nBins = 8
	character(*),parameter :: datDir="../data/",resDir="../results/",	&
	meshExt=".msh",dataExt=".dat"
	character(72) :: meshFile = "a"

	contains

!-----------------------------------------------------------------------!
!	The following subroutine and the subroutines it calls depend
!	on only using one OPEN and one CLOSE statement on the mesh
!	file and serially reading the mesh data into corresponding
!	module variables
!-----------------------------------------------------------------------!

	subroutine readMesh()
		integer,parameter :: fNum=100
		integer :: mDets(7)
		character(72) :: meshFileName

		meshFileName = datDir//trim(adjustl(meshFile))//meshExt
        call openmeshfile(fNum,meshFileName)
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

    subroutine openMeshFile(fNum,meshFileName)
		integer :: openStat
        integer,intent(in) :: fNum
		character(*),intent(in) :: meshFileName

        open(unit=fNum,file=meshFileName,form='formatted',status='old',	&
		action='read',iostat=openStat)
		call checkIoError(openStat,fNum)
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
		integer :: i,j,k,doSize,cSize,temp(10)
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
			allocate(meshSurfs(i)%elNum(numSfFaces))
			allocate(meshSurfs(i)%fcNum(numSfFaces))
			if(mod(numSfFaces,5) == 0) then
                cs = numSfFaces/5
            else
                cs = 1 + numSfFaces/5
            end if
			ct = 0
            do j=1,cs
                read(fNum,'(5(2x,i8,1x,i1))') temp
                do k=1,5
                    if(temp(2*k-1) .ne. 0) then
						ct = ct+1
                        meshSurfs(i)%elNum(ct) = temp(2*k-1)
						meshSurfs(i)%fcNum(ct) = temp(2*k)
						if(l) then
							meshElems(temp(2*k-1))%neighbours(temp(2*k),:) = 	&
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
		integer :: i,j
		type(elementBin),allocatable :: elBins(:,:,:)

		call binElements(elBins)
		call getBinAdjacency(elBins)
		do i=1,nBins
			call findNeighboursWithinBin(elBins(i))
			call findNeighboursAcrossBins(i,elBins)
		end do
	end subroutine getElementNeighbours

	subroutine binElements(elBins)
		integer :: i,j,nBPE,binNum,sz,cRs(3),elNodes(4)
		integer,allocatable :: tempBin(:)
		real(8) :: dmin(3),dmax(3),edges(3),elCent(3),elVerts(4,3)
		type(elementBin),allocatable,intent(out) :: elBins(:,:,:)

		nBPE = (nBins+1)**(1.0d0/3.0d0)
		allocate(elBins(nBPE,nBPE,nBPE))
		dmin = minval(meshVerts,1)
		dmax = maxval(meshVerts,1)
		edges = dmax-dmin
		do i=1,meshNumElems
			elNodes = meshElems(i)%nodes
			elVerts = meshVerts(elNodes,:)
			elCent = sum(elVerts,1)/4.0d0
			meshElems(i)%centroid = elCent
			cRs = ceiling(((elCent-dmin)/edges)*nBPE)
			where(cRs == 0)
				cRs = 1
			end where
			sz = size(elBins(cRs(1),cRs(2),cRs(3))%bin,1)
			if(sz.gt.0) then
				allocate(tempBin(sz+1))
				tempBin(1:sz) = elBins(cRs(1),cRs(2),cRs(3))%bin
				tempBin(sz+1) = i
				call move_alloc(tempBin,elBins(cRs(1),cRs(2),cRs(3))%bin)
			else
				allocate(elBins(cRs(1),cRs(2),cRs(3))%bin(1))
				elBins(cRs(1),cRs(2),cRs(3))%bin = i
			end if
		end do
	end subroutine binElements

	subroutine getBinAdjacency(elBins)
		integer,parameter :: nAdjBins=27,adjFilNum=200
		integer :: i,j,ct,ind,tempAdj(nAdjBins)
		character(*),parameter :: adjFile=datDir//"binadjacency.dat"
		type(elementBin),intent(inout) :: elBins(:)

		open(adjFilNum,file=adjFile)
		do i=1,nBins
			read(adjFilNum,*) tempAdj
			ct = count(tempAdj .ne. 0)
			allocate(elBins(i)%adjacency(ct))
			ind = 0
			do j=1,nAdjBins
 				if(tempAdj(j) .ne. 0) then
					ind = ind + 1
					elBins(i)%adjacency(ind) = tempAdj(j)
				end if
			end do
		end do
		close(adjFilNum)
	end subroutine getBinAdjacency

	subroutine findNeighboursWithinBin(singleBin)
		integer :: i,j,binSize,ct,el1,el2,fc1,fc2,ref,elNo1(4),elNo2(4)
		logical :: shared
		type(elementBin),intent(in) :: singleBin

		binSize = size(singleBin%bin,1)
		ref = 0
		do i=1,binSize-1
			el1 = singleBin%bin(i)
			ct = count(meshElems(el1)%neighbours(:,1) > 0)
			if(ct == 4) cycle
			elNo1 = meshElems(el1)%nodes
			shared = .false.
			do j=i+1,binSize
				el2 = singleBin%bin(j)
				elNo2 = meshElems(el2)%nodes
				call checkForSharedFace(elNo1,elNo2,shared,fc1,fc2)
				if(shared) then
					call addNeighbour(el1,el2,fc1,fc2)
					ref = ref+1
				end if
			end do
		end do
	end subroutine findNeighboursWithinBin

	subroutine findNeighboursAcrossBins(binNum,elBins)
		integer,intent(in) :: binNum
		integer :: i,j,k,ct,binSz1,binSz2,nAdj,adjBinNum,el1,el2,fc1,	&
		fc2,elNo1(4),elNo2(4)
		logical :: shared
		type(elementBin) :: givenBin,adjBin
		type(elementBin),intent(in) :: elBins(:)

		givenBin = elBins(binNum)
		binSz1 = size(givenBin%bin,1)
		nAdj = size(givenBin%adjacency,1)
		do i=1,binSz1
			el1 = givenBin%bin(i)
			ct = count(meshElems(el1)%neighbours(:,1) > 0)
			if(ct == 4) cycle
			elNo1 = meshElems(el1)%nodes
			shared = .false.
			do j=1,nAdj
				adjBinNum = elBins(binNum)%adjacency(j)
				if(adjBinNum == binNum) cycle
				adjBin = elBins(adjBinNum)
				binSz2 = size(adjBin%bin,1)
				do k=1,binSz2
					el2 = adjBin%bin(k)
					ct = count(meshElems(el2)%neighbours(:,1) > 0)
					if(ct == 4) cycle
					elNo2 = meshElems(el2)%nodes
					call checkForSharedFace(elNo1,elNo2,shared,fc1,fc2)
					if(shared) then
						call addNeighbour(el1,el2,fc1,fc2)
					end if
				end do
				if(shared) then
					exit
				end if
			end do
		end do
	end subroutine findNeighboursAcrossBins

	subroutine addNeighbour(el1,el2,fc1,fc2)
		integer,intent(in) :: el1,el2,fc1,fc2

		if(any((/el1,el2,fc1,fc2/) == 0)) then
			write(*,'(a)') "Wrong arguments passed to add neighbour."
			stop
		end if
		meshElems(el1)%neighbours(fc1,1) = el2
		meshElems(el1)%neighbours(fc1,2) = fc2
		meshElems(el2)%neighbours(fc2,1) = el1
		meshElems(el2)%neighbours(fc2,2) = fc1
	end subroutine addNeighbour

	subroutine checkForSharedFace(elNo1,elNo2,shared,fc1,fc2)
		integer :: p,q,sameCt,same1(3),same2(3)
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
			fc1 = getFaceIndex(same1)
			fc2 = getFaceIndex(same2)
		end if
	end subroutine checkForSharedFace

!-----------------------------------------------------------------------!
!	Element level functions for the mesh
!-----------------------------------------------------------------------!

	function getFaceIndex(fcNodes) result(fcNum)
		integer :: fcNum,indices(3)
		integer,intent(inout) :: fcNodes(3)

		call indexedSort(fcNodes,indices)
		if(all(fcNodes == (/1,2,3/))) then
			fcNum = 1
		elseif(all(fcNodes == (/1,2,4/))) then
			fcNum = 2
		elseif(all(fcNodes == (/2,3,4/))) then
			fcNum = 3
		elseif(all(fcNodes == (/1,3,4/))) then
			fcNum = 4
		else
			write(*,'(a)') "Face node indices not recognised."
			stop
		end if
	end function getFaceIndex

	function getFaceNodes(fcNum) result(fcNodes)
		integer :: fcNodes(3)
		integer,intent(in) :: fcNum

		if(fcNum == 1) then
			fcNodes = (/1,2,3/)
		elseif(fcNum == 2) then
			fcNodes = (/1,2,4/)
		elseif(fcNum == 3) then
			fcNodes = (/2,3,4/)
		elseif(fcNum == 4) then
			fcNodes = (/1,3,4/)
		else
			write(*,'(a)') "Face index not valid, must be between 1-4."
			stop
		end if
	end function getFaceNodes

	function getElementShapeFunctions(elNum) result(spFns)
		integer :: elNodes(4)
		integer,intent(in) :: elNum
		real(8) :: ec(4,3),vm(4,4),spFns(4,4)

		elNodes = meshElems(elNum)%nodes
		ec = meshVerts(elNodes,:)
		vm = getElementJacobian(ec)
		call invertReal4by4(vm)
		spFns = transpose(vm);
	end function getElementShapeFunctions

	function getElementJacobian(ec) result(elJac)
		real(8),intent(in) :: ec(4,3)
		real(8) :: elJac(4,4)

		elJac(:,2:4) = ec
		elJac(:,1) = 1.0d0
	end function getElementJacobian

	function getElementVolume(elNum) result(elVol)
		integer :: elNodes(4)
		integer,intent(in) :: elNum
		real(8) :: elVol,ec(4,3),elJac(4,4)

		elNodes = meshElems(elNum)%nodes
		ec = meshVerts(elNodes,:)
		elJac = getElementJacobian(ec)
		elVol = abs(determinantReal4by4(elJac))/6.0d0
	end function getElementVolume


	function getElementUnitStiffness(elNum) result(elemUnitSt)
		integer,intent(in) :: elNum
		real(8) :: ev,b(3,4),bt(4,3),spFns(4,4),elemUnitSt(4,4)

		spFns = getElementShapeFunctions(elNum)
		ev = getElementVolume(elNum)
		bt = spFns(:,2:4)
		b = transpose(bt)
		elemUnitSt = (matmul(bt,b))*ev
	end function getElementUnitStiffness

	function getElementCentroid(elNum) result(centroid)
		integer :: elNodes(4)
		integer,intent(in) :: elNum
		real(8) :: centroid(3),ec(4,3)

		elNodes = meshElems(elNum)%nodes
		ec = meshVerts(elNodes,:)
		centroid = sum(ec,1)/4.d0
	end function getElementCentroid

!-----------------------------------------------------------------------!
!	End of the element level functions for the mesh
!-----------------------------------------------------------------------!

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
		integer,intent(in) :: stCols(:),stRowPtr(:)
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
