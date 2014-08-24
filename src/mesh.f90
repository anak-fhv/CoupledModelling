! Mesh module: contains everything grid-related

module mesh

	use commonRoutines, only: commDatDir,commResDir,commMeshExt,		&
							  commDatExt,indexedSortInteger,			&
							  checkIoError,getFaceIndex,getFaceNodes,	&
							  invertReal4by4,determinantReal4by4,		&
							  triangleArea

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
		real(8) :: totalArea
		real(8),allocatable :: fcArea(:)
	end type surface

	type elementBin
		integer,allocatable :: bin(:)
	end type elementBin

!	Module variables
	integer :: meshNumNodes,meshNumElems,meshNumDoms,meshNumSurfs
	integer,allocatable :: meshStColPtr(:),meshStRowPtr(:)
	real(8),allocatable :: meshTemperatures(:),meshSources(:),			&
	meshVerts(:,:),meshStVals(:)
	type(tetraElement),allocatable :: meshElems(:)
	type(surface),allocatable :: meshSurfs(:)
!	Inputs
	integer :: meshNBins(3) = (/2,2,2/)
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

		meshFileName = commDatDir//trim(adjustl(meshFile))//commMeshExt
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
		integer :: i,j,k,elNum,fcNum,cs,ct,numSfFaces,temp(10)
        integer,intent(in) :: fNum
		character(16) :: surfName
		logical :: l
		
		if(.not.(allocated(meshSurfs))) then
	        allocate(meshSurfs(meshNumSurfs))
		end if
		do i=1,meshNumSurfs
			meshSurfs(i)%sfId = i
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
						elNum = temp(2*k-1)
						fcNum = temp(2*k)
                        meshSurfs(i)%elNum(ct) = elNum
						meshSurfs(i)%fcNum(ct) = fcNum
						if(l) then
							meshElems(elNum)%neighbours(fcNum,:) = 		&
							(/elNum,-meshSurfs(i)%sfId/)
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
		integer,parameter :: datFileNum=101
		integer :: i,j,k,temp(8)
		logical :: datFileExist
		character(72) :: fDat
		type(elementBin),allocatable :: elBins(:,:,:)

		fDat = commDatDir//trim(adjustl(meshFile))//commDatExt
		inquire(file=trim(adjustl(fDat)),exist=datFileExist)
		if(datFileExist) then
			do i=1,meshNumElems
				read(datFileNum,'(1x,4(i8,1x,i4))') temp
				do j=1,4
					meshElems(i)%neighbours(j,:) = temp(2*j-1:2*j)
				end do
			end do
			return
		end if
		call binElements(elBins)
		do k=1,meshNBins(3)
			do j=1,meshNBins(2)
				do i=1,meshNBins(1)
					call findNeighboursWithinBin(elBins(i,j,k))
					call findNeighboursAcrossBins(i,j,k,elBins)
				end do
			end do
		end do
	end subroutine getElementNeighbours

	subroutine binElements(elBins)
		integer :: i,j,sz,cRs(3),elNodes(4)
		integer,allocatable :: temp(:)
		real(8) :: dmin(3),dmax(3),edges(3),elCent(3),elVerts(4,3)
		type(elementBin),allocatable,intent(out) :: elBins(:,:,:)

		allocate(elBins(meshNBins(1),meshNBins(2),meshNBins(3)))
		dmin = minval(meshVerts,1)
		dmax = maxval(meshVerts,1)
		edges = dmax-dmin
		do i=1,meshNumElems
			elNodes = meshElems(i)%nodes
			elVerts = meshVerts(elNodes,:)
			elCent = sum(elVerts,1)/4.0d0
			meshElems(i)%centroid = elCent
			cRs = ceiling(((elCent-dmin)/edges)*meshNBins)
			where(cRs == 0)
				cRs = 1
			end where
			sz = size(elBins(cRs(1),cRs(2),cRs(3))%bin,1)
			if(sz.gt.0) then
				allocate(temp(sz+1))
				temp(1:sz) = elBins(cRs(1),cRs(2),cRs(3))%bin
				temp(sz+1) = i
				call move_alloc(temp,elBins(cRs(1),cRs(2),cRs(3))%bin)
			else
				allocate(elBins(cRs(1),cRs(2),cRs(3))%bin(1))
				elBins(cRs(1),cRs(2),cRs(3))%bin = i
			end if
		end do
	end subroutine binElements

	subroutine findNeighboursWithinBin(singleBin)
		integer :: i,j,binSize,ct1,ct2,el1,el2,fc1,fc2,ref,elNo1(4),	&
		elNo2(4)
		logical :: shared
		type(elementBin),intent(in) :: singleBin

		binSize = size(singleBin%bin,1)
		ref = 0
		do i=1,binSize-1
			el1 = singleBin%bin(i)
			ct1 = count(meshElems(el1)%neighbours(:,1)>0)
			if(ct1 == 4) cycle
			elNo1 = meshElems(el1)%nodes
			do j=i+1,binSize
				shared = .false.
				el2 = singleBin%bin(j)
				elNo2 = meshElems(el2)%nodes
				ct2 = count(meshElems(el2)%neighbours(:,1)>0)
				if(ct2 == 4) cycle
				call checkForSharedFace(elNo1,elNo2,shared,fc1,fc2)
				if(shared) then
					call addNeighbour(el1,el2,fc1,fc2)
					ct1 = count(meshElems(el1)%neighbours(:,1)>0)
					ref = ref+1
				end if
				if(ct1 == 4) exit
			end do
		end do
	end subroutine findNeighboursWithinBin

	subroutine findNeighboursAcrossBins(c1,c2,c3,elBins)
		integer,intent(in) :: c1,c2,c3
		integer :: i,j,k,ol,il,ct1,ct2,binSz1,binSz2,el1,el2,fc1,fc2,	&
		elNo1(4),elNo2(4)
		logical :: shared
		type(elementBin) :: givenBin,adjBin
		type(elementBin),intent(in) :: elBins(:,:,:)

		givenBin = elBins(c1,c2,c3)
		binSz1 = size(givenBin%bin,1)
		currentBin: do ol=1,binSz1
			el1 = givenBin%bin(ol)
			ct1 = count(meshElems(el1)%neighbours(:,1)>0)
			if(ct1 == 4) cycle
			elNo1 = meshElems(el1)%nodes
			zBinIndex: do k=c3-1,c3+1
				if((k==0).or.(k.gt.meshNBins(3))) cycle
				yBinIndex: do j=c2-1,c2+1
					if((j==0).or.(j.gt.meshNBins(2))) cycle
					xBinIndex: do i=c1-1,c1+1
						if((i==0).or.(i.gt.meshNBins(1))) cycle
						if(all((/i,j,k/)==(/c1,c2,c3/))) cycle
						adjBin = elBins(i,j,k)
						binSz2 = size(adjBin%bin,1)
						adjacentBin: do il=1,binSz2
							shared = .false.
							el2 = adjBin%bin(il)
							ct2 = count(meshElems(el2)%neighbours(:,1)>0)
							if(ct2 == 4) cycle
							elNo2 = meshElems(el2)%nodes
							call checkForSharedFace(elNo1,elNo2,shared,	&
							fc1,fc2)
							if(shared) then
								call addNeighbour(el1,el2,fc1,fc2)
								ct1 = count(meshElems(el1)%neighbours(:,1)>0)
							end if
							if(ct1 == 4) exit
						end do adjacentBin
						if(ct1 == 4) exit
					end do xBinIndex
					if(ct1 == 4) exit
				end do yBinIndex
				if(ct1 == 4) exit
			end do zBinIndex
		end do currentBin
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

	subroutine populateElementVolumes()
		integer :: i

		do i=1,meshNumElems
			meshElems(i)%volume = getElementVolume(i)
		end do
	end subroutine populateElementVolumes

	subroutine populateSurfaceFaceAreas()
		integer :: i,j,nSfFaces,elNum,fcNum
		real(8) :: aTot

		do i=1,meshNumSurfs
			nSfFaces = meshSurfs(i)%numFcs
			if(.not.(allocated(meshSurfs(i)%fcArea))) then
				allocate(meshSurfs(i)%fcArea(nSfFaces))
			end if
			aTot = 0.0d0
			do j=1,nSfFaces
				elNum = meshSurfs(i)%elNum(j)
				fcNum = meshSurfs(i)%fcNum(j)
				meshSurfs(i)%fcArea(j) = getFaceArea(elNum,fcNum)
				aTot = aTot + meshSurfs(i)%fcArea(j)
			end do
			meshSurfs(i)%totalArea = aTot
		end do
	end subroutine populateSurfaceFaceAreas
!-----------------------------------------------------------------------!
!	End of the auxiliary routines
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
!	Element level functions for the mesh
!-----------------------------------------------------------------------!

	function getFaceArea(elNum,fcNum) result(fcArea)
		integer :: fcNodes(3),elNodes(4)
		integer,intent(in) :: elNum,fcNum
		real(8) :: fcArea,fcVerts(3,3)

		fcNodes = getFaceNodes(fcNum)
		elNodes = meshElems(elNum)%nodes
		fcVerts = meshVerts(elNodes(fcNodes),:)
		fcArea = triangleArea(fcVerts)
	end function getFaceArea

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


	function getElementCentroid(elNum) result(centroid)
		integer :: elNodes(4)
		integer,intent(in) :: elNum
		real(8) :: centroid(3),ec(4,3)

		elNodes = meshElems(elNum)%nodes
		ec = meshVerts(elNodes,:)
		centroid = sum(ec,1)/4.d0
	end function getElementCentroid

	subroutine getElementUnitStiffness(elNum,elK,elemVol,elemUnitSt)
		integer,intent(in) :: elNum
		real(8),intent(in) :: elK
		real(8) :: ev,b(3,4),bt(4,3),spFns(4,4)
		real(8),intent(out) :: elemVol,elemUnitSt(4,4)

		spFns = getElementShapeFunctions(elNum)
		elemVol = getElementVolume(elNum)
		bt = spFns(:,2:4)
		b = transpose(bt)
		elemUnitSt = elK*matmul(bt,b)*elemVol
	end subroutine getElementUnitStiffness

!-----------------------------------------------------------------------!
!	End of the element level functions for the mesh
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

		if(propNameFlag == "T") then
			if(.not.(allocated(meshTemperatures))) then
				allocate(meshTemperatures(meshNumNodes))
			end if
			meshTemperatures = nodalVals
		elseif(propNameFlag == "S") then
			if(.not.(allocated(meshSources))) then
				allocate(meshSources(meshNumNodes))
			end if
			meshSources = nodalVals
		else
			write(*,'(a)') "Property to be set not recognised."
			stop
		end if

	end subroutine setMeshNodalValues

	subroutine saveMeshNodalStiffness(sySt,stColPtr,stRowPtr)
		integer,intent(in) :: stColPtr(:),stRowPtr(:)
		real(8),intent(in) :: sySt(:)

		if(.not.(allocated(meshStVals))) then
			allocate(meshStVals(size(sySt,1)))
		end if
		if(.not.(allocated(meshStColPtr))) then
			allocate(meshStColPtr(size(stColPtr,1)))
		end if
		if(.not.(allocated(meshStRowPtr))) then
			allocate(meshStRowPtr(size(stRowPtr,1)))
		end if

		meshStVals = sySt
		meshStColPtr = stColPtr
		meshStRowPtr = stRowPtr
	end subroutine saveMeshNodalStiffness

	subroutine setSurfaceConstTemperature(surfId,constT)
		integer :: i,j,el,fc,numFcs,fcNodes(3),elNodes(4)
		integer,intent(in) :: surfId
		real(8),intent(in) :: constT

		if(.not.(allocated(meshTemperatures))) then
			allocate(meshTemperatures(meshNumNodes))
		end if
		do i=1,meshNumSurfs
			if(meshSurfs(i)%sfId == surfId) then
				numFcs = meshSurfs(i)%numFcs
				do j=1,numFcs
					el = meshSurfs(i)%elNum(j)
					fc = meshSurfs(i)%fcNum(j)
					elNodes = meshElems(el)%nodes
					fcNodes = getFaceNodes(fc)
					meshTemperatures(elNodes(fcNodes)) = constT
				end do
				exit
			end if
		end do
	end subroutine setSurfaceConstTemperature

!-----------------------------------------------------------------------!
!	End of the setter subroutines
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
!	Subroutines to write elements of data processed previously.
!-----------------------------------------------------------------------!
	subroutine writeElementNeighbourhoodData()
		integer,parameter :: datFileNum=101
		character(72) :: datFile,wrtFmt

		datFile = commDatDir//trim(adjustl(meshFile))//commDatExt
		wrtFmt = '(1x,4(i8,1x,i4))'
		open(datFileNum,file=datFile)
		do i=1,meshNumElems
			write(datFileNum,wrtFmt) transpose(meshElems(i)%neighbours)
		end do
	end subroutine writeElementNeighbourhoodData
!-----------------------------------------------------------------------!
!	End of the writer subroutines
!-----------------------------------------------------------------------!

end module mesh
