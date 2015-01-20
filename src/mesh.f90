! Mesh module: contains everything grid-related

module mesh

	use utilities, only: commDatDir,commResDir,commMeshExt,commDatExt,	&
						 indexedSortInteger,checkIoError,getFaceIndex,	&
						 getFaceNodes,invertReal4by4,triangleArea,		&
						 determinantReal4by4
							  

	implicit none

!	Declare all relevant types
	type tetraElement
		integer :: domain,nodes(4),neighbours(4,3)
		real(8) :: volume,centroid(3)
	end type tetraElement

	type surface
		character(len=16) :: sfName
		integer :: sfId,numFcs
		integer,allocatable :: elNum(:),fcNum(:)
		real(8) :: totalArea
		real(8),allocatable :: fcArea(:)
		logical :: iFace
	end type surface

	type elementBin
		real(8) :: limLow(3),limHigh(3)
		integer,allocatable :: bin(:)
	end type elementBin

!	Module variables
	integer :: meshNumNodes,meshNumElems,meshNumDoms,meshNumSurfs,		&
	meshNumBys,meshNumIntfcs
	integer,allocatable :: meshBys(:),meshIntfcs(:),meshStColPtr(:),	&
	meshStRowPtr(:)
	real(8),allocatable :: meshTemperatures(:),meshSources(:),			&
	meshVerts(:,:),meshStVals(:)
	type(tetraElement),allocatable :: meshElems(:)
	type(surface),allocatable :: meshSurfs(:)
	type(elementBin),allocatable :: meshElBins(:,:,:)
!	Inputs
	integer :: meshNBins(3) = (/2,2,2/)
	character(72) :: meshFilePre = "a"
	character(*),parameter :: meshNHFileSuff = "_nHood",				&
							  meshBinFileSuff = "_bins"

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

		meshFileName = commDatDir//trim(adjustl(meshFilePre))//commMeshExt
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
!			l = (verify('iface',surfName) .ne. 0)
			l = (count(i==meshIntfcs) .gt. 0)
			meshSurfs(i)%iFace = l
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
							meshElems(elNum)%neighbours(fcNum,3) = i
						else
							meshElems(elNum)%neighbours(fcNum,1) = elNum
							meshElems(elNum)%neighbours(fcNum,2) = -1
							meshElems(elNum)%neighbours(fcNum,3) = i
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
		integer,parameter :: nNhFile=101,nBinFile=102
		integer :: i,j,k,binSize,empBinCt,numTotBins,tempNBins(3),		&
		temp(8),values(8)
		logical :: nhFileExist,binFileExist
		character(72) :: nhFile,binFile
		character :: date*8,time*10,zone*5

		call date_and_time(date,time,zone,values)
		write(*,'(a,1x,4(i4))')"Neighbourhood gathering started: ",values(5:8)		
		nhFile = trim(adjustl(meshFilePre))//meshNHFileSuff
		nhFile = commDatDir//trim(adjustl(nhFile))//commDatExt
		inquire(file=trim(adjustl(nhFile)),exist=nhFileExist)
		binFile = trim(adjustl(meshFilePre))//meshBinFileSuff
		binFile = commDatDir//trim(adjustl(binFile))//commDatExt
		inquire(file=trim(adjustl(binFile)),exist=binFileExist)
		numTotBins = product(meshNBins)
		if(nhFileExist.and.binFileExist) then
			open(nBinFile,file=binfile)
			read(nBinFile,*) tempNBins
			close(nBinFile)
			if(all(tempNBins .eq. meshNBins)) then
				write(*,'(a)')"All neighbourhood data exists. Reading &
				&files now ..."
				call readElementBins(nBinFile,binFile,empBinCt)
			else
				goto 100
			end if
			call readElementNeighbours(nNhFile,nhFile)
			call date_and_time(date,time,zone,values)
			write(*,'(a,1x,4(i4))')"Neighbourhood gathering completed: ",values(5:8)
			return
		elseif(nhFileExist.and.(.not.binFileExist)) then
100			write(*,'(a)')"Legacy neighbourhood data found. New bins &
			&needed."
			call binElements(empBinCt)
			call writeElementNeighbourhoodData('bins')
			call readElementNeighbours(nNhFile,nhFile)
			return
		elseif(binFileExist.and.(.not.nhFileExist)) then
			open(nBinFile,file=binfile)
			read(nBinFile,*) tempNBins
			close(nBinFile)
			if(all(tempNBins .eq. meshNBins)) then
				write(*,'(a)')"Legacy bins found."
				call readElementBins(nBinFile,binFile,empBinCt)
				goto 200
			else
				goto 150
			end if
			return
		else
			write(*,'(a)')"No neighbourhood data found. Generating ..."
			goto 150
			return
		end if

150		call binElements(empBinCt)
		call writeElementNeighbourhoodData('bins')
		
200		do k=1,meshNBins(3)
			do j=1,meshNBins(2)
				do i=1,meshNBins(1)
					binSize = size(meshElBins(i,j,k)%bin,1)
					if(binSize .eq. 0) cycle
					call findNeighboursWithinBin(i,j,k)
					call findNeighboursAcrossBins(i,j,k)
				end do
			end do
		end do
		call date_and_time(date,time,zone,values)
		call writeElementNeighbourhoodData('nhbr')
		write(*,*)"Total number of bins assigned: ",numTotBins
		write(*,*)"Number of empty bins found: ",empBinCt
		write(*,'(a,1x,4(i4))')"Neighbourhood gathering completed: ",values(5:8)
	end subroutine getElementNeighbours

	subroutine readElementNeighbours(nNhFile,nhFile)
		integer,intent(in) :: nNhFile
		integer :: i,j,k,temp(12)
		character(*),intent(in) :: nhFile
		
		open(nNhFile,file=nhFile)
		do i=1,meshNumElems
			read(nNhFile,'(1x,4(i8,1x,i4,1x,i4))') temp
			do j=1,4
				meshElems(i)%neighbours(j,:) = temp(3*j-2:3*j)
			end do
		end do
		close(nNhFile)
	end subroutine readElementNeighbours

	subroutine readElementBins(nBinFile,binFile,empBinCt)
		integer,intent(in) :: nBinFile
		integer,intent(out) :: empBinCt
		integer :: i,j,k,a,b,c,binSz,m,n,p
		character(*),intent(in) :: binFile

		open(nBinFile,file=binFile)
		read(nBinFile,*) meshNBins
		if(.not.(allocated(meshElBins))) then
			allocate(meshElBins(meshNBins(1),meshNBins(2),meshNBins(3)))
		else
			deallocate(meshElBins)
			allocate(meshElBins(meshNBins(1),meshNBins(2),meshNBins(3)))
		end if
		empBinCt = 0
		do k=1,meshNBins(3)
			do j=1,meshNBins(2)
				do i=1,meshNBins(1)
					read(nBinFile,*) a,b,c,binSz
					read(nBinFile,*) meshElBins(i,j,k)%limLow,	&
					meshElBins(i,j,k)%limHigh
					if(binSz .eq. 0) then
						empBinCt = empBinCt+1
						cycle
					end if
					allocate(meshElBins(i,j,k)%bin(binSz))
					n = binSz
					p=0
					do while(n.ne.0)
						m=min(n,10)
						read(nBinFile,*)meshElBins(i,j,k)%bin(p+1:p+m)
						p=p+m
						n=n-m
					enddo				
				end do
			end do
		end do
		close(nBinFile)		
	end subroutine readElementBins

	subroutine binElements(empBinCt)
		integer :: i,j,k,sz,maxPerBin,endPt,check,cRs(3),elNodes(4)
		integer,intent(out) :: empBinCt
		integer,allocatable :: temp(:),lZ(:,:,:)
		real(8) :: dmin(3),dmax(3),edges(3),elCent(3),binEdges(3),		&
		elVerts(4,3)

		if(.not.(allocated(meshElBins))) then
			allocate(meshElBins(meshNBins(1),meshNBins(2),meshNBins(3)))
		end if
		allocate(lZ(meshNBins(1),meshNBins(2),meshNBins(3)))
		lZ = 0
		dmin = minval(meshVerts,1)
		dmax = maxval(meshVerts,1)
		edges = dmax-dmin
		binEdges = edges/(real(meshNbins,8))
		maxPerBin = 4*meshNumElems/product(meshNBins)

		do i=1,meshNumElems
			elNodes = meshElems(i)%nodes
			elVerts = meshVerts(elNodes,:)
			elCent = sum(elVerts,1)/4.0d0
			meshElems(i)%centroid = elCent
			cRs = ceiling(((elCent-dmin)/edges)*real(meshNBins,8))
			where(cRs == 0)
				cRs = 1
			end where
			sz = size(meshElBins(cRs(1),cRs(2),cRs(3))%bin,1)
			if(sz.gt.0) then
				lZ(cRs(1),cRs(2),cRs(3)) = lZ(cRs(1),cRs(2),cRs(3)) + 1
				meshElBins(cRs(1),cRs(2),cRs(3))%bin(lZ(cRs(1),cRs(2),cRs(3))) = i
			else
				allocate(meshElBins(cRs(1),cRs(2),cRs(3))%bin(4*maxPerBin))
				meshElBins(cRs(1),cRs(2),cRs(3))%bin = 0
				meshElBins(cRs(1),cRs(2),cRs(3))%bin(1) = i
				lZ(cRs(1),cRs(2),cRs(3)) = 1
			end if
		end do

		check = 0
		do i=1,meshNBins(1)
			do j=1,meshNBins(2)
				do k=1,meshNBins(3)
					meshElBins(i,j,k)%limLow = real((/i-1,j-1,k-1/),8)*binEdges
					meshElBins(i,j,k)%limHigh = real((/i,j,k/),8)*binEdges
					if(lZ(i,j,k) .eq. 0) then
						empBinCt = empBinCt + 1
						cycle
					end if
					allocate(temp(lZ(i,j,k)))
					check = check + lZ(i,j,k)
					temp = meshElBins(i,j,k)%bin(1:lZ(i,j,k))
					deallocate(meshElBins(i,j,k)%bin)
					call move_alloc(temp,meshElBins(i,j,k)%bin)
				end do
			end do
		end do
		write(*,'(i8,1x,a)')check, "elements binned"
	end subroutine binElements

	subroutine findNeighboursWithinBin(xin,yin,zin)
		integer,intent(in) :: xin,yin,zin
		integer :: i,j,binSize,ct1,ct2,el1,el2,fc1,fc2,ref,elNo1(4),	&
		elNo2(4)
		logical :: shared
		type(elementBin) :: singleBin

		singleBin = meshElBins(xin,yin,zin)
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

	subroutine findNeighboursAcrossBins(c1,c2,c3)
		integer,intent(in) :: c1,c2,c3
		integer :: i,j,k,ol,il,ct1,ct2,binSz1,binSz2,el1,el2,fc1,fc2,	&
		elNo1(4),elNo2(4)
		integer,allocatable :: givenBin(:),adjBin(:)
		logical :: shared

		binSz1 = size(meshElBins(c1,c2,c3)%bin,1)
		allocate(givenBin(binSz1))
		givenBin = meshElBins(c1,c2,c3)%bin
		currentBin: do ol=1,binSz1
			el1 = givenBin(ol)
			ct1 = count(meshElems(el1)%neighbours(:,1)>0)
			if(ct1 == 4) cycle
			elNo1 = meshElems(el1)%nodes
			zBinIndex: do k=c3-1,c3+1
				if((k==0).or.(k.gt.meshNBins(3))) cycle zBinIndex
				yBinIndex: do j=c2-1,c2+1
					if((j==0).or.(j.gt.meshNBins(2))) cycle yBinIndex
					xBinIndex: do i=c1-1,c1+1
						if((i==0).or.(i.gt.meshNBins(1))) cycle xBinIndex
						if(all((/i,j,k/)==(/c1,c2,c3/))) cycle xBinIndex
						binSz2 = size(meshElBins(i,j,k)%bin,1)
						if(binSz2 .eq. 0) then
							cycle xBinIndex
						else
							if(allocated(adjBin)) deallocate(adjBin)
							allocate(adjBin(binSz2))
						end if
						adjBin = meshElBins(i,j,k)%bin
						adjacentBin: do il=1,binSz2
							shared = .false.
							el2 = adjBin(il)
							ct2 = count(meshElems(el2)%neighbours(:,1)>0)
							if(ct2 == 4) cycle adjacentBin
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

	subroutine shapeFunctionsAtPoint(elNum,pt,spFnVals)
		integer,intent(in) :: elNum
		real(8) :: spFns(4,4)
		real(8),intent(in) :: pt(3)
		real(8),intent(out) :: spFnVals(4)

		spFns = getElementShapeFunctions(elNum)
		spFnVals = (/1.0d0,pt/)
		spFnVals = matmul(spFns,spFnVals)
	end subroutine shapeFunctionsAtPoint

	subroutine getPointBin(pt,numBins,dmin,dmax,pointBin)
		integer,intent(in) :: numBins(3)
		integer,intent(out) :: pointBin(3)
		real(8),intent(in) :: pt(3),dmin(3),dmax(3)
		real(8) :: edges(3)

		edges = dmax-dmin
		pointBin = ceiling(((pt-dmin)/edges)*real(numBins,8))
		where(pointBin == 0)
			pointBin = 1
		end where
	end subroutine getPointBin
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
	subroutine writeElementNeighbourhoodData(wrJob)
		character(*),intent(in) :: wrJob

		if(wrJob .eq. 'bins') then
			call writeElementBins()
		elseif(wrJob .eq. 'nhbr') then
			call writeElementNeighbours()
		elseif(wrJob .eq. 'both') then
			call writeElementBins()
			call writeElementNeighbours()
		else
			write(*,'(a)')"Wrong write job specification, please review."
		end if

	end subroutine writeElementNeighbourhoodData

	subroutine writeElementNeighbours()
		integer,parameter :: nNhFile=101
		integer :: i
		character(*),parameter :: wFmt='(1x,4(i8,1x,i4,1x,i4))'
		character(72) :: nhFile

		nhFile = trim(adjustl(meshFilePre))//meshNHFileSuff
		nhFile = commDatDir//trim(adjustl(nhFile))//commDatExt
		open(nNhFile,file=nhFile)
		do i=1,meshNumElems
			write(nNhFile,wFmt) transpose(meshElems(i)%neighbours)
		end do
		close(nNhFile)
	end subroutine writeElementNeighbours

	subroutine writeElementBins()
		integer,parameter :: nBinFile=102
		integer ::  i,j,k,binSz,n,p,m
		character(*),parameter :: wFmt='(10(1x,i8))'
		character(72) :: binFile

		binFile = trim(adjustl(meshFilePre))//meshBinFileSuff
		binFile = commDatDir//trim(adjustl(binFile))//commDatExt
		open(nBinFile,file=binFile)
		write(nBinFile,'(3(i4,2x))') meshNBins
		do k=1,meshNBins(3)
			do j=1,meshNBins(2)
				do i=1,meshNBins(1)
					binSz = size(meshElBins(i,j,k)%bin,1)
					write(nBinFile,'(3(i4,2x),i8)') i,j,k,binSz
					write(nBinFile,'(6(e14.6,2x))')meshElBins(i,j,k)%limLow,&
					meshElBins(i,j,k)%limHigh
					if(binSz .eq. 0) cycle
					n = binSz
					p=0
					do while(n.ne.0)
						m=MIN(n,10)
						write(nBinFile,wFmt)meshElBins(i,j,k)%bin(p+1:p+m)
						p=p+m
						n=n-m
					enddo				
				end do
			end do
		end do
		close(nBinFile)
	end subroutine writeElementBins
!-----------------------------------------------------------------------!
!	End of the writer subroutines
!-----------------------------------------------------------------------!

end module mesh
