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
		integer :: sfNumber
		integer,allocatable :: elNum(:),fcNum(:)
	end type surface

	contains

!	Module variables
	integer :: meshNumNodes,meshNumElems,meshNumDoms,meshNumSurfs
	integer,allocatable :: meshConnTable(:,:)
	real(8),allocatable :: meshT(:),meshF(:),meshVertices(:,:)
	type(tetraElement),allocatable :: meshElems(:)
	type(surface),allocatable :: meshSurfs(:)
	character(*),intent(in) :: meshFile


!	Main mesh reading routine. Note that the entire file is read
!	serially and variables are populated as they are available,
!	since carriage control to read specific data is quite complex
!	for the mesh file format
	subroutine readMesh()
		integer,parameter :: fNo=100
		integer :: openStat,mDets(7)

        call openmeshfile(fno,meshFile)
        call readmeshdetails(fno,mDets)
		meshNumNodes = mDets(1)
		meshNumElems = mDets(2)
		meshNumDoms = mDets(6)
		meshNumSurfs = mDets(7)
        call readmeshvertices(fno,meshNumNodes,meshVertices)
        call readmeshconnectivity(fno,meshNumElems,meshElems)
        call readmeshdomains(fno,meshNumDoms,meshElems)
        call readmeshsurfaces(fno,meshNumSurfs,meshSurfs)
        call closemeshfile(fno)
	end subroutine readMesh

    subroutine openmeshfile(fno,meshFile)
        integer :: fno
        character(*),intent(in) :: meshFile

        open(unit=fno,file=fname,form='formatted',status='old',		&
		action='read',iostat=openStat)
		call check_io_error(openStat)
    end subroutine openmeshfile



end module mesh
