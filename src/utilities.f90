module utilities

	implicit none

	integer,parameter :: commLogErr=1,commLogInf=2,commLogDeb=3
	real(8),parameter :: pi=3.14159265358979d0,sigB=5.673d-8,			&
	LARGE=1.0d3,MEGA=1.0d6,GIGA=1.0d9,SMALL=1.0d-3,MICRO=1.0d-6,		&
	NANO=1.0d-9,ANGSTROM=1.d-10,PICO=1.0d-12,DPICO=1.d-13
	character(*),parameter :: commDatDir="../data/",					&
	commResDir="../results/",commMeshExt=".msh",commDatExt=".dat",		&
	commOutExt=".out",commVTKExt=".vtk"

	contains

!-----------------------------------------------------------------------!
!	The following routines work with intrinsic operations, mainly
!	as an aid to traceback and debugging.
!-----------------------------------------------------------------------!
    subroutine checkIoError(openStat,unitNum,message)
        integer,intent(in) :: openStat,unitNum
        character(*),intent(in),optional :: message
        
        if ((openStat > 0) .and. present(message)) then
            write(*,'(A)') "io-error: "//message
            close(unitNum)
            stop
        else if (openStat > 0) then
            write(*,'(a)') "error in IO read/write!"
            close(unitNum)
            stop
        end if
        return
    end subroutine checkIoError

	subroutine initRandom()
		integer,parameter :: rnFilNum=100
		integer :: i
		character(*),parameter :: rnFile="/dev/urandom"

		open(rnFilNum,file=rnFile,access='stream',form='UNFORMATTED')
		read(rnFilNum) i
		close(rnFilNum)
		call random_seed(put=(/i/))
	end subroutine initRandom

!-----------------------------------------------------------------------!
!	End of intrinsic procedure auxiliaries
!-----------------------------------------------------------------------!

!-----------------------------------------------------------------------!
!	The following subroutines and functions are mathematical operations
!	required in one or multiple places in the modelling code. These may
!	sometimes be specific to the tetrahedral elements, however, most are
!	general routines.
!-----------------------------------------------------------------------!

	subroutine solveQuadratic(a,b,c,x)
		real(8) :: q
		real(8),intent(in) :: a,b,c
		real(8),intent(out):: x(2)

		q = -(1.d0/2.d0)*(b+sign(sqrt(b**2.d0 - 4*a*c),b))
		x(1) = q/a
		x(2) = c/q
	end subroutine solveQuadratic

	subroutine getVolFracFromMassFrac(rho1,rho2,fm1,fv1)
		real(8) :: rhoRatio,mfInv
		real(8),intent(in) :: rho1,rho2,fm1
		real(8),intent(out) :: fv1

		rhoRatio = rho1/rho2
		mfInv = 1.0d0/fm1 - 1
		fv1 = 1.0d0/(rhoRatio*mfInv + 1)
	end subroutine getVolFracFromMassFrac

	function getPPVConstParticleSize(fV,dP) result(Nt)
		integer :: Nt
		real(8),intent(in) :: fV,dP

		Nt = nint(fV/((pi/6.0d0)*dP**3.0d0))
	end function getPPVConstParticleSize

	subroutine indexedSortInteger(a,b)
		integer,intent(inout) :: a(:)
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
	end subroutine indexedSortInteger

	subroutine invertReal4by4(vm)
		integer,parameter :: m=4,n=4,lda=4,lwork=256
		integer :: info,ipiv(4)
		real(8) :: work(256),vm(4,4)

		call dgetrf(m,n,vm,lda,ipiv,info)
		call dgetri(m,vm,lda,ipiv,work,lwork,info)
	end subroutine invertReal4by4

	function selTetraPoint(ec) result(pt)
		real(8) :: a,s,t,u,temp,pt(3)
		real(8),intent(in) :: ec(4,3)

		call random_number(s)
		call random_number(t)
		call random_number(u)
		if(s+t > 1.0d0) then
			s = 1.0d0 - s
			t = 1.0d0 - t
		end if
		if(s+t+u > 1.0d0) then
			if(t+u > 1.0d0) then
				temp = u
				u = 1.0d0 - s - t
				t = 1.0d0 - temp
			else
				temp = u
				u = s + t + u - 1.0d0
				s = 1.0d0 - t - temp
			end if
		end if
		a = 1.0d0-s-t-u
		pt = a*ec(1,:)+s*ec(2,:)+t*ec(3,:)+u*ec(4,:)
	end function selTetraPoint

!	function selTriPoint(fcVerts) result(pt)
!		real(8) :: r1,r2,a,b,c,pt(3)
!		real(8),intent(in) :: fcVerts(3,3)

!		call random_number(r1)
!		call random_number(r2)
!		a = 1.0d0-sqrt(1-r1)
!    	b = (1.0d0-a)*r2
!		c = 1-a-b
!		pt = a*fcVerts(1,:)+b*fcVerts(2,:)+c*fcVerts(3,:)
!	end function selTriPoint

	function selTriPoint(fcVerts) result(pt)
		real(8) :: r1,r2,a,b,c,pt(3)
		real(8),intent(in) :: fcVerts(3,3)

		call random_number(r1)
		call random_number(r2)
		do while((r1+r2).gt. 1.d0)
			r1 = 1.d0 - r1
			r2 = 1.d0 - r2
			if((r1+r2).gt. 1.d0) then
				call random_number(r1)
				call random_number(r2)
			end if
		end do
		a = 1.d0 - r1 - r2
		b = r1
		c = r2 
		pt = a*fcVerts(1,:)+b*fcVerts(2,:)+c*fcVerts(3,:)
	end function selTriPoint

	function getFaceNorm(fcVerts,remVert,inward) result(fcNorm)
		integer :: i
		real(8),parameter :: small=0.000000000001d0
		real(8) :: remPtDist,edge1(3),edge2(3),fcNorm(3)
		real(8),intent(in) :: fcVerts(3,3),remVert(3)
		logical,intent(in) :: inward

		edge1 = fcVerts(2,:)-fcVerts(1,:)
		edge2 = fcVerts(3,:)-fcVerts(1,:)
		fcNorm = cross_product_3(edge1,edge2)
		fcNorm = fcNorm/norm2(fcNorm)
		remPtDist = dot_product(fcNorm,(remVert-fcVerts(1,:)))
		if(inward) then
			if(remPtDist<0) then
				do i=1,3
					if(abs(fcNorm(i)) .gt. small) then
						fcNorm(i) = -fcNorm(i)
					end if
				end do
			end if
		else
			if(remPtDist>0) then
				do i=1,3
					if(abs(fcNorm(i)) .gt. small) then
						fcNorm(i) = -fcNorm(i)
					end if
				end do
			end if
		end if
	end function getFaceNorm

	function getFaceIndex(fcNodes) result(fcNum)
		integer :: fcNum,indices(3)
		integer,intent(inout) :: fcNodes(3)

		call indexedSortInteger(fcNodes,indices)
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
			write(*,'(a)')"Face index not valid, must be between 1-4."
			stop
		end if
	end function getFaceNodes

    function insideFaceCheck(fcVerts,pt) result(ptInside)
		integer :: i
        real(8),intent(in) :: pt(3),fcVerts(3,3)
        logical:: ptInside
        real(8) :: dpna,na(3),nx(3),ny(3),nz(3),bc(3)
        
        na = cross_product_3(fcVerts(2,:)-fcVerts(1,:),					&
		fcVerts(3,:)-fcVerts(1,:))
        nx = cross_product_3(fcVerts(3,:)-fcVerts(2,:),pt-fcVerts(2,:))
        ny = cross_product_3(fcVerts(1,:)-fcVerts(3,:),pt-fcVerts(3,:))
        nz = cross_product_3(fcVerts(2,:)-fcVerts(1,:),pt-fcVerts(1,:))
        dpna = dot_product(na,na)
        bc(1) = dot_product(na,nx)/dpna
        bc(2) = dot_product(na,ny)/dpna
        bc(3) = dot_product(na,nz)/dpna
        ptInside = (all(bc.gt. 0.0d0).and.(abs(sum(bc)-1.0d0).le.PICO))
!        if (.not.(ptInside)) then
!			write(*,*)"Point found not lying on face. Values: "
!			write(*,*) "point: "
!			write(*,'(3(f15.12,2x))') pt
!			write(*,*) "fcVerts: "
!			do i=1,3
!				write(*,'(3(f15.12,2x))') fcVerts(i,:)
!			end do
!			write(*,*) "bc: "
!	        write(*,'(3(f15.12,2x))') bc
!	        write(*,*) abs(sum(bc)-1.0d0)
!        end if
    end function insideFaceCheck

	function triangleArea(trVerts) result(trArea)
		real(8) :: trArea,v(3),w(3),aVec(3)
		real(8),intent(in) :: trVerts(3,3)

		v = trVerts(2,:)-trVerts(1,:)
		w = trVerts(3,:)-trVerts(1,:)
		aVec = cross_product_3(v,w)
		trArea = 0.5d0*norm2(aVec)
	end function triangleArea

	function areaIntShapeFuncs(fcNodes) result(surfint)
		integer :: i,j,n1,n2
		integer,intent(in) :: fcNodes(3)
		real(8) :: surfint(4,4)

		surfint = 0.0d0
		do i=1,3
			n1 = fcnodes(i)
			do j=1,3
				n2 = fcnodes(j)
				if(n1==n2) then
					surfint(n1,n2) = 2.0d0
				else
					surfint(n1,n2) = 1.0d0
				end if
			end do
		end do
	end function areaIntShapeFuncs

	subroutine createLocalCS(s,nx,ny)
		real(8),parameter :: xAx(3) = (/1.d0,0.d0,0.d0/), 					&
		yAx(3) = (/0.d0,1.d0,0.d0/), zAx(3) = (/0.d0,0.d0,1.d0/)
		real(8),intent(in) :: s(3)
		real(8),intent(out) :: nx(3),ny(3)
		real(8) :: v(3)

		v = s/norm2(s)	! Just in case it is not normalised
		ny = cross_product_3(v,zAx)
		if(norm2(ny) .lt. 1.d-9) then
			nx = xAx
			ny = cross_product_3(v,nx)
		else
			ny = ny/norm2(ny)
			nx = cross_product_3(ny,v)
		end if
	end subroutine createLocalCS

	function cross_product_3(v,w) result(cP)
		real(8) :: cP(3)
		real(8),intent(in) :: v(3),w(3)

		cP(1) = v(2)*w(3)-w(2)*v(3)
		cP(2) = -(v(1)*w(3)-w(1)*v(3))
		cP(3) = v(1)*w(2)-w(1)*v(2)
	end function cross_product_3

	function determinantReal4by4(A) result(d)
		real(8) :: d,A(4,4)

		d = A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+				&
			A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+						&
			A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))-						&
			A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+				&
			A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+ 						&
			A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))+						&
			A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+				&
			A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+						&
			A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))-						&
			A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+ 				&
			A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+						&
			A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))

	end function determinantReal4by4

	subroutine cubicSplineNaturalFit(x,y,yy)
		integer :: i,n
		real(8),allocatable :: rhs(:),temp(:),lhs(:,:)
		real(8),intent(in) :: x(:),y(:)
		real(8),allocatable,intent(out) :: yy(:)

		n = size(x,1)
		allocate(lhs(n,n))
		allocate(rhs(n))
		allocate(temp(n))
		allocate(yy(n-2))
		lhs = 0.0d0
		rhs = 0.0d0
		do i=2,n-1
			lhs(i,i-1) = (x(i)-x(i-1))/(6.d0)
			lhs(i,i) = (x(i+1)-x(i-1))/(3.d0)
			lhs(i,i+1) = (x(i+1)-x(i))/(6.d0)
			rhs(i) = (y(i+1)-y(i))/(x(i+1)-x(i)) - 						&
			(y(i)-y(i-1))/(x(i)-x(i-1))
		end do

		call triDiagonalSolver(lhs(2:n-1,2:n-1),rhs(2:n-1),yy)
		temp = 0.0d0
		temp(2:n-1) = yy
		call move_alloc(temp,yy)
	end subroutine cubicSplineNaturalFit

	subroutine triDiagonalSolver(lhs,rhs,sols)
		integer :: i,n
		real(8),intent(in) :: lhs(:,:),rhs(:)
		real(8),allocatable :: l(:,:),r(:)
		real(8),allocatable,intent(out) :: sols(:)

		n = size(rhs,1)
		allocate(r(n))
		allocate(l(n,n))
		allocate(sols(n))
		l = lhs
		r = rhs
		do i=1,n-1
			if(i==1) then
				l(i,i+1) = l(i,i+1)/l(i,i)
				r(i) = r(i)/l(i,i)
			else
				l(i,i+1) = l(i,i+1)/(l(i,i)-l(i,i-1)*l(i-1,i))
				r(i) = (r(i)-l(i,i-1)*r(i-1))
				r(i) = r(i)/(l(i,i)-l(i,i-1)*l(i-1,i))
			end if
		end do
		r(n) = (r(n)-l(n,n-1)*r(n-1))/(l(n,n)-l(n,n-1)*l(n-1,n))
		sols(n) = r(n)
		do i=n-1,1,-1
			sols(i) = r(i) - l(i,i+1)*sols(i+1)
		end do
		deallocate(r)
		deallocate(l)
	end subroutine triDiagonalSolver

	subroutine skipReadLines(unitNum,numLines)
		integer :: i
		integer,intent(in) :: unitNum,numLines

		do i=1,numLines
			read(unitNum,*)
		end do
	end subroutine skipReadLines

	subroutine getFileNumLines(unitNum,nLines)
		integer :: error
		integer,intent(in) :: unitNum
		integer,intent(out) :: nLines

		nLines = 0
		do
			read(unitNum,*,iostat = error)
			if (error .lt. 0) exit
			nLines = nLines + 1
		end do
		close(unitNum)
	end subroutine getFileNumLines

	subroutine bisLocReal(A,val,ind1,ind2)
		integer :: k,st,en,mid
		real(8),intent(in) :: val,A(:)
		integer,intent(out) :: ind1,ind2
		logical :: search

		k = size(A,1)
		if(val .lt. A(1)) then
			ind1 = 1
			ind2 = 1
			return
		end if
		if(val .gt. A(k-1)) then
			ind1 = k-1
			ind2 = k
			return
		end if
		st = 1
		mid = k/2
		en = k
		search = .true.
		do while(search)
			if(val < A(mid)) then
				en = mid
				mid = (st + en)/2
				if((A(mid).gt.val).and.(A(mid-1).lt.val)) then
					ind1 = mid-1
					ind2 = mid
					search = .false.
				end if
			else
				st = mid
				mid = (st + en)/2
				if((A(mid).lt.val).and.(A(mid+1).gt.val)) then
					ind1 = mid
					ind2 = mid+1
					search = .false.
				end if				
			end if
		end do
	end subroutine

!-----------------------------------------------------------------------!
!	End of mathematical operation routines
!-----------------------------------------------------------------------!

end module utilities
