module commonRoutines

	implicit none

	real(8),parameter :: pi=3.14159265358979d0

	contains

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

!-----------------------------------------------------------------------!
!	The following subroutines and functions are mathematical operations
!	required in one or multiple places in the modelling code. These may
!	sometimes be specific to the tetrahedral elements, however, most are
!	general routines.
!-----------------------------------------------------------------------!

	subroutine indexedSort(a,b)
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
	end subroutine indexedsort

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

	function getRaySphDir() result(dir)
		real(8) :: r1,r2,th,ph,dir(3)

		call random_number(r1)
		call random_number(r2)
		ph = 2*pi*r1
		th = acos(1.0d0 - 2.0d0*r2)
		dir(1) = sin(th)*cos(ph)
		dir(2) = sin(th)*sin(ph)
		dir(3) = cos(th)
	end function getRaySphDir

	function selTriPoint(fcVerts) result(pt)
		real(8) :: r1,r2,a,b,c,pt(3)
		real(8),intent(in) :: fcVerts(3,3)

		call random_number(r1)
		call random_number(r2)
		a = 1.0d0-sqrt(1-r1)
    	b = (1.0d0-a)*r2
		c = 1-a-b
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
			write(*,'(a)')"Face index not valid, must be between 1-4."
			stop
		end if
	end function getFaceNodes

	function triangleArea(trVerts) result(trArea)
		real(8) :: trArea,v(3),w(3),aVec(3)
		real(8),intent(in) :: trVerts(3,3)

		v = trVerts(2,:)-trVerts(1,:)
		w = trVerts(3,:)-trVerts(1,:)
		aVec = cross_product_3(v,w)
		trArea = 0.5d0*norm2(aVec)
	end function triangleArea

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

end module commonRoutines
