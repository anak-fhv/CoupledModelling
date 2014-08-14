module commonRoutines

	implicit none

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
