program runModel

	use rt, only: rtSimple

	implicit none

	integer :: i,values1(8),values2(8),nHood(4,2)
	character :: date*8,time*10,zone*5

	call date_and_time(date,time,zone,values1)
	write(*,'(a,(3i8,2x),i8)') "Start time: ", values1(5:8)
	call rtSimple()
	call date_and_time(date,time,zone,values2)
	write(*,'(a,(3i8,2x),i8)') "End time: ",values2(5:8)
	write(*,'(a,(3i8,2x),i8)') "Runtime: ",values2(5:8)-values1(5:8)
end program runModel
