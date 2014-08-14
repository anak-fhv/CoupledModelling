module rtHelper

	use commonRoutines
	implicit none

	contains

	subroutine getNextFace(ec,pt,dir,newFc,lToFc)
		integer :: fc,remNo,newFc,fcNodes(3)
		integer,intent(out) :: newFc
		real(8) :: ml,remVert(3),fcNorm(3),cent(3),fcVerts(3,3)
		real(8),intent(in) :: pt(3),dir(3),ec(4,3)
		real(8),intent(out) :: lToFc
		logical :: inward

		inward = .true.
		ml = 1000.0d0
		do fc =1,4
			fcNodes = getFaceNodes(fc)
			remNode = 10-sum(fcNodes)
			fcVerts = ec(fcNodes,:)
			remVert = ec(remNode,:)
			fcNorm = getFaceInwardNorm(fcVerts,remVert,inward)
			ptToFc = dot_product(fcNorm,(pt-fcVerts(1,:)))
			lToFc = ptToFc/dot_product(dir,-fcNorm)
			if(lToFc < 0) cycle
			if (lToFc < ml) then
    			ml = lToFc
    			newFc = fc
			end if
		end do
		lToFc = ml
	end subroutine getNextFace

end module rtHelper
