subroutine setup(pdbin,maxd,n,diff,distance,connect)
integer, intent(in) :: n
real(8), intent(in) :: pdbin(n,3)
real(8), intent(in) :: maxd
integer i,j
real(8), intent(out) :: diff(n,n,3)
real(8), intent(out) :: distance(n,n)
real(8), intent(out) :: connect(n,n)

!----------------------------------------
! calculate diffrence and distance between residues i and j
! connect two nodes if close enough

do i=1,n
  do j=1,n
    diff(i,j,:)=pdbin(i,:)-pdbin(j,:)
    distance(i,j)=sqrt(dot_product(diff(i,j,:),diff(i,j,:)))
    if(distance(i,j)<maxd) then
      connect(i,j)=-1
    else 
      connect(i,j)=0
    end if
  end do
end do
end subroutine


subroutine makeHess(diff,distance,connect,HESS,n)
integer, intent(in) :: n
real(8), intent(in) :: diff(n,n,3)
real(8), intent(in) :: distance(n,n)
real(8), intent(in) :: connect(n,n)
real(8), intent(out) :: HESS(3*n,3*n)
integer i,j,k,l

!----------------------------------------
! Compute off-diagonal super elements of HESS
! super elemnts are 3x3 matrics

do i=1,n
  do j=1,n
    do k=1,3
      do l=1,3
        if(i == j) then
          HESS(3*(i-1)+k,3*(j-1)+l)=0.0d0
        else
          HESS(3*(i-1)+k,3*(j-1)+l)=dble(connect(i,j)) &
          *diff(i,j,k)*diff(i,j,l)/distance(i,j)
        end if
        !         write(*,"(4i5,f10.4)") i,j,k,l,HESS(3*(i-1)+k,3*(j-1)+l)
      end do
    end do
  end do
end do

!----------------------------------------
! Compute diagonal super elements of HESS
! super elemnts are 3x3 matrics

do i=1,n
  do j=1,n
    do k=1,3
      do l=1,3
        if( i /=j) then
          HESS(3*(i-1)+k,3*(i-1)+l)= &
          HESS(3*(i-1)+k,3*(i-1)+l)-HESS(3*(i-1)+k,3*(j-1)+l)
        endif
        !         write(*,"(4i5,f10.4)") i,j,k,l,HESS(3*(i-1)+k,3*(j-1)+l)
      enddo
    enddo
  enddo
enddo

end subroutine


subroutine calcBfactor(eig_vec,eig_val,neig,nres,B)
real(8), intent(in) :: eig_vec(neig,nres,3)
real(8), intent(in) :: eig_val(neig)
integer, intent(in) :: neig,nres
real(8), intent(out) :: B(nres)
integer i,j

B(:) = 0

do i = 1,nres
  do j = 1,neig
    B(i) = B(i) + dot_product(eig_vec(j,i,:),eig_vec(j,i,:)) / eig_val(j)
  end do
end do

end subroutine
