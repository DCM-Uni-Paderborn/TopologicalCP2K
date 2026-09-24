program signature_reallocation
  implicit none
  integer, allocatable :: matrix(:,:), weights(:), expected(:)
  integer :: p, a, b
  allocate(matrix(3,1), weights(1))
  matrix(:,1) = [2,1,1]
  weights = [1]
  do p=1,2
    a=p
    b=2*p-1
    expected=matmul(matrix(a:b,:),weights)
    print *, p, a, b, size(expected), expected
    if (size(expected)/=b-a+1) error stop 'Result extent incorrect'
  end do
end program signature_reallocation
