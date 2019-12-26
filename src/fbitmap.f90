

subroutine symmetric_jaccard_kernel(n_items, vec1, kernel)

    ! lengths = np.sum(vectors, axis=1)
    ! shape = lengths.shape[0]
    ! lengths = lengths.reshape((shape, 1))
    !
    ! kernel = np.dot(vectors, vectors.T)
    ! kernel = kernel / (lengths + lengths.T - kernel)

    implicit none

    integer, intent(in) :: n_items
    double precision, intent(in), dimension(:,:) :: vec1
    double precision, intent(out), dimension(n_items,n_items) :: kernel
    integer :: i, j
    double precision :: dp

    ! n_items = shape(vec1)

    integer, dimension(:), allocatable :: lengths

    kernel = 0.0

    write(*,*) "jck", n_items

    allocate(lengths(n_items))

    !$omp parallel do
    do i = 1, n_items
        lengths(i) = sum(vec1(i,:))
    end do
    !$omp end parallel do

    !$omp parallel do private(dp)
    do j = 1, n_items
        do i = j, n_items

            dp = DOT_PRODUCT(vec1(i,:), vec1(j,:))

            kernel(i,j) = dp / (lengths(i) + lengths(j) - dp)
            kernel(j,i) = kernel(i,j)

        end do
    end do
    !$omp end parallel do

    deallocate(lengths)

end subroutine symmetric_jaccard_kernel

