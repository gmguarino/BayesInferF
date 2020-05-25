program LinearData
    
    use linear_data
    
    implicit none

    
    integer :: i, m, b, n_elements
    real, dimension(6) :: x, y, sigmay
    
    m = 2
    b = 0
    n_elements = 6
    
    call generate_data(x, y, sigmay, m, b, n_elements)
    
    do i = 1,n_elements
        print *, x(i), " ", y(i), " ", sigmay(i)
    end do
    
end program


