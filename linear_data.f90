module linear_data
    contains
        subroutine generate_data(x, y, sigmay, m, b, n_elements)
            integer, intent(in) :: n_elements, m, b
            real, intent(out) :: x(n_elements), y(n_elements), sigmay(n_elements)
            integer(kind=2) :: i
            real :: yrange
            
            call RANDOM_NUMBER(x)
            sigmay = gaussian_sampling(n_elements)
            x = x * 75 + 2
            
            y = m * x + b
            y_range = MAXVAL(y) - MINVAL(y)
            
            sigmay = sigmay * (0.1 * y_range + 2 * sigma / ABS(sigma))
            
            y = y + sigmay
            
        end subroutine generate_data

        function gaussian_sampling(n_elements) result(normal)
            integer, intent(in) :: n_elements
            real, dimension(n_elements) :: normal, uniform
            real, dimension(2) :: temp_n, temp_u
            integer(kind=2) :: i, half_elements
            
            half_elements = n_elements / 2
            
            call RANDOM_NUMBER(uniform)
            
            do i = 1, half_elements
                temp_u = (/uniform(i), uniform(i + half_elements) /)
                temp_n = (/normal(i), normal(i + half_elements) /)
                call box_miller(temp_u, temp_n)
                normal(i) = temp_n(1)
                normal(i + half_elements) = temp_n(2)
            end do
            
            
        end function gaussian_sampling

        subroutine box_miller(uniform, normal)
            real, intent(in) :: uniform(2)
            real, intent(out) :: normal(2)
            real(16), parameter :: PI_16 = 4 * atan (1.0_16)
             
            normal(1) = SQRT(-2 * LOG(uniform(1))) * COS(2 * PI_16 * uniform(2))
            normal(2) = SQRT(-2 * LOG(uniform(2))) * COS(2 * PI_16 * uniform(1))

            
        end subroutine box_miller
end module
