
module data_generator
    contains
    subroutine generate_data(x, y, sigmay, m, b, n_elements)
        integer, intent(in) :: n_elements
        real, intent(in) :: m, b
        real, intent(out) :: x(n_elements), y(n_elements), sigmay(n_elements)

        real :: y_range
        
        call RANDOM_NUMBER(x)
        x = x * 75 + 2
        
        y = m * x + b
        y_range = MAXVAL(y) - MINVAL(y)
        
        sigmay = gaussian_sampling(n_elements) * 0.1 * y_range 
        sigmay = sigmay + 2 * (sigmay / ABS(sigmay))
        
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

module probability

    real(16), parameter :: PI_16 = 4 * atan (1.0_16)

    contains
    
    subroutine log_likelihood(parameters, x, y, error, n_elements, res)
        integer, intent(in) :: n_elements
        real, intent(in) :: parameters(2), x(n_elements), y(n_elements), error(n_elements)
        real, intent(out) :: res

        res = SUM(LOG(1/(SQRT(2 * PI_16) * ABS(error)))) + &
                SUM(-0.5 * (y - (parameters(1) * x + parameters(2)))**2 / (error ** 2))

    end subroutine log_likelihood
    
    subroutine ln_uniform_prior_dist(parameters, limits, res)
        real, intent(in) :: parameters(2), limits(2, 2)
        real, intent(out) :: res
        real :: prior_m, prior_b, m, b, mlimits(2), blimits(2)
        
        m = parameters(1)
        b = parameters(2)
        
        mlimits = limits(1, :)
        blimits = limits(2, :)
        
        if ((m >= mlimits(1)).OR.(m <= mlimits(2))) then
            prior_m = LOG(1.0 / (mlimits(2) - mlimits(1)))
        else
            prior_m = (- HUGE(0.0))
        end if
            
        if ((b >= blimits(1)).OR.(b <= blimits(2))) then
            prior_b = LOG(1.0 / (blimits(2) - blimits(1)))
        else
            prior_b = (- HUGE(0.0))
        end if
        
        res = prior_m + prior_b
    
    end subroutine ln_uniform_prior_dist
    
    subroutine ln_posterior_dist(parameters, x, y, error, limits, n_elements, res)
        
        integer, intent(in) :: n_elements
        real, intent(in) :: x(n_elements), y(n_elements), error(n_elements)
        real, intent(in) :: parameters(2), limits(2, 2)
        real, intent(out) :: res
        real :: prior, likelihood
        
        call log_likelihood(parameters, x, y, error, n_elements, likelihood)
        call ln_uniform_prior_dist(parameters, limits, prior)
        
        res = prior + likelihood
        
    end subroutine ln_posterior_dist
        
end module

module MCMC

    contains
    
    subroutine Metropolis(data_block, parameters, limits, stepsize, &
                        n_iter, n_elements, chain_res, log_array, acceptance)
        
        use data_generator
        use probability
        
        integer, intent(in) :: n_elements, n_iter
        real, intent(in) :: data_block(3, n_elements), parameters(2), limits(2,2)
        real, intent(in) :: stepsize(2)
        real, intent(out) :: chain_res(n_iter, 2), log_array(n_iter), acceptance
        
        real :: x(n_elements), y(n_elements), y_err(n_elements)
        real :: dummy_parameters(2), new_parameters(2)
        integer :: n_accepted, iter
        real :: log_post, new_log_post

        call SRAND(time())    
        x = data_block(1, :)
        y = data_block(2, :)
        y_err = data_block(3, :)
        
        n_accepted = 0
        
        call ln_posterior_dist( &
                parameters, &
                 x,&
                 y, &
                 y_err, &
                 limits,&
                 n_elements, &
                 log_post)
        
        dummy_parameters = parameters
        do iter = 1,n_iter
            new_parameters = dummy_parameters + stepsize * gaussian_sampling(2)
            
            call ln_posterior_dist( &
                new_parameters, &
                 x,&
                 y, &
                 y_err, &
                 limits,&
                 n_elements, &
                 new_log_post)

            if (LOG(RAND()) < (new_log_post - log_post)) then
                dummy_parameters = new_parameters
                log_post = new_log_post
                n_accepted = n_accepted + 1
            end if
            chain_res(iter, :) = dummy_parameters
            log_array(iter) = log_post
        end do

        acceptance = real(n_accepted) / real(n_iter)
        
    end subroutine Metropolis

end module






















