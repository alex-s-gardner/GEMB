function D = fast_divisors(N)

    % Find divisors up to the square root
    K = 1:ceil(sqrt(N));
    D = K(rem(N,K)==0);

    % Find corresponding divisors > sqrt(N) and combine
    D = [D sort(N./D)];
end