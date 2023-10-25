function [ S, U, V ] = SVD_BiDiag_ImpShift( B )
    maxiter = 500000;
    m = size( B, 1 );  
    U = eye( m );
    V = eye( m );
    
    % Add deflation when honing in on the singular values of B.
    for i = 0: abs( 2 - m )
        iter = 0; % debug

        % current m (cm) iteration
        cm = m - i;

        while Check_Run_Criteria( B( cm-1, cm ), B( cm-1, cm-1 ), ...
                B( cm, cm ) ) && iter < maxiter
            [ B( 1:cm, 1:cm ), U( :, 1:cm ), ...
                V( :, 1 : cm ) ] = ...
                Implicitly_Shifted_Bidiag_QR( B( 1:cm, 1:cm ), ...
                    U( :, 1:cm ), V ( :, 1:cm) );
            iter = iter + 1; % debug
        end

        %fprintf("cm: %d, iteration(s): %d \n", cm, iter); % debug
    end

    %pre_sorted_B = B % debug
    [ B, U, V ] = Sort_SVD_B( B, U, V );
    S = diag( B );
end