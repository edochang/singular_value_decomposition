function [ S, U, V ] = SVD_ALAFF( A )
    m = size( A, 1 );

    % Create vector in which to store the scalars tau from the Householder
    % transformations
    t = rand( m, 1 );

    % Create vector in which to store the scalars rho from the Householder
    % transformations
    r = rand( m, 1 );

    % Compute the reduction to tridiagonal form
    [ B_With_HouseV, t, r ] = BiRed( A, t, r );

    % Extract the bidiagonal matrix from B.
    B = BiFromB( B_With_HouseV );

    % Recommendation from Piazza Post @417 to use FormQ( A, t ) to create
    % the initial U and V.  t is m and r is m-1.
    % Note: Per TA we can pad t or r with 1's to generally FormQ.  
    % Assuming he means using the Identity per 3.3.5. 
    U = FormQ( B_With_HouseV, t );

    % See 3.3.5 in understanding the theory of Forming Q.  Given the
    % horizontal Householder Vectors only goes up to 1 to m-1, when we
    % transpose B to put it in the form that FormQ can use, we are missing
    % a householder vector to the left of the bidiagonal values and as a
    % result is given a m-1 x m-1 V from FormQ().
    % For a m x m V, start with an identity per 3.3.5. knowing that there's
    % no Horizontal Hoursholder vector for the first column and row of this
    % matrix to update.
    V = eye( m );
    V_FormQ = FormQ( B_With_HouseV( 1:m-1, 2:m )', r );
    V( 2:m, 2:m ) = V_FormQ;

    %display_UBiVt = U * B * V' % debug
    %L2_norm_error = norm( U * B * V' - A, 2 ) % debug

    [ S, U, V ] = Decompose_B_with_UA_VA( B, U, V );
end

function [ B, U, V ] = Decompose_B_with_UA_VA( B, U, V )
maxiter = 500000;
    m = size( B, 1 );  
    
    % Add deflation when honing in on the singular values of B.
    for i = 0: abs( 2 - m )
        iter = 0; % debug

        % current m (cm) iteration
        cm = m - i;

        while Check_Run_Criteria( B( cm-1,cm ), B( cm-1,cm-1 ), ...
                B( cm,cm ) ) && iter < maxiter
            [ B( 1 : cm, 1 : cm ), U( :, 1 : cm ), ...
                V( :, 1 : cm ) ] = ...
                Implicitly_Shifted_Bidiag_QR( B( 1 : cm, 1 : cm ), ...
                    U( :, 1 : cm ), V ( :, 1 : cm) );
            iter = iter + 1; % debug
        end

        %fprintf("cm: %d, iteration(s): %d \n", cm, iter); % debug
    end

    [ B, U, V ] = Sort_SVD_B( B, U, V );
    B = diag( B );
end
