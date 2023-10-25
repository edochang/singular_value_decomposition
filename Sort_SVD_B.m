function [ B, U, V ] = Sort_SVD_B( B, U, V )
    %disp( "Called Sort_SVD_B()" ) % debug
    %B_pre_absolute_value_and_permutation = B % debug
    abs_B_diag = abs( diag(B) );

    % Sort by column in descending order for B
    [ ~, Ir ] = sort( abs_B_diag, 'descend' );

    % Instantiate with an Identity matrix of size m x m.
    P = eye( size( abs_B_diag, 1 ) );
    % Create the left permutation matrix that permutes the rows.
    P = P( Ir, : );

    % Per Chapter 11.2.1, to permute the spectral decomposition
    % (QP')(PDP')(QP')'
    % Permute U
    U = U * P';

    % Permute B
    B = P * B * P';

    % Permute V
    % Note: Double check if the transpose is accounted for here or in the
    % SVD calculation further down in the process.  (QP')' vs. (QP')
    V = V * P';
end