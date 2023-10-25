function [ B, U, V ] = Implicitly_Shifted_Bidiag_QR( B, U, V )
    m = size( B, 1 ); 
    
    %B % debug

    % Introduce the bulge
    % Compute the first Givens' rotation
    % Use 11.2.1 where 
    % ( t( 1,1 ) - t( m,m ) 
    %       t( 2,1 )        )
    % We know the following...
    %   t(2,1) = B( 1,1 ) * B( 1,2 )
    %   
    tmm = ( B( m-1,m ) * B( m-1,m ) ) + ( B( m,m ) * B( m,m ) );
    t11 = B( 1,1 ) * B( 1,1 );
    t21 = B( 1,1 ) * B( 1,2 );
    G = Givens_rotation( [ t11 - tmm 
                           t21 ] ); 
    
    B_temp = [ B( 1,1 ) B( 1,2 ) 
               B( 2,1 ) B( 2,2 ) ] * G;

    B( 1:2, 1:2 ) = B_temp( 1:2, 1:2 );

    V_temp = V( :,1:2 ) * G;
    V( :, 1:2 ) = V_temp( :, 1:2 );

    %B % debug

    % Chase the bulge that will appear above super diagonal and below the
    % diagonal.  This only happens if the matrix is bigger than 2x2.
    % Otherwise get rid of one bulge.
    if m > 2
        for i = 1 : m - 2 
            % Bulge above the superdiagonal (horizontal partition of 2x3)
            G_hat = Givens_rotation( [ B( i,i )
                                       B( i+1,i ) ] );
    
            B_temp = G_hat' * [ B( i,i )   B( i,i+1 )   B( i,i+2 )
                                B( i+1,i ) B( i+1,i+1 ) B( i+1,i+2 ) ];
    
            % row 1 of horizontal partition
            B( i,i ) = B_temp( 1,1 ); B( i,i+1 ) = B_temp( 1,2 ); 
            B( i,i+2 ) = B_temp( 1,3 );
            % row 2 of horizontal partition
            B( i+1,i ) = B_temp( 2,1 ); B( i+1,i+1 ) = B_temp( 2,2 ); 
            B( i+1,i+2 ) = B_temp( 2,3 );

            U_temp = U( :, i:i+1 ) * G_hat;
            U( :, i:i+1 ) = U_temp( :, 1:2 );
    
            %B % debug
    
            % Bulge below the diagonal (vertical partition of 3x2)
            G = Givens_rotation( [ B( i,i+1 )
                                   B( i,i+2 ) ] );
    
            B_temp = [ B( i,i+1 )   B( i,i+2 ) 
                       B( i+1,i+1 ) B( i+1,i+2 )
                       B( i+2,i+1 ) B( i+2,i+2 ) ] * G;
            
            B( i,i+1 ) = B_temp( 1,1 ); B( i,i+2 ) = B_temp( 1,2 ); 
            B( i+1,i+1 ) = B_temp( 2,1 ); B( i+1,i+2 ) = B_temp( 2,2 ); 
            B( i+2,i+1 ) = B_temp( 3,1 ); B( i+2,i+2 ) = B_temp( 3,2 ); 

            V_temp = V( :, i+1:i+2 ) * G;
            V( :, i+1:i+2 ) = V_temp( :, 1:2 );
    
            %B % debug
        end
    end

    % Remove the bulge on the last row 
    G_hat = Givens_rotation( [ B( m-1,m-1 )
                               B( m,m-1 ) ] );

    B_temp = G_hat' * [ B( m-1,m-1 ) B( m-1,m )
                        B( m,m-1 )   B( m,m ) ];

    B( m-1:m, m-1:m ) = B_temp( 1:2, 1:2 );

    U_temp = U( :, m-1:m ) * G_hat; 
    U( :, m-1:m ) = U_temp( :, 1:2 );

    %B % debug
end