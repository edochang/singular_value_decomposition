function [ A_out, t_out, r_out ] = BiRed( A, t, r )

  [ ATL, ATR, ...
    ABL, ABR ] = FLA_Part_2x2( A, ...
                               0, 0, 'FLA_TL' );

  [ tT, ...
    tB ] = FLA_Part_2x1( t, ...
                         0, 'FLA_TOP' );

  [ rT, ...
    rB ] = FLA_Part_2x1( r, ...
                         0, 'FLA_TOP' );
                     
  while ( size( ATL, 1 ) < size( A, 1 ) )

    [ A00,  a01,     A02,  ...
      a10t, alpha11, a12t, ...
      A20,  a21,     A22 ] = FLA_Repart_2x2_to_3x3( ATL, ATR, ...
                                                    ABL, ABR, ...
                                                    1, 1, 'FLA_BR' );

    [ t0, ...
      tau1, ...
      t2 ] = FLA_Repart_2x1_to_3x1( tT, ...
                                    tB, ...
                                    1, 'FLA_BOTTOM' );

    [ r0, ...
      rho1, ...
      r2 ] = FLA_Repart_2x1_to_3x1( rT, ...
                                    rB, ...
                                    1, 'FLA_BOTTOM' );
                                
    %------------------------------------------------------------%

    tau1 = 1;
    %if size( a21, 1 ) ~= 0 % Note: Per HQR, test_HQR, and test_FormQ, 
    % this is not needed.  We do want the .5 tau even though it produces 
    % a 0 Householder vector and carry through ||x||_2 to the last element.
        
    %a12t % debug
    %A22 % debug
    
    % Create the 0's on the first column
    % Calculates the Householder Vector u and tau
    [ alpha11, ...
      a21, tau1 ] = Housev( alpha11, ...
                              a21 );
    u = [ 1
          a21 ];
    A_sub = [ a12t
              A22 ]; 

    % Calculate the Householder Transformation: I - (1/tau)uu^H
    H = A_sub - (  u * ( u' ) * A_sub / tau1 );
    %disp("New a12t") % debug
    a12t = H( 1, : );
    %disp("New A22") % debug
    A22 = H( 2 : end, : );
    
    %tau1 % debug

    rho1 = 0;
    % Create the 0's in the first row of this updated matrix
    %if size( a12t, 2 ) > 1 % Note: Same adjustments must be made here
    % given the note above.  Thus, t is m and r is m-1.
    if size( a12t, 2 ) > 0
        [ u12, rho1 ] = Housev1(a12t');
        a12t = u12';
    
        % Update A22
        u2 = [ 1
               u12( 2:end, 1 ) ];
        %disp("New A22") % debug
        A22 = A22 - ( A22 * u2 * ( u2' ) / rho1 );
    end

    %rho1 % debug
    
    %------------------------------------------------------------%

    [ ATL, ATR, ...
      ABL, ABR ] = FLA_Cont_with_3x3_to_2x2( A00,  a01,     A02,  ...
                                             a10t, alpha11, a12t, ...
                                             A20,  a21,     A22, ...
                                             'FLA_TL' );

    [ tT, ...
      tB ] = FLA_Cont_with_3x1_to_2x1( t0, ...
                                       tau1, ...
                                       t2, ...
                                       'FLA_TOP' );

    [ rT, ...
      rB ] = FLA_Cont_with_3x1_to_2x1( r0, ...
                                       rho1, ...
                                       r2, ...
                                       'FLA_TOP' );
                                   
  end

  A_out = [ ATL, ATR
            ABL, ABR ];

  t_out = [ tT
            tB ];
        
  r_out = [ rT
            rB ];

return