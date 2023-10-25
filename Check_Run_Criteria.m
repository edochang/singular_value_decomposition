% From 10.3.6, when abs( m-1, m ) is small relative to ( m - 1, 
% m - 1 ) and ( m, m ): abs( T1( m-1, m ) ) <= eps * sqrt( abs( 
% T1( m-1, m-1 ) ) + abs( T1( m, m ) ) ).  For the while loop 
% to keep running, invert the inquality to be greater than (>).
function [ criteria_not_met ] = Check_Run_Criteria( b_1m, b_m1m1, b_mm )
    criteria_not_met = abs( b_1m ) > eps * sqrt( abs( b_m1m1 ) + ...
        abs( b_mm ) );
end