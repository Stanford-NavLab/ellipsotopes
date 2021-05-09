function E_out = convhull(E_1,E_2)
% E = convhull(E_1,E_2)
%
% Compute the convex hull of two ellipsotopes
%
% Authors: Adam Dai and Shreyas Kousik
% Created: 21 Apr 2021
% Updated: 5 May 2021 (messing with I_extra)

% extract properties
[p_1,c_1,G_1,A_1,b_1,I_1] = E_1.get_properties ;
[p_2,c_2,G_2,A_2,b_2,I_2] = E_2.get_properties ;

% sanity check
if p_1 ~= p_2
    error(['We do not yet support operations on ellipsotopes with ',...
        'different p-norms!'])
end
p = p_1;

% retrieve sizes of things
n = size(c_1,1);
m_1 = size(G_1,2);
m_2 = size(G_2,2);
q_1 = size(A_1,1);
q_2 = size(A_2,1);

% create convex hull ellipsotope properties
c_CH = 0.5*(c_1+c_2);
G_CH = [G_1, G_2, 0.5*(c_1-c_2), zeros(n,2*(m_1+m_2))];

A_31 = [eye(m_1); -eye(m_1); zeros(2*m_2,m_1)];

A_32 = [zeros(2*m_1,m_2) ; eye(m_2); -eye(m_2)];

A_30 = [-0.5*ones(2*m_1,1); 0.5*ones(2*m_2,1)];

A_CH = [A_1,            zeros(q_1,m_2), -b_1/2, zeros(q_1,2*(m_1+m_2));
        zeros(q_2,m_1), A_2,             b_2/2, zeros(q_2,2*(m_1+m_2));
        A_31,           A_32,            A_30,  eye(2*(m_1+m_2))];

b_CH = [b_1/2; b_2/2; -0.5*ones(2*(m_1+m_2),1)] ;

m_3 = m_1 + m_2 ;
J_extra = (m_3 + 2):(3*m_3 + 1) ;
I_extra = num2cell(J_extra) ;

% I_11 = I_1 ;
% I_12 = shift_index_set(I_1,max(cell2mat(I_11))) ;
% I_21 = shift_index_set(I_2,max(cell2mat(I_12))) ;
% I_22 = shift_index_set(I_2,max(cell2mat(I_21))) ;
% 
% I_extra = shift_index_set([I_11,I_12,I_21,I_22],m_3+1) ;

I_2_shifted = shift_index_set(I_2,m_1) ;

I_CH = [I_1,I_2_shifted,{m_3+1},I_extra] ;

% create output
E_out = ellipsotope(p,c_CH,G_CH,A_CH,b_CH,I_CH);

end