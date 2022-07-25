syms A B C D
syms a b c d n1 n2

Fcurly = [a*A-d*B, b*C+d*A-d*C, C,c*D]
Vcurly = [A-b*A, b*B+c*B, c*C, D-a*D]

FF=jacobian(Fcurly, [A B C D]);
VV=jacobian(Vcurly, [A B C D]);
FF
VV



MatrixF=subs(FF, [A B C D], [0, 0, 0, 0])
MatrixV=subs(VV, [A B C D], [0, 0, 0, 0])

MatrixV = subs(MatrixV, [b+c 1-a], [n1 n2]);
MatrixV

RR=-MatrixF*inv(MatrixV)

% Find eigenvalue of RR, largest eigenvalue=R0
syms lambda
pp=det(RR-lambda*eye(4));
pp_factors=factor(pp);
num_factors=length(pp_factors);
eigen_values=[];
for i=1:num_factors
eigen_values=[eigen_values; solve(pp_factors(i)==0, lambda,'MaxDegree', 4)];
end

a = 0.5
b = 0.3
c = 0.4
d = 0.6

n1 = b+c
n2 = 1-a

%eivenvalues are copied from  eig=solve(p, lambda)
sol1=double(subs(eigen_values(1)));
sol2=double(subs(eigen_values(2)));
sol3=double(subs(eigen_values(3)));
sol4=double(subs(eigen_values(4)));

sol=[sol1 sol2 sol3 sol4]

% %%%% largest eigenvalue is R0
[subR0,max_index]=max(sol)
R0=eigen_values(max_index)


syms n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 j1 j2 j3 j4 j5 j6 j7 j8

