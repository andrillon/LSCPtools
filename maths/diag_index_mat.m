function idx=diag_index_mat(M)

% from https://math.stackexchange.com/questions/1392491/measure-of-how-much-diagonal-a-matrix-is

j=ones(1,size(M,1));
r=1:size(M,1);
r2=r.^2;

n=j*M*j';
Sx=r*M*j';
Sy=j*M*r';
Sx2=r2*M*j';
Sy2=j*M*r2';
Sxy=r*M*r';

idx=(n*Sxy-Sx*Sy)/(sqrt(n*Sx2-(Sx)^2)*sqrt(n*Sy2-(Sy)^2));