function [Xf]=nspds(A)



% first, check if the matrix is SPD by attempting to take the choleski factor
[~,test1]=chol(A,'lower');


% if it is SPD in the first place (test value=0), no need to apply any change
if test1<1
Xf=A;


% if it is not SPD, fix it
elseif test1>=1

% obtain the matrix Xf by applying theorem 2.1

% as a prelimineary, correct symmetry
C=tril(A)+tril(A,-1)';

% create B, as defined in theorem 2.1
B=(C+C')/2;
% produce the polar decomposition of B
% first, implement a singular value decomposition of B
[U,S,V]=svds(B);
% then recover the polar matrix H from the singular value decomposition formula: H=V*S*V'
H=V*S*V';
% finally, recover Xf from the formula in theorem 2.1: Xf=(B+H)/2
Xf=(B+H)/2;
% ensure symmetry

Xf=(Xf+Xf')/2;

% now the matrix should be symmetric positive definite: test again
[~,test2]=chol(Xf,'lower');

   % if the test is passed and Choleski factorisation was possible (test2=0), then Xf needs no further work;
   if test2<1
   
   % if the test is not passed, it must be because of some very small numerical disturbance
   % hence keep modifying Xf very slightly until it becomes SPD
   elseif test2>=1
      n=1;
      while test2>=1
      mineig=min(eig(Xf))+0.0000000000000001;
      Xf=Xf+(-mineig*n.^2+eps(mineig))*eye(size(Xf)); 
      % test again
      [~,test2]=chol(Xf,'lower');
      n=n+1;
      end

   end

end



























