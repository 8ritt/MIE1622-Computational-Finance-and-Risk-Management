function gval = computeGradERC (x)

global QQQ

n = size(QQQ,1) ;

if(size(x,1)==1)
    x = x';
end
y = x.*(QQQ*x);
% Insert your gradiant computations here
% You can use finite differences to check the gradient

gval = zeros(n,1);
 for i = 1:n
     for j = 1:n
         diff1 = QQQ(i,:)*x + QQQ(i,i) * x(i);
         diff2 = QQQ(i,j) * x(i);
         g = (y(i) - y(j)) * (diff1 -diff2);
         gval(i,:) = gval(i,:) + g;
     end
     gval(i,:) =4 * gval(i,:);
 end
% gval = gradient(x);
end
