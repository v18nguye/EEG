function s_h = SISSY(x,A,T,lambda,alpha)

% Estimate the source signal by using SISSY
rho = 1;
iter = 60;

P=sparse(rho*(T.'*T+speye(size(T,2))));
APi=A/P;

L=chol(eye(size(A,1))+APi*A.','lower');
s1=A.'*x;

% Intialize parameters u,v,y,z,s
z=zeros(size(T,1),size(x,2));
u=zeros(size(T,1),size(x,2));
s=zeros(size(T,2),size(x,2));
v=zeros(size(T,2),size(x,2));
y=zeros(size(T,2),size(x,2));

for i =1:iter
    
    % update s
    b=s1+rho*(T.'*(z+u/rho)+y+v/rho);
    s=P\b-APi.'*(L.'\(L\(APi*b)));
    
    % update z
    z_ = T*s - 1/rho*u;
    mask1 = (z_ > lambda/rho);
    z(mask1) = z_(mask1) - lambda/rho;
    mask2 = (z_ < -lambda/rho);
    z(mask2) = z_(mask2) + lambda/rho;
    mask3 =(-lambda/rho<= z_ <=lambda/rho);
    z(mask3) = 0;
    
    % update y
    y_ = s - 1/rho*v; 
    mask1 = y_ > lambda*alpha/rho;
    y(mask1) = y_(mask1) - lambda*alpha/rho;
    mask2 = y_ < -lambda*alpha/rho;
    y(mask2) = y_(mask2) + lambda*alpha/rho;
    mask3 =(-lambda*alpha/rho<= y_ <=lambda*alpha/rho);
    y(mask3) = 0;
    
    % update u
    u = u + rho*(z-T*s);
    
    % update v
    v = v + rho*(y-s);
   
end

s_h = s;
end