function dle_value=DLE(idx,idx_est,r_grid)

 
r_k = r_grid(idx,:); %r_k is the  positions of dipoles.
r_l = r_grid(idx_est,:); %r_l is the positions of dipoles.


L = length(r_k);
L_est = length(r_l);

DLE1 = 0;
for i = 1:L %The indices of dipoles of original sources 
    
    suqreDifference_kl=(repmat( r_k(i,:) , L_est , 1 ) - r_l).^2;
    DLE1 = DLE1 + min(sqrt(sum(suqreDifference_kl, 2)));  % only add the minimum norm value   
end

DLE2 = 0;
for j = 1:L_est %The indices of dipoles of estimated sources 
    
    
    suqreDifference_lk=(repmat( r_l(j,:) , L, 1 ) - r_k).^2;
    DLE2 = DLE2 + min(sqrt(sum(suqreDifference_lk, 2))); 
end

dle_value = 1/(2*L) * DLE1 + 1/(2*L_est) * DLE2;

end

