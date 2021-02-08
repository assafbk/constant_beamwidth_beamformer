%
% the function takes the [weights x freqs] matrix of the produced array's vector
% and transforms it to the equivalent [weights x freqs] matrix of the physical
% array (i.e. the equivalent filter should make the physical beampattern 
% equal the produced beampattern)
%
% w_prod - the produced array's weigths
% M1 - 1st virtual array's size
% M2 - 2nd virtual array's size
%
function w_phy = produced_array_weights_to_equivalent_physical_array_weights(w_prod, M1, M2)

M = M1+M2-1; % physical array's size
num_of_freqs = size(w_prod);    
num_of_freqs = num_of_freqs(2);
w_phy = zeros([M num_of_freqs]);

for p=0:(M2-1)
    for q=1:M1
        w_phy(p+q,:)= w_phy(p+q,:) + w_prod(q+M1*p,:);
    end
end

end

