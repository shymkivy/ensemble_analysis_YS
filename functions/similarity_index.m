function [SI_cosine, SI_hammilarity] = similarity_index(vec1, vec2)
% calculates Similarity Index
% we use this instead cosine of angle because this accounts for magnitude
% and cosine doesn't
% SI = (Ca * Cb) / ((||Ca||^2 + ||Cb||^2)/2) - accounts for magnitude diff
% best for highly skewed distributions
% cos(ang) = (Ca * Cb) / (||Ca|| * ||Cb||) - ok for binary
    dot_prod = vec1*vec2';
    
%     vec1_norm = vec1./vecnorm(vec1')';
%     vec2_norm = vec2./vecnorm(vec2')';
%     si1 = vec1_norm*vec1_norm';
%     
%     magnitudes1 = sqrt(sum(vec1.*vec1,2));
%     magnitudes2 = sqrt(sum(vec2.*vec2,2));
    
    magnitudes1 = vecnorm(vec1')';
    magnitudes2 = vecnorm(vec2')';
    
    %SI_hammilarity = dot_prod./((magnitudes1.^2 + magnitudes2.^2)/2);
    
    SI_hammilarity = dot_prod./((magnitudes1.^2 + magnitudes2.^2')/2);
    SI_cosine =  dot_prod./(magnitudes1 .* magnitudes2');
    
    
    
    
end