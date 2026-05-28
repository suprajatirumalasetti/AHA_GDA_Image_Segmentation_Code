
function [jaccardIdx,jaccardDist] = jaccard_coefficient(img_Orig,img_Seg)

inter_image = img_Orig & img_Seg;
union_image = img_Orig | img_Seg;
jaccardIdx = sum(inter_image(:))/sum(union_image(:));
jaccardDist = 1 - jaccardIdx;