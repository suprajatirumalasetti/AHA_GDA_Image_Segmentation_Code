function [RI,GCE,VI,BDE,MHD,HD,JE,DICE,RFP,RFN] = segmentation_parameters(reference_image, seg_image)

sampleLabels1 = reference_image;
sampleLabels2 = seg_image;


if ndims(sampleLabels1) == 2
    sampleLabels1 = reshape(sampleLabels1, size(sampleLabels1,1), size(sampleLabels1,2), 1);
end
if ndims(sampleLabels2) == 2
    sampleLabels2 = reshape(sampleLabels2, size(sampleLabels2,1), size(sampleLabels2,2), 1);
end

channels = size(sampleLabels1, 3);

for i = 1:channels
    [ri(i), gce(i), vi(i)] = compare_segmentations(sampleLabels1(:,:,i), sampleLabels2(:,:,i));
    averageError(i) = compare_image_boundary_error(double(sampleLabels1(:,:,i)), double(sampleLabels2(:,:,i)));
    mhd(i) = ModHausdorffDist(double(sampleLabels1(:,:,i)), double(sampleLabels2(:,:,i)));
    hd(i) = HausdorffDist(double(sampleLabels1(:,:,i)), double(sampleLabels2(:,:,i)));
    [jaccardIdx(i), jaccardDist(i)] = jaccard_coefficient(double(sampleLabels1(:,:,i)), double(sampleLabels2(:,:,i)));
    [Jaccard(i), Dice(i), rfp(i), rfn(i)] = sevaluate(double(sampleLabels1(:,:,i)), double(sampleLabels2(:,:,i)));
end

RI = mean(ri);
GCE = mean(gce);
VI = mean(vi);
BDE = mean(averageError);
MHD = mean(mhd);
HD = mean(hd);
JE = mean(jaccardDist);
DICE = mean(Dice);
RFP = mean(rfp);
RFN = mean(rfn);

end
