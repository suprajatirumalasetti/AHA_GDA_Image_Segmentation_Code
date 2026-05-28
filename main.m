
clc;
clear all;
close all;
  f1=dir('D:\BSDS500-master\BSDS500-master\BSDS500\image segmentation data set\images\BSDS AHA 2 IMAGES\*.jpg');
 fil={f1.name};
 loc={f1.folder};
for sri=1:numel(fil)
for SRI=2:8
     file=fil{sri};
     loca=loc{sri};
     location=fullfile(loca,file);
 P = imread(location);
 [rows, cols, ~] = size(P);                
     [filepath,in_name,in_ext] = fileparts(file);
      destinationFolder = 'D:\VIT\IMAGE SEGMENTATION\AHA-2\AHA PAPER 2\AHA2\FINAL AHA 2\scientific reports revision\final\scientific reprots\MINOR REVISION\AHA2 CODE\Hybrid Artificial Hummingbird code\OUTPUT';
         if ~exist(destinationFolder, 'dir')
         mkdir(destinationFolder);
     end
       destinationFolder1 = 'D:\VIT\IMAGE SEGMENTATION\AHA-2\AHA PAPER 2\AHA2\FINAL AHA 2\scientific reports revision\final\scientific reprots\MINOR REVISION\AHA2 CODE\Hybrid Artificial Hummingbird code\OUTPUT';
          if ~exist(destinationFolder1, 'dir')
         mkdir(destinationFolder1);
      end
 SRI1 = sprintf( '%s',num2str(SRI) );
     fullDestinationFileName = fullfile(destinationFolder,sprintf('%s_TWDOfusion_%s.%s', in_name,SRI1, in_ext));
    fullDestinationFileName1 = fullfile(destinationFolder1,sprintf('%s_TWDOwithoutfusion_%s.%s', in_name,SRI1, in_ext));
      filename = 'D:\VIT\IMAGE SEGMENTATION\AHA-2\AHA PAPER 2\AHA2\FINAL AHA 2\scientific reports revision\final\scientific reprots\MINOR REVISION\AHA2 CODE\Hybrid Artificial Hummingbird code\OUTPUT\AHA IMAGES.xlsx';
for ix=1:3
I=P(:,:,ix);
     Dim=SRI;
MaxIteration=2000;
PopSize=30;
 FunIndex=1;
  [m n o]=size(I);
  sgrays = double(I);
 caddterm = m*n*8;
  for q=1:o
 for i=1:256
      bp = sgrays(:,:,q)>i-1;      
      negmat = bp-1;      
      b = negmat+bp;      
 
   b(1,:) = 0;
     b(end,:) = 0;
      b(:,1) = 0;
      b(:,end) = 0;
 
      pienergy = 0;          
      for j=2:m-1
          for k=2:n-1
              pienergy = pienergy + (b(j,k).^2) - sum(sum(b(j,k).*b(j-1:j+1,k-1:k+1)));
          end
      end
      E(q,i) = pienergy + caddterm;
  end
  end
   I_energy = E';
  tic
 [BestX,BestF,HisBestF]=AHA_GDA(FunIndex,MaxIteration,PopSize,I_energy,Dim);
thresh = BestX;
segout_ = multilevel_segment(I, BestX); 
  segout_ = uint8(255 * mat2gray(segout_));     
     segout_ = imresize(segout_, [rows cols]);     
    segout_RGB(:,:,ix) = segout_;                  
    threshold_values(:,ix)=thresh;
    output_image_name = sprintf('%s_threshold_%d_RGB_%d.jpg', in_name, SRI, ix);
    imwrite(segout_, fullfile(destinationFolder, output_image_name));
 hFig = figure;
 probI = I_energy / (sum(I_energy) + eps);
 intensity = BestX;
  plot(probI);
 vmax = max(probI);          
for i = 1:length(BestX)
   line([intensity(i), intensity(i)], [0 vmax], 'Color', 'r','LineWidth', 2,'Marker', '.', 'LineStyle', '-');
end
title('Histogram');
xlabel('Intensity');
ylabel('Probability');
xlim([0 255])
set(gca, 'LineWidth', 2);
set(gca, 'Box', 'on');  
             drawnow; 
                output_image_names = sprintf('%s_histogram_%d_round_%d.jpg', in_name, SRI, ix);
                saveas(hFig, fullfile(destinationFolder, output_image_names));
                close (hFig)
               TIME(ix)=toc;
end
T = segout_RGB;  
final_output_name = sprintf('%s_Final_RGB_threshold_%d.jpg', in_name, SRI);
fullDestinationFileName = fullfile(destinationFolder, final_output_name);
imwrite(uint8(T), fullDestinationFileName);  
seg_image1=uint8(T);
T=uint8(T);
for e= 1:3
A=P(:,:,e);
B=T(:,:,e);
 D = abs(uint8(A) - uint8(B)).^2;
 mse(e) = sum(D(:))/numel(A);
 psnr(e) = 10*log10(255*255/mse(e));
K = [0.6 0.6];
 [mssim, ssim_map] = ssim_index(A,B,K);
 MS(e)=mssim;
 [FSIM, FSIMc] = FeatureSIM(A, B);
 FS(e)=FSIM;
muA = mean2(A);
muB = mean2(B);
sigmaA = std2(A);
sigmaB = std2(B);
sigmaDiff = std2(A(:) - B(:)); 
term1 = (2 * muA * muB) / (muA^2 + muB^2 + eps);
term2 = (2 * sigmaA * sigmaB) / (sigmaA^2 + sigmaB^2 + eps);
term3 = 1 - (sigmaDiff / (sigmaA * sigmaB + eps));
QILV(e) = term1 * term2 * term3;
CORR(e) = corr2(A, B);
    GI = edge(uint8(A), 'sobel');
    GS = edge(uint8(B), 'sobel');
    numerator = sum(sum(GI .* GS));
    denominator = sqrt(sum(GI(:).^2) * sum(GS(:).^2)) + eps;
    EPI(e) = numerator / denominator;
    H = joint_histogram(uint8(A), uint8(B));
    MI(e) = mutual_information(H);
end
FS;
MS;
pixel = imhist(P);
Entropy=entropy(uint8(T));
ME = 1-(nnz(~P & ~T) + nnz(P & T))/sum(pixel);
MSE=sum(mse)/length(mse);
PSNR=sum(psnr)/length(psnr);
SSIM=sum(MS)/length(MS);
FSIM=sum(FS)/length(FS);
QILVf = mean(QILV);
CORRf = mean(CORR);
EPIf  = mean(EPI);
MIf   = mean(MI);
[RI,GCE,VI,BDE,MHD,HD,JE,DICE,RFP,RFN]=segmentation_parameters(P,T);
thrrrr=threshold_values(:)';
TT = table(cellstr(in_name),RI,GCE,VI,BDE,MHD,HD,JE,DICE,RFP,RFN,ME,SSIM,FSIM,PSNR,QILVf,CORRf,EPIf,MIf,TIME,BestF,thrrrr);
    TT(1,:)
    sheet_name = sprintf( '%s',num2str(SRI) );
    xlRange= sprintf( 'A%s',num2str(sri) ); 
    writetable(TT,filename,'Sheet',sheet_name,'WriteVariableNames',false,'Range',xlRange);
    clear threshold_values

end

end