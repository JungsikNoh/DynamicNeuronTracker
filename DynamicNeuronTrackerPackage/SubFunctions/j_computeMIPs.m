function [maxXY,maxZY,maxZX,three]=j_computeMIPs(vol,ZXRatio,minInt,maxInt,varargin)
% j_computeMIPs

ip = inputParser;
ip.addParameter('stripeSize',4,@isnumeric);
ip.addParameter('tracks',[],@(t) isa(t,'Tracks'));
ip.addParameter('frameIdx',[],@isnumeric);
ip.addParameter('usePrctile',false);
ip.addParameter('raw',false);
ip.addParameter('compThree',true);
ip.parse(varargin{:});
p=ip.Results;

% set other parameters
stripeSize = p.stripeSize; % the width of the stripes in the image that combines all three maximum intensity projections
stripeColor = 0; %the stripe color, a number between 0 (black) and 1 (white).  (If you're not using all of the bit depth then white might be much lower than 1, and black might be higher than 0.)

ScaledZ=ceil(size(vol,3)*ZXRatio);
% find the maximum intensity projections
maxXY = (max(vol, [], 3));
maxZY = imresize((squeeze(max(vol, [], 2))),[size(vol,1) ScaledZ]);
maxZX = imresize((squeeze(max(vol, [], 1))),[size(vol,2) ScaledZ]);

if(~p.raw)
	maxXY = uint8((2^8-1)*mat2gray(maxXY,double([minInt,maxInt])));
	maxZY = uint8((2^8-1)*mat2gray(maxZY,double([minInt,maxInt])));
	maxZX = uint8((2^8-1)*mat2gray(maxZX,double([minInt,maxInt])));
else
	maxXY=uint16(maxXY);
	maxZY=uint16(maxZY);
	maxZX=uint16(maxZX);
end

three=[];
% generate a single image with all three projections
if(p.compThree)
    if (~p.raw)
        three=uint8(zeros(size(maxXY,1)+stripeSize+size(maxZX',1),size(maxXY,2)+stripeSize+size(maxZY,2)));
    else
        three=uint16(zeros(size(maxXY,1)+stripeSize+size(maxZX',1),size(maxXY,2)+stripeSize+size(maxZY,2)));
    end
 
  three(1:size(maxXY,1),1:size(maxXY,2),:)=maxXY;
  three((size(maxXY,1)+stripeSize)+(1:size(maxZX',1)),1:size(maxZX',2),:)=maxZX';
  three((1:size(maxZY,1)),size(maxXY,2)+stripeSize+(1:size(maxZY,2)),:)=maxZY;
        
end
