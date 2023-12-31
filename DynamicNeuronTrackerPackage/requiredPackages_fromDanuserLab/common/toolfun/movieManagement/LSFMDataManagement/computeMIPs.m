function [maxXY,maxZY,maxZX,three]=computeMIPs(vol,ZXRatio,minInt,maxInt,varargin)
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

% if(p.usePrctile)
%     
% end
if(~p.raw)
	maxXY = uint8((2^8-1)*mat2gray(maxXY,double([minInt,maxInt])));
	maxZY = uint8((2^8-1)*mat2gray(maxZY,double([minInt,maxInt])));
	maxZX = uint8((2^8-1)*mat2gray(maxZX,double([minInt,maxInt])));
else
	maxXY=uint16(maxXY);
	maxZY=uint16(maxZY);
	maxZX=uint16(maxZX);
end
% F=figure();
% imshow((maxXY));
% F=figure();
% hist(double(maxXY(:)));
% waitfor(F);

three=[];
% generate a single image with all three projections
if(p.compThree)
	three=projMontage(maxXY,maxZX',maxZY',false);
end

% three = uint8((2^8-1)*mat2gray(three,double([minInt,maxInt])));
