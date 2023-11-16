function [ varargout ] = polarplot_noneg( varargin )
%polarplot_noneg Make a polar plot with no negative values

ip = inputParser;
ip.addOptional('theta',[],@isnumeric);
ip.addRequired('rho',@isnumeric);
ip.KeepUnmatched = true;
ip.parse(varargin{:});

if(isempty(ip.Results.theta))
    rho = ip.Results.rho;
    theta =  linspace(0,2*pi,size(rho,1)).';
    otherArgs = varargin(2:end);
else
    theta = ip.Results.theta;
    rho = ip.Results.rho;
    otherArgs = varargin(3:end);
end

rho(rho < 0) = 0;

[varargout{1:nargout}] = polarplot(theta,rho,otherArgs{:});


end

