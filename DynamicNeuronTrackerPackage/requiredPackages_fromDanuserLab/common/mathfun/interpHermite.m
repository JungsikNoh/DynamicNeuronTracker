function [ varargout ] = interpHermite( x, y, varargin )
%interpHermite Interpolate with known derivative values
%
% INPUT
% x - abscissa where values f(x) are known
% y - values of f(x)
% s - values of f'(x)
% c - values of f''(x)
% ... values of f^(n)(x)
% xq - unknown values, f(xq), wanted (not specified for sp or pp output)
%
% PARAMETERS
% 'output' - One of 'values' (default), 'sp', or 'pp'
%            values: Values at xq
%            sp: spline output (xq not expected)
%            pp: piecewise polynomial (xq not expected)
% 'knots'  - Default: Use augknt with breaks at the given x, see spapi
%                     'optknt' will use optknt to alculate knots instead
%
% OUTPUT
% depends on output parameter
%  first output represents f(x)
% second output represents f'(x)
% etc...
%
% EXAMPLE
% x = (0:5)/5*2*pi;
% Y = vertcat(sin(5*x),5*cos(5*x),-25*sin(5*x),-125*cos(5*x));
% % Plot using default
% interpHermite(x,Y(1,:),Y(2,:),Y(3,:),Y(4,:),'output','sp')
% subplot(4,1,1)
% title('Plotting output of interpHermite with nargout = 0');
%
% figure;
% xx = (0:360)/360*2*pi;
% %Obtain values and derivatives
% [vq,vqd,vqdd,vqddd] = interpHermite(x,Y(1,:),Y(2,:),Y(3,:),Y(4,:),xx)
% plot(xx,sin(5*xx),'-');
% hold on;
% plot(xx,vq,'o');
% grid on;
% legend({'sin(5*xx)','vq = interpHermite'})
% title('Interpolation at xx')
%
% % Compare to pwch
% x = (0:20)/20*2*pi;
% Y = vertcat(sin(5*x),5*cos(5*x));
% ppInterpHermite = interpHermite(x,Y(1,:),Y(2,:),'output','pp')
% ppPWCH = pwch(x,Y(1,:),Y(2,:));
% figure;
% plot(xx,sin(5*xx),'k','LineWidth',5);
% hold on; fnplt(ppPWCH);
% fnplt(ppInterpHermite,'--');
% legend({'sin(5*xx)','pwch','interpHermite'});
% title('interpHermite vs pwch (First derivative only)');
% grid on
%
% ppInterpHermite.coefs - ppPWCH.coefs
% norm(ppInterpHermite.coefs - ppPWCH.coefs)
%
%
% See also pchip, pwch, spapi, augknt, optknt

% Mark Kittisopikul, February 2018
% Goldman Lab
% Northwestern University

%% Process input
nderiv = 0;
for nderiv=1:length(varargin)
    if(ischar(varargin{nderiv}))
        nderiv = nderiv - 1;
        break;
    end
end
if(nderiv < length(varargin))
    ip = inputParser;
    ip.addParameter('output','values',@(x) any(x(1) == 'vsp'));
    ip.addParameter('knots',[]);
    ip.parse(varargin{nderiv+1:end});
    params = ip.Results;
else
    params.output = 'values';
    params.knots = [];
end

switch(params.output(1))
    case 's'
    case 'p'
    otherwise
        xq = varargin{nderiv};
        nderiv = nderiv - 1;
end

%% Specify x nderiv times, (Could have used brk2knt)
% Ensure x is a row vector
x = x(:).';
x = repmat(x,nderiv+1,1);
x = x(:).';

%% Ensure y is a row vector
y = y(:).';

%% Pack known values and derivatives into Y
derivs = cell(1,nderiv);
for i=1:nderiv
    derivs{i} = varargin{i}(:).';
end
Y = vertcat(y,derivs{:});
Y = Y(:).';

%% Setup knots using augknt, consider optknt
if(isempty(params.knots))
    % This is the true hermite interpolation with breaks at the x specified
    knots = augknt(x,nderiv*2+2);
elseif(params.knots(1) == 'o')
    % This is 'optimal'
    knots = optknt(x,nderiv*2+2);
else
    knots = params.knots;
end

%% The core of the function, use spapi
sp = spapi(knots,x,Y);

%% Plot if no output requested
if(nargout == 0)
    % Plot a figure showing the function and it's derivatives
    h = figure;
    m = nderiv+1;
    s = 1:m:length(x);
    if(exist('xq','var'))
        f = @(sp) plot(xq,spval(sp,xq));
    else
        f = @fnplt;
    end
    xx = x(s);
    
    subplot(m,1,1);
    f(sp);
    hold on;
    plot(xx,y,'o');
    ylabel('f(x)');
    grid on;
    
    for d=1:nderiv
        subplot(m,1,d+1);
        sp = fnder(sp);
        f(sp);
        hold on;
        plot(xx,derivs{d},'o');
        ylabel(['f^{(' num2str(d) ')} (x)']);
        grid on;
    end
end

%% Actual output
varargout = cell(1,nargout);
switch(params.output(1))
    case 's'
        % spline output
        varargout{1} = sp;
        for i=2:nargout
            varargout{i} = fnder(varargout{i-1});
        end
    case 'p'
        % piecewise polynomial output
        pp = sp2pp(sp);
        varargout{1} = pp;
        for i=2:nargout
            varargout{i} = fnder(varargout{i-1});
        end
    otherwise
        % value output
        for i=1:nargout
            varargout{i} = spval(sp,xq);
            sp = fnder(sp);
        end
end

end