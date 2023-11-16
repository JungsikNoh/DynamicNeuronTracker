classdef DetectPolesProcess < DetectionProcess
%PointSourceDetectionProcess3D is a concrete class of a point source
%detection process for 3d
    
    methods (Access = public)
        function obj = DetectPolesProcess(owner, outputDir, funParams)
            % Constructor of the SubResolutionProcess
            super_args{1} = owner;
            super_args{2} = DetectPolesProcess.getName;
            super_args{3} = @detectPolesProcess;
            
            if nargin < 3 || isempty(funParams)  % Default funParams
                if nargin <2, outputDir = owner.outputDirectory_; end
                funParams = DetectPolesProcess.getDefaultParams(owner,outputDir);
            end
            
            super_args{4} = funParams;
            
            obj = obj@DetectionProcess(super_args{:});
            obj.is3Dcompatible_ = true;

        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'detectPolesZ', 'detectPolesMIP'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
            ip.addParameter('useCache',false,@islogical);
            ip.addParameter('iZ',[], @(x) ismember(x,1:obj.owner_.zSize_)); 
            ip.addParameter('output',outputList{1}, @(x) all(ismember(x,outputList)));
            ip.addParameter('projectionAxis3D','Z', @(x) ismember(x,{'Z','X','Y','three'}));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            projAxis3D = ip.Results.projectionAxis3D;
            iZ = ip.Results.iZ;
            if ischar(output),output={output}; end
            
            if ~obj.funParams_.isoOutput
                ZXRatio = obj.owner_.pixelSizeZ_/obj.owner_.pixelSize_;
                ZXRatio_z = 1;
            else
                ZXRatio = 1;
                ZXRatio_z = obj.owner_.pixelSizeZ_/obj.owner_.pixelSize_;
            end
%             XRatio=size(img,2)/(XLimit(2)-XLimit(1));
%             YRatio=size(img,1)/(YLimit(2)-YLimit(1));
            
            for iout = 1:numel(output)
                switch output{iout}
                    case {'detectPolesZ'}
                        s = cached.load(obj.outFilePaths_{1,iChan}, '-useCache', ip.Results.useCache, 'poleMovieInfo');

                        if numel(ip.Results.iFrame)>1
                            v1 = s.poleMovieInfo;
                        else
                            v1 = s.poleMovieInfo(iFrame);
                        end
                        if ~isempty(v1.xCoord) && ~isempty(iZ)
                            % Only show Detections in Z. 
                            zThick = 1;
                            tt = table(v1.xCoord(:,1), v1.yCoord(:,1), v1.zCoord(:,1), 'VariableNames', {'xCoord','yCoord','zCoord'});
                            valid_states = (tt.zCoord>=(iZ-zThick)*ZXRatio_z & tt.zCoord<=(iZ+zThick)*ZXRatio_z);
                            dataOut = tt{valid_states, :};

                            if isempty(dataOut) || numel(dataOut) <1 || ~any(valid_states)
                                dataOut = [];
                            end
                        else
                            dataOut = [];
                        end
                        if ~isempty(dataOut)
                            dataOutz = obj.convertProjection3D(dataOut, projAxis3D, ZXRatio);
                        end
                        varargout{iout} = dataOut;

                    case {'detectPolesMIP'}
                    
                        s = cached.load(obj.outFilePaths_{1, iChan}, '-useCache', ip.Results.useCache, 'poleMovieInfo');

                        if numel(ip.Results.iFrame)>1
                            v1 = s.poleMovieInfo;
                        else
                            v1 = s.poleMovieInfo(iFrame);
                        end
                        if ~isempty(v1.xCoord) && ~isempty(iZ)

                            tt = table(v1.xCoord(:,1), v1.yCoord(:,1), v1.zCoord(:,1),...
                                       'VariableNames', {'xCoord','yCoord','zCoord'});
                            valid_states = (tt.zCoord>=1 & tt.zCoord<=(obj.owner_.zSize_*ZXRatio_z));
                            
                            dataOut = tt{:, :};

                            if isempty(dataOut) || numel(dataOut) <1 || ~any(valid_states)
                                dataOut = [];
                            end
                        else
                            dataOut = [];
                        end
                        if ~isempty(dataOut)
                            dataOutz = obj.convertProjection3D(dataOut, projAxis3D, ZXRatio);
                        end
                        varargout{iout} = dataOutz;
                    otherwise
                end

            end
        end
        function output = getDrawableOutput(obj)
            n = 1;
            output = getDrawableOutput@DetectionProcess(obj);
            output(n).name='Poles by z';
            output(n).var = 'detectPolesZ';
            output(n).formatData=@DetectionProcess.formatOutput3D;
            output(n).defaultDisplayMethod=@(x) LineDisplay('Marker','.',...
                'LineStyle','none','MarkerSize',25,'Color', 'y');
            n = n + 1;
            colors = jet(numel(obj.owner_.channels_));
            output(n) = getDrawableOutput@DetectionProcess(obj);
            output(n).name='Poles';
            output(n).var = 'detectPolesMIP';
            output(n).formatData=@DetectionProcess.formatOutput3D;
            output(n).defaultDisplayMethod=@(x) LineDisplay('Marker','.',...
                'LineStyle','none','MarkerSize',35,'Color', 'y');
        end  
        

    end
    methods (Static)
        
        function name = getName()
            name = 'Pole Detection';
        end
        
        function h = GUI()
            h = @detectPolesProcessGUI;
        end
        
        function funParams = getDefaultParams(owner, varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir = ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1;
            funParams.InputImageProcessIndex = 0; % ?? (can we add some way to check what is availble.)
            funParams.OutputDirectory = [outputDir  filesep 'poles'];
            funParams.processFrames = 1:owner.nFrames_;
            
            funParams.Alpha = 0.05;
            funParams.scales = 3;
            funParams.showAll = false;
            funParams.printAll = false;
            funParams.isoOutput = false;
            funParams.type = 'simplex';
           
        end
    end    
end