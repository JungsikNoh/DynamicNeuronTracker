classdef LCHCellSegmentationProcess < SegmentationProcess
    % A concrete process for segmenting Live Cell Histology Phase Contrast Images
    % This is quick hack/rapid prototype, far from perfect...
    % may be place holder for better algorithm
    % this may also serve as just a prep process for converting to .png image format
    
    % Andrew R. Jamieson - Dec 2017
    
    methods
        function obj = LCHCellSegmentationProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = LCHCellSegmentationProcess.getName;
                super_args{3} = @segCellLCH_MD;
                if isempty(funParams)
                    funParams=LCHCellSegmentationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@SegmentationProcess(super_args{:});
        end
        
    end
    methods (Static)
        function name = getName()
            name = 'LCH Cell Segmentation';
        end
        function h = GUI()
            h= @LCHCellSegmentationProcessGUI;
        end
    
        function funParams = getDefaultParams(owner, varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner', @(x) isa(x, 'MovieData'));
            ip.addOptional('outputDir', owner.outputDirectory_, @ischar);
            ip.parse(owner, varargin{:})
            outputDir = ip.Results.outputDir;

            [~, imgFileNamePrefix] = fileparts(owner.outputDirectory_);
            funParams.imageFileNamePrefix = imgFileNamePrefix;

            % Set default parameters
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir];
            funParams.ProcessIndex = []; %Default is to use raw images

            funParams.preview = false ;
            funParams.usephasecontrastSeg = false; % even more unvetted version from previous lab members
            funParams.algorithm = 'basic'; % {'basic','phasecontrast','lever'}

            % verify all these....
            % Specific to the LEVER algorithm
            funParams.leverParams.CONSTANTS.imageData.PixelPhysicalSize = 0.325;
            funParams.leverParams.minimumRadius_um = 5;
            funParams.leverParams.maximumRadius_um = 200;
            funParams.leverParams.channels = 1; 
            funParams.leverParams.alphaLevels = 1;
            funParams.leverParams.convexRatio = 1.1;
        end
    end
end