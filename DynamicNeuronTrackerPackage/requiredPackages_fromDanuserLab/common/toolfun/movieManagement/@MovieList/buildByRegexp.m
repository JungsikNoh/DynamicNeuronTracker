function ML = buildByRegexp(filter,outputDirectory)
% Build a MovieList using a regular expression.
%
% INPUT
% filter - regular expression for directories in the current directory
% which should contain .mat files of the same name as their directory
% outputDirectory - outputDirectory for the MovieList constructor
%
% OUTPUT
% ML - A MovieList containing MovieData objects that matches filter

% Mark Kittisopikul, March 2018
% Goldman Lab
% Northwestern University

    if(nargin < 2)
        outputDirectory = pwd;
    end
    D = dir;
    D = D([D.isdir]);
    D = D(~cellfun(@(x) isempty(regexp(x,filter, 'once')),{D.name}));
    movieDataFileNames = strcat(pwd,filesep,{D.name},filesep,{D.name},'.mat');
    ML = MovieList(movieDataFileNames,outputDirectory);
    ML.sanityCheck;
end