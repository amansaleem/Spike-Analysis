function [outfilename] = BuildFileName(varargin)

% [outfilename] = BuildFileName(varargin) is a function that takes a
% variable input and builds an appropriate string file name separated by
% underscores

% 2010-10 written by ND.

nargin = size(varargin,2);
if nargin<1, error('No file structure input. Can''t build a filename.'); end
outfilename = sprintf('%s',varargin{1});
for iargin = 2:nargin
    if isnumeric(varargin{iargin})
        if iargin == 2
            outfilename = sprintf('%s\\%d',outfilename,varargin{iargin});
        else
            outfilename = sprintf('%s_%d',outfilename,varargin{iargin});
        end
    else % string
        if iargin == 2
            outfilename = sprintf('%s\\%s',outfilename,varargin{iargin});
        else
            outfilename = sprintf('%s_%s',outfilename,varargin{iargin});
        end
    end
end


