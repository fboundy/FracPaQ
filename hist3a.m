function [N, c] = hist3a(XY, Nxy)
%% Copyright
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the
% "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
% 
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF COnTracesACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.
%
%   This script provides the basic functionality of the hist3 XY binning
%   functionality of the hist3 function without the graphics output. hist3
%   requires the Statistics and Machine Learning Toolbox
%
%   Created by Francis Boundy - April 2021
%

    X = XY(:,1);    % vector of X coordinates
    Y = XY(:,2);    % vector of X coordinates
    
    Nx = Nxy(:,1);  % number of X bins
    Ny = Nxy(:,2);  % number of Y bins

    % vectors of the bin that each point occupies
    I = floor(min((X - min(X)) / (max(X) - min(X)) * Nx + 1, Nx));
    J = floor(min((Y - min(Y)) / (max(Y) - min(Y)) * Ny + 1, Ny));
    
    % initialise the ouput matrix
    N = zeros(Nx,Ny);
    
    % count the number of points in each bin
    for i = 1:length(I)
        N(I(i),J(i)) = N(I(i),J(i)) + 1;
    end
    
    % create the matrix of centre points
    dx = (max(X) - min(X)) / Nx;
    cx = min(X)+dx/2:dx:max(X);
    dy = (max(Y) - min(Y)) / Ny;
    cy = min(Y)+dy/2:dy:max(Y);
    c = {cx, cy};

end