function guiFracPaQ2Dpattern(traces, numPixelsPerMetre, ...
                             nBlocks, ... 
                             flag_intensitymap, flag_densitymap, ... 
                             flag_triangle, flag_showcircles, ...
                             nCircles, flag_revY, flag_revX, sColour, nPixelsItoY, grid_angle)
%   guiFracPaQ2Dpattern.m
%       calculates and plots statistics of trace segment patterns
%
%   David Healy
%   July 2014
%   d.healy@abdn.ac.uk
%
%   Modified by:
%   Nikolai Andrianov 
%   January 2019 
%
%   Modified by:
%   David Healy 
%   February 2019 
%
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
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

global sTag ; 

if numPixelsPerMetre > 0
    sUnits = ' metres' ; 
else 
    sUnits = ' pixels' ; 
end 


[xMin, xMax, yMin, yMax] = getMapLimits(traces) ; 
%iLength = xMax - xMin ;
%jLength = yMax - yMin ;

TraceLimits = getRotatedMapLimits(traces, grid_angle);

x00 = TraceLimits(1,1);
x01 = TraceLimits(1,2);
x10 = TraceLimits(1,3);
x11 = TraceLimits(1,4);

y00 = TraceLimits(2,1);
y01 = TraceLimits(2,2);
y10 = TraceLimits(2,3);
y11 = TraceLimits(2,4);

disp('Rotated');
disp('-------');
disp(['P00: ', num2str(x00, '%12.0F'), ', ', num2str(y00, '%12.0f')]) ;
disp(['P10: ', num2str(x10, '%12.0F'), ', ', num2str(y10, '%12.0f')]) ;
disp(['P01: ', num2str(x01, '%12.0F'), ', ', num2str(y01, '%12.0f')]) ;
disp(['P11: ', num2str(x11, '%12.0F'), ', ', num2str(y11, '%12.0f')]) ;
disp(' ');
xv = [x10, x00, x01];
yv = [y10, y00, y01];


jLength = hypot((x01 - x00),(y01 - y00));
iLength = hypot((x10 - x00),(y10 - y00));

numTraces = length(traces) ;
traceLengths = [ traces.segmentLength ] ;
% mapArea = (xMax - xMin) * (yMax - yMin) ; 
mapArea = iLength * jLength;

if flag_intensitymap || flag_densitymap || flag_showcircles
    
%       % apply circular scan lines and windows
%       % calculate I, D for a selected set of points
%       %   define circle centres
%         xNumCircle = nCircles ;
%         yNumCircle = nCircles ;
%         xDeltaCircle = ( xMax - xMin ) / ( xNumCircle + 1 ) ;
%         yDeltaCircle = ( yMax - yMin ) / ( yNumCircle + 1 ) ;
%     
%         %   set circle radius, as function of x and y increments
%         rCircle = 0.99 * min(xDeltaCircle, yDeltaCircle) / 2 ;
    
    %% X dimension is longer than y
    if iLength >= jLength        
        jNumCircle = nCircles ;

        jDeltaCircle = jLength / ( jNumCircle + 1 ) ;
        %   set circle radius, as function of y increments
        rCircle = 0.99 * jDeltaCircle / 2 ;

        % Calculate number of circles in x to match the image dimensions
        % based on the radius determined by the number of circles in y
        iDeltaCircle=(rCircle*2)/0.99;
        iNumCircle=floor((iLength / iDeltaCircle) - 1);
        
    %% Y dimension is longer than x
%    elseif ( xMax - xMin ) < ( yMax - yMin )
    else       
        iNumCircle = nCircles ;
        
        iDeltaCircle =  iLength  / ( iNumCircle + 1 ) ;
        %   set circle radius, as function of x increments
        rCircle = 0.99 * iDeltaCircle / 2 ;
        % Calculate number of circles in y to match the image dimensions
        % based on the radius determined by the number of circles in x
        jDeltaCircle=(rCircle*2)/0.99;
        jNumCircle=floor((jLength / jDeltaCircle) - 1);
        
    %% X and y dimensions are equal
%    elseif ( xMax - xMin ) == ( yMax - yMin )
%        %   define circle centres
%        iNumCircle = nCircles ;
%        jNumCircle = nCircles ;
%        iDeltaCircle = ( iLength ) / ( iNumCircle + 1 ) ;
%        jDeltaCircle = ( jLength ) / ( jNumCircle + 1 ) ;
%        %   set circle radius, as function of x and y increments
%        rCircle = 0.99 * min(iDeltaCircle, jDeltaCircle) / 2 ;
        
    end
    
    disp(' ') ;
    disp('Circular scan windows...') ;
    disp(['Circle increment in I: ', num2str(iDeltaCircle, '%8.2E'), sUnits]) ;
    disp(['Circle increment in J: ', num2str(jDeltaCircle, '%8.2E'), sUnits]) ;
    disp(['Circle radius: ', num2str(rCircle, '%8.2E'), sUnits]) ;
    rCircleMetres = rCircle ;
    I = zeros(jNumCircle, iNumCircle) ;
    D = zeros(jNumCircle, iNumCircle) ;
    
end ;

if flag_showcircles
    
    f = figure ;
    %   mapped lines and circle centres
    hold on ;
    for k = 1:numTraces
        plot( [ traces(k).Node.x ]', [ traces(k).Node.y ]', 'LineWidth', 0.75, 'Color', sColour) ;
    end ;

    line(xv,yv,'LineWidth', 0.25, 'Color', 'r');
    %   for each circle centre
    for i = 1:iNumCircle
 %       xCentreCircle = xMin + iDeltaCircle * i ;
        for j = 1:jNumCircle
            
 %           yCentreCircle = yMin + jDeltaCircle * j ;
            [xCentreCircle, yCentreCircle] =  ij2xy(i, j, x00, x01, x10, y00, y01, y10, jDeltaCircle, jDeltaCircle);
            %   *** need to draw the circles actual size here
            %   use rectangle() with 'Curvature' at [1 1] 
            pos = [ xCentreCircle - rCircle, yCentreCircle - rCircle, 2*rCircle, 2*rCircle ] ; 
            rectangle('Position', pos, 'Curvature', [1 1], 'LineWidth', 0.25, 'EdgeColor', 'r') ; 
        end ;
    end ;
    hold off ;
    axis on equal ;
    box on ;
 
    disp(['P00: ', num2str(x00, '%12.8E'), ', ', num2str(y00, '%12.8E')]) ;
    disp(['P10: ', num2str(x10, '%12.8E'), ', ', num2str(y10, '%12.8E')]) ;
    disp(['P01: ', num2str(x01, '%12.8E'), ', ', num2str(y01, '%12.8E')]) ;
 
    
    if x00 < x01
        xlim([x00 x10+x01-x00]) ;
        ylim([y01 y10]);
    else
        xlim([x01 x10]);
        ylim([y00 y10+y01-y00]) ;
    end;
    
    if flag_revX
        set(gca, 'XDir', 'reverse') ;
    end ;
    if flag_revY
        set(gca, 'YDir', 'reverse') ;
    end ;
    xlabel(['X,', sUnits]) ;
    ylabel(['Y,', sUnits]) ;
    title({['Mapped trace segments, n=', num2str(length(traceLengths))];''}) ;
    
    %   save to file
    guiPrint(f, 'FracPaQ2D_scancircle') ;
    
end ;


if flag_intensitymap || flag_densitymap
    
    hWait = waitbar(0, 'Calculating scan circle intersections...', 'Name', 'Intensity/Density maps') ;
    nCircle = 0 ;
    disp(['p00:' , num2str(x00,'%10.0f'), ',', num2str(y00,'%10.0f'),' p01:', num2str(x01,'%10.0f'), ',',  num2str(y01,'%10.0f'), ' p10:', num2str(x10,'%10.0f'), ',', num2str(y10,'%10.0f')]);
    disp(['Ni: ' , num2str(iNumCircle,'%4.0f'), ' Nj: ', num2str(jNumCircle,'%4.0f'),' Li: ', num2str(iLength,'%10.0f'), ' Lj: ',  num2str(jLength,'%10.0f'), ' p10:', num2str(x10,'%10.0f'), ',', num2str(y10,'%10.0f')]);

    %   for each circle centre
    for i = 1:iNumCircle
        
 %       xCentreCircle = xMin + iDeltaCircle * i ;
        
        for j = 1:jNumCircle
            
            nCircle = nCircle + 1 ;
            
            waitbar(nCircle/(iNumCircle*jNumCircle), hWait, 'Calculating scan circle intersections...') ;
            
 %           yCentreCircle = yMin + jDeltaCircle * j ;
            [xCentreCircle, yCentreCircle] =  ij2xy(i, j, x00, x01, x10, y00, y01, y10, iDeltaCircle, jDeltaCircle);
            disp(['i: ', num2str(i,'%4.0f'), ' j: ', num2str(j,'%4.0f'), ' xy:' , num2str(xCentreCircle,'%10.0f'), ', ', num2str(yCentreCircle,'%10.0f')]);
            n = 0 ;
            m = 0 ;

            for k = 1:numTraces
                
                for s = 1:traces(k).nSegments
                    
                    bPoint1Inside = false ;
                    bPoint2Inside = false ;
                    
                    %           first end of line
                    rPoint = sqrt( ( traces(k).Segment(s).Point1(1) - xCentreCircle )^2 ...
                        + ( traces(k).Segment(s).Point1(2) - yCentreCircle )^2 ) ;
                    if rPoint < rCircleMetres
                        m = m + 1 ;
                        bPoint1Inside = true ;
                    end ;
                    
                    %           second end of line
                    rPoint = sqrt( ( traces(k).Segment(s).Point2(1) - xCentreCircle )^2 ...
                        + ( traces(k).Segment(s).Point2(2) - yCentreCircle )^2 ) ;
                    if rPoint < rCircleMetres
                        m = m + 1 ;
                        bPoint2Inside = true ;
                    end ;
                    
                    %           find any intersections of line with circle
                    if ( bPoint1Inside && bPoint2Inside )
                        continue ;
                        
                    elseif bPoint1Inside
                        n = n + 1 ;
                        
                    elseif bPoint2Inside
                        n = n + 1 ;
                        
                    else
                        
                        dx = traces(k).Segment(s).Point2(1) - traces(k).Segment(s).Point1(1) ;
                        dy = traces(k).Segment(s).Point2(2) - traces(k).Segment(s).Point1(2) ;
                        dr = sqrt( dx^2 + dy^2 ) ;
%                         Det = traces(k).Segment(s).Point1(1) * traces(k).Segment(s).Point2(2) - ...
%                             traces(k).Segment(s).Point2(1) * traces(k).Segment(s).Point1(2) - ...
%                             xCentreCircle * ( traces(k).Segment(s).Point2(2) - traces(k).Segment(s).Point1(2) ) + ...
%                             yCentreCircle * ( traces(k).Segment(s).Point2(1) - traces(k).Segment(s).Point1(1) ) ;
                         
                        x1new = traces(k).Segment(s).Point1(1) - xCentreCircle ; 
                        x2new = traces(k).Segment(s).Point2(1) - xCentreCircle ;
                        y1new = traces(k).Segment(s).Point1(2) - yCentreCircle ; 
                        y2new = traces(k).Segment(s).Point2(2) - yCentreCircle ; 
                        
                        Det = x1new * y2new - x2new * y1new ;
                        
                        deltaSecant = rCircleMetres^2 * dr^2 - Det^2 ;
                        
                        if deltaSecant > 0
                            
                            %   this test checks to see if the adjusted line end
                            %   points are either side of the circle centre
                            %   (to correct a bug where finite line
                            %   segments far away from circle, but
                            %   aligned to intersect were being reported as
                            %   intersections)
                            if ( ( x1new < 0 && x2new >= 0 ) || ...
                                 ( x1new > 0 && x2new <= 0 ) ) && ...    
                               ( ( y1new < 0 && y2new >= 0 ) || ...
                                 ( y1new > 0 && y2new <= 0 ) )    

                                n = n + 2 ;

                            end ; 
                                
                        end ;
                        
                    end ;
                    
                end ;
                
            end ;
            
            %       calculate I
            I(j, i) = n / ( 4 * rCircleMetres ) ;
            
            %       calculate D
            D(j, i) = m / ( 2 * pi * rCircleMetres^2 ) ;
            
            % %       calculate MTL
            %         if m > 0
            %             MTL(j, i) = 0.5 * pi * rCircleMetres * ( n / m ) ;
            %         else
            %             MTL(j, i) = sqrt(xMax^2 + yMax^2) ;
            %         end ;
            
        end ;
        
    end ;
    
    close(hWait) ;
    
end ;

if flag_triangle
    
    %   connectivity plot
    k = 0 ;
    for i = 1:numTraces
        
        for j = 1:traces(i).nSegments
            
            k = k + 1 ;
            segmentsxy(k, :) = [ traces(i).Segment(j).Point1(1), traces(i).Segment(j).Point1(2), ...
                traces(i).Segment(j).Point2(1), traces(i).Segment(j).Point2(2) ] ;
            
        end ;
        
    end ;
    
    [cY, cX, cI] = getConnectivity(traces, segmentsxy, xMin, yMin, xMax, yMax, ...
                                    nBlocks, flag_revX, flag_revY, ...
                                    numPixelsPerMetre, nPixelsItoY) ;
    cTot = cY + cX + cI ;
    disp('Connectivity...') ;
    disp(['Y:X:I = ', ...
        num2str(cY/cTot, '%5.2f'), ':', ...
        num2str(cX/cTot, '%5.2f'), ':', ...
        num2str(cI/cTot, '%5.2f')]) ;
    disp(['Total number of I-Y-X nodes: ', num2str(cTot, '%4.0f')]) ;
    xX = 0.5 * ( 2 * cX + cI ) / cTot ;
    yX = sqrt(3) * 0.5 * cI / cTot ;
    triData = [ (cY/cTot)*100, (cX/cTot)*100, (cI/cTot)*100 ] ;

    %   equations 5c in Sanderson & Nixon, 2015 
    triLinesCL2 = [ 50, 0, 50 ; 0, 33, 67 ] ; 
    triLinesCL357 = [ (8.3/9.3)*100, 0, (1/9.3)*100 ; 0, (0.8925/1.8925)*100, (1/1.8925)*100 ] ; 
    
    f = figure ;
    hold on ;
    triplot() ; 
    tripts(triData, '', 0, sColour) ; 
    tripts(triLinesCL2, '', 1, sColour) ; 
    tripts(triLinesCL357, '', 1, sColour) ; 
    tripts(triLinesCL2(2,:), '\it C_L\rm = 2', 2, sColour) ; 
    tripts(triLinesCL357(2,:), '\it C_L\rm = 3.57', 2, sColour) ; 
    hold off ;
    
    % Title with fractions of Y, X, I
    title({'Connectivity of traces'; ...
           ['Y:X:I = ', ...
            num2str(cY/cTot, '%5.2f'), ':', ...
            num2str(cX/cTot, '%5.2f'), ':', ...
            num2str(cI/cTot, '%5.2f')];''}) ;
    
    %   save to file
    guiPrint(f, 'FracPaQ2D_IYXtriangle') ;
    
    %   write a file of counts
    fn = ['FracPaQ2Dconnectivity', sTag, '.txt'] ; 
    fidConn = fopen(fn, 'wt') ;
    fprintf(fidConn, '%s %i\n', 'I', cI) ;
    fprintf(fidConn, '%s %i\n', 'X', cX) ;
    fprintf(fidConn, '%s %i\n', 'Y', cY) ;
    fclose(fidConn) ;
    
end ;

if flag_intensitymap || flag_densitymap
    
    if flag_intensitymap
        %   plot trace intensity, I
        f = figure ;
%         contourf(X2, Y2, Inew, nContours) ;
        xv = (xMin+iDeltaCircle):iDeltaCircle:(xMax-iDeltaCircle) ; 
        yv = (yMin+jDeltaCircle):jDeltaCircle:(yMax-jDeltaCircle) ; 
        
        imagesc(xv, yv, I);
                
        ax = gca ; 
        ax.YDir = 'normal' ;
        axis on equal ;
        box on ;
        xlim([xMin xMax]) ;
        ylim([yMin yMax]) ;
        if flag_revX
            set(gca, 'XDir', 'reverse') ;
        end ;
        if flag_revY
            set(gca, 'YDir', 'reverse') ;
        end ;
        title({'Estimated Intensity of segments (P21)';''}) ;
        xlabel(['X,', sUnits]) ;
        ylabel(['Y,', sUnits]) ;
        c = colorbar('location', 'southoutside') ;
        if numPixelsPerMetre > 0
            c.Label.String = 'Intensity of segments (P21), metre^{-1}' ; 
        else 
            c.Label.String = 'Intensity of segments (P21), pixel^{-1}' ; 
        end 
        cmap = colormap('hot') ; 
        cmap = flipud(cmap) ; 
        colormap(cmap) ; 
        caxis([0 max(max(I))]) ; 

        %   save to file
        guiPrint(f, 'FracPaQ2D_intensityP21') ;
        
        %   print Intensity for whole map area
        disp(' ') ; 
        disp('Intensity for whole map:') ; 
        if numPixelsPerMetre > 0
            disp([num2str(sum(traceLengths)/mapArea), ' metre^-1']) ; 
        else 
            disp([num2str(sum(traceLengths)/mapArea), ' pixel^-1']) ; 
        end ; 
            
    end ;
    
    if flag_densitymap
        %   plot trace density, D
        f = figure ;
%         contourf(X2, Y2, Dnew, nContours) ;
        
        xv = (xMin+iDeltaCircle):iDeltaCircle:(xMax-iDeltaCircle) ; 
        yv = (yMin+jDeltaCircle):jDeltaCircle:(yMax-jDeltaCircle) ; 
        
        imagesc(xv, yv, D) ; 
        
        ax = gca ; 
        ax.YDir = 'normal' ;

        axis on equal ;
        box on ;
        xlim([xMin xMax]) ;
        ylim([yMin yMax]) ;
        if flag_revX
            set(gca, 'XDir', 'reverse') ;
        end ;
        if flag_revY
            set(gca, 'YDir', 'reverse') ;
        end ;
        title({'Estimated Density of segments (P20)';''}) ;
        xlabel(['X,', sUnits]) ;
        ylabel(['Y,', sUnits]) ;
        c = colorbar('southoutside') ; 
        if numPixelsPerMetre > 0
            c.Label.String = 'Density of segments (P20), metre^{-2}' ; 
        else
            c.Label.String = 'Density of segments (P20), pixel^{-2}' ; 
        end
        cmap = colormap('hot') ; 
        cmap = flipud(cmap) ; 
        colormap(cmap) ; 
        caxis([0 max(max(D))]) ; 
        
        %   save to file
        guiPrint(f, 'FracPaQ2D_densityP20') ;
        
        %   print Density for whole map area
        disp(' ') ; 
        disp('Density for whole map:') ; 
        if numPixelsPerMetre > 0
            disp([num2str(length(traceLengths)/mapArea), ' metre^-2']) ; 
        else 
            disp([num2str(length(traceLengths)/mapArea), ' pixel^-2']) ; 
        end ; 

    end ;
    
    %   write I and D data to file
    fn1 = ['FracPaQ2Dintensity', sTag, '.txt'] ; 
    fn2 = ['FracPaQ2Ddensity', sTag, '.txt'] ; 
    fidIntensity = fopen(fn1, 'wt') ;
    fidDensity = fopen(fn2, 'wt') ;
    for i = 1:iNumCircle
        
        xCentreCircle = xMin + iDeltaCircle * i ;
        
        for j = 1:jNumCircle
            
            yCentreCircle = yMin + jDeltaCircle * j ;
            
            fprintf(fidIntensity, '%8.1f %8.1f %14.8f\n', xCentreCircle, yCentreCircle, I(j, i)) ;
            fprintf(fidDensity, '%8.1f %8.1f %14.8f\n', xCentreCircle, yCentreCircle, D(j, i)) ;
            
        end ;
        
    end ;
    fclose(fidIntensity) ;
    fclose(fidDensity) ;
    
end ;


end 

function [x00, x01, x10, y00, y01, y10] = rotateGrid(x0, x1, y0, y1, theta)

    lx = x1 - x0;
    ly = y1 - y0;

    if theta < 0
        phi = deg2rad(90 - theta);
    else
        phi = deg2rad(theta);
    end;
    
    dx0 = ly * sin(phi) * cos(phi);
    dx1 = lx * sin(phi)^2;
    dy0 = ly * sin(phi)^2;
    dy1 = lx * sin(phi) * cos(phi);
    
    x00 = x0 - dx0 ;
    y00 = y0 + dy0 ;
    x01 = x0 + dx1 ;
    y01 = y1 + dy1 ;
    x10 = x1 - dx1 ;
    y10 = y0 - dy1 ;
end

function [x, y] = ij2xy(i, j, x00, x01, x10, y00, y01, y10, di, dj)
% returns a pair of x, y coordinates for grid index i, j in a grid defined
% with origin (x00, y00), x axis (x10, y10) and y axis (x01, y01) and
% resolution di and dj
    lj = hypot((x01 - x00),(y01 - y00));
    li = hypot((x10 - x00),(y10 - y00));

    dxi = di * (x10 - x00) / li;
    dxj = dj * (x01 - x00) / lj;
    
    dyi = di * (y10 - y00) / li;
    dyj = dj * (y01 - y00) / lj;
    
    x = x00 + i * dxi + j * dxj;
    y = y00 + i * dyi + j * dyj;
end

