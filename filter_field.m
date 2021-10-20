function filtered_field = filter_field(field,window,type)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This is a global filtering function which can perform a
% horizontal, vertical or temporal filter.
%
% INPUTS:
%
% field = 1,2,3 or 4D field.
%
% window = size of filtering window.
%
% type = '-t' :> temporal filter, will simply filter the last
%               dimension of field using a box filter of window
%               size window (must be odd). The end points that are
%               not filtered will be output as NaNs.
%
%       '-s' :> spatial 'circle' filter, will filter over two
%               horizontal dimensions using a box/circle filter
%               with number of points window. The field should be
%               two or three dimensions.
%
%       '-v' :> vertical filter, will filter in the vertical using
%               a box of number of points (at the moment; but this
%               should be improved since the dimension is not
%               linear). The field should be three dimensions
%
% OUTPUTS:
% 
% filtered_field = filtered_field of same size as field.
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 17/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

if (window>1)
switch type
  case '-t'  

    %%%%%%%%%%Time Filter%%%%%%%%%%%%%%%%%
    if (mod(window,2)==0)
        warning('time filter window must be odd...');
        warning('field not filtered!!!!!!!');
        filtered_field = field;
        return;
    end
    OneSide = (window-1)/2;
    filtered_field = zeros(size(field));
    filtered_field = NaN*filtered_field;
    if (length(find(size(field) == 1))>0)
        tL = length(filtered_field);
        for t=(OneSide+1):(tL-OneSide)
            filtered_field(t) = nanmean(field((t-OneSide):(t+OneSide)));%/window;
        end
        
    elseif (length(size(field)) == 2)
        tL = length(filtered_field(1,:));
        for t=(OneSide+1):(tL-OneSide)
            filtered_field(:,t) = nanmean(field(:,((t-OneSide):(t+OneSide))),2);%/window;
        end
        
    elseif (length(size(field)) == 3)
        tL = length(filtered_field(1,1,:));
        for t=(OneSide+1):(tL-OneSide)
            filtered_field(:,:,t) = nanmean(field(:,:,((t-OneSide):(t+OneSide))),3);%/window;
        end
    elseif (length(size(field)) == 4)
        tL = length(filtered_field(1,1,1,:));
        
        for t=(OneSide+1):(tL-OneSide)
            filtered_field(:,:,:,t) = nanmean(field(:,:,:,((t-OneSide):(t+OneSide))),4);%/window;
        end
    end
    
  case '-s'  
    
    %%%%%%%%%%Horizontal Filter%%%%%%%%%%%%%%%%%
if (length(find(size(field) == 1))>0)
        lx = length(field);
        diags = 1/window.*ones(max(size(field)),2*window);
        wdiags = -(window-1)/2:(window-1)/2;
        if (circ == 1)
            wdiags = [-(max(size(field))-1):(-(max(size(field))-1)+((window-1)/2-1)) wdiags ((max(size(field))-1)-((window-1)/2-1)):(max(size(field))-1)] 
        end
        filtx = spdiags(diags,wdiags,lx,lx);
        if (circ ~= 1)
        for i=1:(window-1)/2
            filtx(i,1) =0.5*(window+3-2*i)/window;
        end
        end
        filtered_field = filtx*field;
    elseif (length(size(field)) == 2)
        [lx ly]=size(field);
        diags = 1/window.*ones(max(size(field)),window);
        wdiags = -(window-1)/2:(window-1)/2;
        filtx = spdiags(diags, wdiags,lx,lx);
        filty = spdiags(diags, wdiags,ly,ly);
        for i=1:(window-1)/2
            filtx(i,1) =0.5*(window+3-2*i)/window;
            filty(i,1) =0.5*(window+3-2*i)/window;
        end
        filtered_field = filtx*field;
        filtered_field = filtered_field*(filty');
    elseif (length(size(field)) == 3)
        [lx,ly,lz]=size(field);
        diags = 1/window.*ones(max([lx ly]),window);
        wdiags = -(window-1)/2:(window-1)/2;
        filtx = spdiags(diags, wdiags,lx,lx);
        filty = spdiags(diags, wdiags,ly,ly);
        for i=1:(window-1)/2
            filtx(i,1) =0.5*(window+3-2*i)/window;
            filty(i,1) =0.5*(window+3-2*i)/window;
        end
        
        filtered_field = zeros(size(field));
        for z = 1:lz
        filtered_field(:,:,z) = filtx*field(:,:,z);
        filtered_field(:,:,z) = filtered_field(:,:,z)*(filty');
        end
    elseif (length(size(field)) == 4)
        [lx,ly,lz,lt]=size(field);
        diags = 1/window.*ones(max([lx ly]),window);
        wdiags = -(window-1)/2:(window-1)/2;
        filtx = spdiags(diags, wdiags,lx,lx);
        filty = spdiags(diags, wdiags,ly,ly);
        for i=1:(window-1)/2
            filtx(i,1) =0.5*(window+3-2*i)/window;
            filty(i,1) =0.5*(window+3-2*i)/window;
        end
        
        filtered_field = zeros(size(field));
        for t=1:lt
            for z = 1:lz
                filtered_field(:,:,z,t) = filtx*field(:,:,z,t);
                filtered_field(:,:,z,t) = filtered_field(:,:,z,t)*(filty');
            end    
        end
        
    end

  case '-v'  

    %%%%%%%%%%Vertical Filter%%%%%%%%%%%%%%%%%
        
        [lx,ly,lz]=size(field);
        diags = 1/window.*ones(max([lx ly lz]),window);
        wdiags = -(window-1)/2:(window-1)/2;
        filtz = spdiags(diags, wdiags,lz,lz);
        for i=1:(window-1)/2
            filtz(i,1) =0.5*(window+3-2*i)/window;
        end        
        filtered_field = zeros(size(field));
        for x = 1:lx
        filtered_field(x,:,:) = permute(filtz*permute(field(x,:,:),[3 ...
                            2 1]),[3 2 1]);
        end

  otherwise
    warning('Unknown filter type...');
end
else
filtered_field = field;
end
end

