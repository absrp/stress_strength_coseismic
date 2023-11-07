clear; close all;
% List all KMZ files in the directory
fileList = dir('*.kmz');


% Loop through the KMZ files
for i = 1:length(fileList)
    % Initialize an empty structure to store the kmz_struct data
    kmzData = struct();
    lines_struct = struct('Name', {}, 'Geometry', {}, 'Lat', [], 'Lon', [], 'Description',{},'BoundingBox',[],'Folder',{}, ...
   'Color',{});
    kmzFileName = fileList(i).name;
    kmz_struct = kmz2struct(kmzFileName);
    % Loop through the elements of kmz_struct
    for p = 1:length(kmz_struct)
        element = kmz_struct(p);
        
        % Check if the Geometry field is 'Line'
        if strcmp(element.Geometry, 'Line') %& strcmp(element.Folder,'/Ruptures - Principal')
            lines_struct(end + 1).Name = element.Name;
            lines_struct(end).Geometry = element.Geometry;
            lines_struct(end).Lat = element.Lat;
            lines_struct(end).Lon = element.Lon;
            lines_struct(end).Description = element.Description;
            lines_struct(end).BoundingBox = element.BoundingBox;
            lines_struct(end).Folder = element.Folder;
            lines_struct(end).Color = element.Color;
        end
    end

    % Filter out distributed and primary lines to leave ECS line only
    mask = ~strcmp({lines_struct.Folder}, '/Ruptures - Distributed');
    lines_struct = lines_struct(mask);

    mask = ~strcmp({lines_struct.Folder}, '/Ruptures - Principal');
    lines_struct = lines_struct(mask);
    
    % remove unnecessary fields for shapefile
    lines_struct = rmfield(lines_struct, 'Color');
    lines_struct = rmfield(lines_struct, 'Description');

    % prepare shapefile name
    firstString = lines_struct(1).Name;
    secondString = '_ECS.shp';
    shapename = [firstString, secondString];
   
     if strcmp(lines_struct(1).Name, 'Ridgecrest1')
        mask_secondary =  ~strcmp({lines_struct.Folder}, '/Ruptures (RUP_DS_ID=145) - Distributed');
        lines_struct = lines_struct(mask_secondary);
    else
     end

         if strcmp(lines_struct(1).Name, 'Ridgecrest1')
        mask_secondary =  ~strcmp({lines_struct.Folder}, '/Ruptures (RUP_DS_ID=145) - Principal');
        lines_struct = lines_struct(mask_secondary);
    else
         end

    if strcmp(lines_struct(1).Name, 'Ridgecrest1')
        mask_secondary =  ~strcmp({lines_struct.Folder}, '/Ruptures (RUP_DS_ID=132) - Principal');
        lines_struct = lines_struct(mask_secondary);
    else
    end

    if strcmp(lines_struct(1).Name, 'Ridgecrest1')
        mask_secondary =  ~strcmp({lines_struct.Folder}, '/Ruptures (RUP_DS_ID=132) - Distributed');
        lines_struct = lines_struct(mask_secondary);
    else
    end

        if strcmp(lines_struct(1).Name, 'Ridgecrest2')
        mask_secondary =  ~strcmp({lines_struct.Folder}, '/Ruptures (RUP_DS_ID=132) - Principal');
        lines_struct = lines_struct(mask_secondary);
    else
    end

    if strcmp(lines_struct(1).Name, 'Ridgecrest2')
        mask_secondary =  ~strcmp({lines_struct.Folder}, '/Ruptures (RUP_DS_ID=132) - Distributed');
        lines_struct = lines_struct(mask_secondary);
    else
    end


    

     if strcmp(lines_struct(1).Name, 'Ridgecrest2')
        mask_secondary =  ~strcmp({lines_struct.Folder}, '/Ruptures (RUP_DS_ID=145) - Distributed');
        lines_struct = lines_struct(mask_secondary);
    else
     end


         if strcmp(lines_struct(1).Name, 'Ridgecrest2')
        mask_secondary =  ~strcmp({lines_struct.Folder}, '/Ruptures (RUP_DS_ID=145) - Principal');
        lines_struct = lines_struct(mask_secondary);
    else
    end

 shapewrite(lines_struct,shapename)
end

