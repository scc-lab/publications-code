function [Header,Data] = LoadData_Amira(Pathname)
% LoadData_Amira - Loads 3D Amira data set
%
% USAGE:
% [Header,Data] = LoadData_Amira(Pathname)
%
% VARIABLES:
% Pathname = full path name for a '.am' file
% Header = relevant hear information
% Data = 3D matrix of image data
%
% NOTES: (1) Test for v5.2.2 files
%        (2) Doesn't support ascii, but this should be simple to implement
%        (3) Tested on a limited number of, but no warranty; use at your own risk.
%        (4) Shawn Stapleton, Nov 2012
%        (5) Modified to support multiple components by Rushikesh
%            Kamalapurkar, Dec 2023

    Header = [];
    Data = [];
    isUniform = 0;
    
    fid = fopen(Pathname);
    
    if( fid<=0 )
        disp(['Unable to open file: ' Pathname]);
        return
    end
    
    % Load Ascii header data
    AsciiData = strtrim(fgetl(fid));
    
    if(strcmp(AsciiData, '# AmiraMesh BINARY-LITTLE-ENDIAN 2.1')==0 )
        disp(['Not a valid Amira File']);
        return
    end
    
    % Load remainder of header Ascii Header and store pertinent information    
    AsciiData = strtrim(strtrim(fgetl(fid)));
    while(strcmp(AsciiData,'# Data section follows')==0)
        
        % Ignore blank lines
        if( length(AsciiData) == 0 )
            AsciiData = strtrim(fgetl(fid));
            continue;
        end
        
        % Get Lattice Dimensions
        if( findstr(AsciiData,'define Lattice')==1)
            vals = sscanf(AsciiData,'define Lattice %d %d %d');
            Header.xSize = vals(1);
            Header.ySize = vals(2);
            Header.zSize = vals(3);
            clear vals
            
        % if Materials File type, load material properties
        elseif( findstr(AsciiData,'Materials {')==1)
            openBracketCount = 1;
            MaterialCount = 1;

            while(openBracketCount > 0 )                
                AsciiData = strtrim(fgetl(fid));
                
               % Start of a new material
               if( findstr(AsciiData,'{') >= 1 )                    
                    openBracketCount = openBracketCount + 1;
                    
                    % Save Name
                    Header.Materials{MaterialCount}.Name = sscanf(AsciiData,'%s {');

               % Material Colour    
               elseif( findstr(AsciiData,'Color') >= 1 )
                    Colours = sscanf(AsciiData,'Color %f %f %f');
                    Header.Materials{MaterialCount}.Colours  = Colours';
                    clear Colours 
               % End of a material
               elseif( strcmp(AsciiData,'}') == 1 )
                    openBracketCount = openBracketCount - 1;
                    MaterialCount = MaterialCount + 1;                
                end
                
                
            end
            
            clear MaterialCount openBracketCount
                        
        % Get Bounding Box
        elseif( findstr(AsciiData,'BoundingBox')==1)
            vals = sscanf(AsciiData,'BoundingBox %g %g %g %g %g %g');
            Header.BoundingBox.xMin = vals(1);
            Header.BoundingBox.xMax = vals(2);
            Header.BoundingBox.yMin = vals(3);
            Header.BoundingBox.yMax = vals(4);
            Header.BoundingBox.zMin = vals(5);
            Header.BoundingBox.zMax = vals(6);
            
        % is Uniform Grid?
        elseif( findstr(AsciiData,'CoordType "uniform"')==1)
            isUniform = 1;
            
        elseif( findstr(AsciiData,'Lattice') >= 1 )
            
            % Determine Data Type
            if(findstr(AsciiData,'byte') >=1 ) 
                DataType = 'uint8';
            elseif(findstr(AsciiData,'float') >=1 ) 
                DataType = 'single';
            elseif(findstr(AsciiData,'Labels') >=1 ) 
                DataType = 'single';
            elseif(findstr(AsciiData,'short') >=1 ) 
                DataType = 'int16';
            elseif(findstr(AsciiData,'ushort') >=1 ) 
                DataType = 'uint16';
            elseif(findstr(AsciiData,'int') >=1 ) 
                DataType = 'int32';
            elseif(findstr(AsciiData,'double') >=1 ) 
                DataType = 'double';
            else
                DataType = 'int8';
            end
                
             % Determine if compressed
            if(findstr(AsciiData,'HxByteRLE') >=1 )
                Encoding = 'HxByteRLE';
                ind1 = findstr(AsciiData,'HxByteRLE') + length('HxByteRLE');
                ind2 = findstr(AsciiData(ind1:end),')');
                NumBytes = str2num(AsciiData(ind1+1:ind1+ind2-2));
            elseif(findstr(AsciiData,'HxZip') >=1 )
                Encoding = 'HxZip';
                ind1 = findstr(AsciiData,'HxZip') + length('HxZip');
                ind2 = findstr(AsciiData(ind1:end),')');
                NumBytes = str2num(AsciiData(ind1+1:ind1+ind2-2));
            else
                Encoding = 'None';
                NumBytes = inf;
            end
            
        elseif( findstr(AsciiData,'Content') >= 1 )
            
            % Determine Data Type
            if(findstr(AsciiData,'byte') >=1 ) 
                DataType = 'uint8';
            elseif(findstr(AsciiData,'float') >=1 ) 
                DataType = 'single';
            elseif(findstr(AsciiData,'Labels') >=1 ) 
                DataType = 'single';
            elseif(findstr(AsciiData,'short') >=1 ) 
                DataType = 'int16';
            elseif(findstr(AsciiData,'ushort') >=1 ) 
                DataType = 'uint16';
            elseif(findstr(AsciiData,'int') >=1 ) 
                DataType = 'int32';
            elseif(findstr(AsciiData,'double') >=1 ) 
                DataType = 'double';
            else
                DataType = 'int8';
            end
                
        end

        % Find number of components
        if(findstr(AsciiData,'[2]') >=1 )
            Components = 2;
        elseif(findstr(AsciiData,'[3]') >=1 )
            Components = 3;
        else
            Components = 1;
        end
        AsciiData = strtrim(fgetl(fid));
    end

    % Read the @1 marker
	AsciiData = strtrim(strtrim(fgetl(fid)));
     
    
    % Decode the data
    if( strcmp(Encoding,'HxByteRLE') ) 
        
        Data = fread( fid,NumBytes, DataType);    
    	Data = RLEDecodeEx(Data,Header.xSize*Header.ySize*Header.zSize,DataType);
        
    elseif( strcmp(Encoding,'HxZip') ) 
         Data = fread( fid,NumBytes, 'int8');
         Data = dunzip(Data,DataType);
    else  
        Data = fread( fid,NumBytes, DataType);
    end
    
    eval(['Data=',DataType,'(Data);']);
     
    % Close the file
    fclose(fid);
    

    
    % Calculate resolution
    Header.PixelSpacing(1)   = (Header.BoundingBox.xMax-Header.BoundingBox.xMin)./(Header.xSize-1);
    Header.PixelSpacing(2)   = (Header.BoundingBox.yMax-Header.BoundingBox.yMin)./(Header.ySize-1);
    Header.SliceThickness   = (Header.BoundingBox.zMax-Header.BoundingBox.zMin)./(Header.zSize-1);
   
    
    % Removes eof marker if needed
    if( Data(end) == 10 )
        Data(end) = [];
    
	% not sure why this is needed, but sometimes the data has one to few pixels
    elseif( length(Data) < Header.xSize*Header.ySize*Header.zSize )
        Data(end+1) = 0;
    end
    
    % format the data dimensions x,y,z for matlab use
    if Components == 2
        Data1 = permute(reshape(Data(1:2:end),Header.xSize,Header.ySize,Header.zSize,1),[2 1 3 4]);
        Data2 = permute(reshape(Data(2:2:end),Header.xSize,Header.ySize,Header.zSize,1),[2 1 3 4]);
        Data = cat(4,Data1,Data2);
    elseif Components == 3
        Data1 = permute(reshape(Data(1:3:end),Header.xSize,Header.ySize,Header.zSize,1),[2 1 3 4]);
        Data2 = permute(reshape(Data(2:3:end),Header.xSize,Header.ySize,Header.zSize,1),[2 1 3 4]);
        Data3 = permute(reshape(Data(3:3:end),Header.xSize,Header.ySize,Header.zSize,1),[2 1 3 4]);
        Data = cat(4,Data1,Data2,Data3);
    else
        Data = permute(reshape(Data,Header.xSize,Header.ySize,Header.zSize),[2 1 3]);    
    end
return


function M = dunzip(Z, DataType)
% DUNZIP - decompress DZIP output to recover original data
%
% USAGE:
% M = dzip(Z)
%
% VARIABLES:
% Z = compressed variable to decompress
% M = decompressed output
%
% NOTES: (1) The input variable Z is created by the DZIP function and
%            is a vector of type uint8
%        (2) The decompressed output will have the same data type and
%            dimensions as the original data provided to DZIP.
%        (3) See DZIP for other notes.
%        (4) Carefully tested, but no warranty; use at your own risk.
%        (5) Michael Kleder, Nov 2005

    import com.mathworks.mlwidgets.io.InterruptibleStreamCopier

    
    a=java.io.ByteArrayInputStream(int8(Z));
    b=java.util.zip.InflaterInputStream(a);    
    isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
    c = java.io.ByteArrayOutputStream;
    isc.copyStream(b,c);   
    
    M = typecast(c.toByteArray,DataType);
return


% A much faster way to RLE decode
function NewData = RLEDecodeEx(Data, DataSize, DataType)
    
    eval(['NewData = ',DataType,'(zeros(DataSize,1));']);
    eval(['Data = ',DataType,'(Data);']);
    index = 1;
    DataOffset = 1;
    FlagOffset = 1;
    ticker = 1;
    
    % Find the position of all the flags specifying to start of an individual pixels segement
    FlagInd = find(Data > 127 );
    
    % Loop through all the data chunks which should be organised as (RLE | individual pixels)
    for i=1:length(FlagInd)    
        
        % Read in all of the RLE Data up to the first individual pixel segement     
        eval(['PosData = ',DataType,'(Data(FlagOffset:FlagInd(i)-1));']);  
        
        if( FlagInd(i) ~= 1 )
            Lens = PosData(1:2:end);
            Vals = PosData(2:2:end)';
            clear PosData
            
            MaxLens = max(Lens(:));            
            nData = repmat(Vals,[MaxLens 1]); 
            
            clear Vals;
            
            
            ind = uint8(ones(size(nData)));

            for j=1:length(Lens)
                 if( Lens(j) < MaxLens )
                     ind(Lens(j)+1:MaxLens,j) = 0;
                 end
            end
            
            
            % Copy data 
            Num = sum(Lens);
            NewData(DataOffset:DataOffset+Num-1) = nData(find(ind==1));
            DataOffset = DataOffset+Num;
            clear Lens nData MaxLens ind Num
        end
        
        
        % Now read in the individual pixels
        Num = double(Data(FlagInd(i)) - 128);
        NewData(DataOffset:DataOffset+Num-1) = Data(FlagInd(i)+1:FlagInd(i)+Num);
        
        DataOffset = DataOffset+Num;        
        FlagOffset = FlagInd(i)+Num+1;
        
        
    end
    
    % Load the rest of the data in,
    if( FlagOffset ~= DataSize)
        eval(['PosData = ',DataType,'(Data(FlagOffset:end));']);  
        Lens = PosData(1:2:end);
        Vals = PosData(2:2:end)';
        clear PosData

        MaxLens = max(Lens(:));            
        nData = repmat(Vals,[MaxLens 1]); 

        clear Vals;


        ind = uint8(ones(size(nData)));

        for j=1:length(Lens)
             if( Lens(j) < MaxLens )
                 ind(Lens(j)+1:MaxLens,j) = 0;
             end
        end


        % Copy data 
        Num = sum(Lens);
        NewData(DataOffset:DataOffset+Num-1) = nData(find(ind==1));
        DataOffset = DataOffset+Num;
        clear Lens nData MaxLens ind Num
    end
    
             
return 
        
% Very Slow RLE Decode
function NewData = RLEDecode(Data,DataSize,DataType)
    
    eval(['NewData = ',DataType,'(zeros(DataSize,1));']);
    index = 1;
    offset = 1;
    ticker = 1;
    while(index < length(Data) )
        
        if( mod(ticker,10000) == 0 || index == 1 )
            disp(sprintf('Percentage Complete: %3.2f',100.*index./length(Data)));
        end
        
        Num = Data(index);
        
        if( Num < 0 )
            Num = Num+128;
            eval(['NewData(offset:offset+Num-1) = ',DataType,'(Data(index:index+Num-1));']);
            index = index + Num+1;
        else
            eval(['NewData(offset:offset+Num-1) = ',DataType,'(repmat(Data(index+1),[1,Num]));']);
            index = index + 2;
        end
        offset = offset+Num;
        
        ticker = ticker+1;
         
        
    end
    
    disp(num2str(offset));
%     ind1 = [2:2:length(Data)];
%     ind = find(Data(ind1)>0);
%     
%     lens = Data(ind1(ind)-1);
%     newData(
%     for i=1:2:length(Data)
%         
%         newData = [newData repmat(Data(i+1),[1,Data(i)])];
%     end
%     
        