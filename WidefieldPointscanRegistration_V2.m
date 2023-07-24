%% version 2 intended to build a batch input, automatic registration using intensity-based or
% 2021/10/20 created, no modification yet
% 2022/03/23 first modified for 3 group SRS-FISH-Drug including sub folder
% search
% 2022/09/26 check it is already modified compatible with two colors
% conditions

%% image file input

destination_whole = 'F:\Data\Experiment Local\20230113 SRS-FISH-drug_drug distribution amend\20230113 Tiff_SRS-FISH-drug_drug distribution amend\';
FileName = '';
Key = '';
destination = append(destination_whole, FileName);
cd(destination)
SearchSubFolder_flag = 1;

if SearchSubFolder_flag == 1
    file_keyword = [Key,'*'];  % folder keyword
    files = dir(file_keyword);
    dirFlags = [files.isdir];
    subFolders = files(dirFlags);
    subFolders_name = strings([length(subFolders), 1]);
    for i = 1:length(subFolders)
        fprintf('Subfolders #%d = %s: found\n',i, subFolders(i).name);
        subFolders_name(i,1) = subFolders(i).name;
    end
    
    Process = input('Register? (1 or 0)\n');
    if Process ~= 1
        error('Abort processing\n');
    end
    
%     keyword = '*.tif';    % txt file keyword
    
    for i = 4:5%(length(subFolders))
%     for i = 3
        cd(append(destination,'/',subFolders_name(i)));
        No_FOV_key = ['*_FOV','*CY*.tif'];
        listing = dir(No_FOV_key);
        listing = listing(~startsWith({listing.name}, 'Registered'));

        name = string(extractfield(listing, 'name'));
        No_FOV = length(name);
        
        for FOV = 4:No_FOV+3
%         for FOV = 4
            PointScan_file = ['Step size*FOV',num2str(FOV),'_CH.txt.tif'];
            listing = dir(PointScan_file);
            name = string(extractfield(listing, 'name'));
            PointScan_file = name;
            fprintf('Image set:\n');
            name_seq = [[1:length(name)]', name'];
            disp(name_seq);
          
             % image preprocessing for visualization
            WidefieldT_file = [Key,'*_FOV',num2str(FOV),'*TRANS.tif'];
            listing_T = dir(WidefieldT_file);
            if ~isempty(listing_T)
            name = string(extractfield(listing_T, 'name'));
            WidefieldT_file = name;
            end
            
            Widefield_file = [Key,'*_FOV',num2str(FOV),'*CY*.tif'];
            listing = dir(Widefield_file);
            name = string(extractfield(listing, 'name'));
            Widefield_file = name;
            
            Widefield2_file = [Key,'*_FOV',num2str(FOV),'*CY5X.tif'];
            listing_2 = dir(Widefield2_file);
            if ~isempty(listing_2)
            name = string(extractfield(listing_2, 'name'));
            Widefield2_file = name;
            end
            
            
            
            PS = imread(append('./', PointScan_file));
            PS_n = (PS-min(PS,[],'all'))/(max(PS,[],'all')-min(PS,[],'all'));
            
            
            if ~isempty(listing_T)
            WFt = imread(append('./', WidefieldT_file));
            WFt = medfilt2(single(WFt(:,:,2)),[3,3]);  % only green channel
%             WFt = flip(medfilt2(single(WFt(:,:,2)),[3,3]),1);  % only green channel
            WFt_n = (WFt-min(WFt,[],'all'))/(max(WFt,[],'all')-min(WFt,[],'all'));
            end
            
            WF = imread(append('./', Widefield_file));
            WF = medfilt2(single(WF(:,:,1)),[4,4]); % only green channel
%             WF = flip(medfilt2(single(WF(:,:,2)),[4,4]),1); % only green channel
            
            if ~isempty(listing_2)
            WF2 = imread(append('./', Widefield2_file));
            WF2 = medfilt2(single(WF2(:,:,1)),[3,3]);  % only red channel
%             WF2 = flip(medfilt2(single(WF2(:,:,1)),[3,3]),1);  % only red channel
            WF2_n = (WF2-min(WF2,[],'all'))/(max(WF2,[],'all')-min(WF2,[],'all'));
            end
            
            if ~isempty(listing_T)
            figure(1);
            subplot(1,2,1);
            imagesc(WFt_n)
            subplot(1,2,2);
            imagesc(PS_n)
            title('Raw');
            end
            
            %% select aligned point
%             clear movingPoints
%             clear fixedPoints
            if exist('movingPoints','var') == 0
                cpselect(WFt_n,PS_n);    % save movingpoints and fixed points
                % save selected point to workspace file
                RegistrationPoints_file = 'RegistrationPoints_Proj.mat';
                m = matfile(RegistrationPoints_file,'Writable',true); % load as handle of .mat
                % tmp = load(RegistrationPoints_file); % load as structure format
                % % rmmatvar(PureChemicalMap_file, 'MICROBIOME'); % clear .mat function is at
                % % the bottom
                movingPointsAdjusted = cpcorr(movingPoints,fixedPoints,WFt_n,PS_n);
                movingPointsAdjusted = movingPoints;
                movingPoints = movingPointsAdjusted;
                m.movingPoints = movingPoints;
                m.fixedPoints = fixedPoints;
            end
            
            %% adjust image
            movingPointsAdjusted = cpcorr(movingPoints,fixedPoints,WFt_n,PS_n); % weird things happens 
            movingPointsAdjusted = movingPoints;
            tform = estimateGeometricTransform(movingPointsAdjusted,fixedPoints,'projective');%,'NonreflectiveSimilarity','pwl'
            %     tform = fitgeotrans(movingPoints(1:3,:),fixedPoints(1:3,:),'affine');
            
            if ~isempty(listing_T)
            WFt_registered = imwarp(WFt,tform,'OutputView',imref2d(size(PS)));
            end
            
            WF_registered = imwarp(WF,tform,'OutputView',imref2d(size(PS)));
            
            if ~isempty(listing_2)
            WF2_registered = imwarp(WF2,tform,'OutputView',imref2d(size(PS)));
            end
            
%             % extra sift registration: run sift first
%             WFt_registered=imtransform(WFt_registered,tform_sift, 'XData',[1 size(PS,2)], 'YData',[1 size(PS,1)]);
%             WF_registered=imtransform(WF_registered,tform_sift, 'XData',[1 size(PS,2)], 'YData',[1 size(PS,1)]);
            
            if ~isempty(listing_T)
            figure(2);
            subplot(1,2,1);
            imagesc(WFt_registered);
            subplot(1,2,2);
            imagesc(PS);
            title('Registered_trans');
            end
            
            figure(3);
            subplot(1,2,1);
            imagesc(WF_registered);
            subplot(1,2,2);
            imagesc(PS);
            title('Registered_Cy3');
            
            %% export image
            
            % saveas(WF_registered, [FilePath,'\FOV',FOV,'_Registered.tif']); % need
            % handle instead of matrix
            
            % imwrite(WF_registered, [FilePath,'\FOV',FOV,'_Registered.tif']);
            
            % WF_registered = imread('example.tif');
%             t = Tiff(['FOV',FOV,'_Registered.tif'],'w');
%             t.setTag('ImageLength',size(WF_registered,1));
%             t.setTag('ImageWidth', size(WF_registered,2));
%             t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
%             t.setTag('BitsPerSample', 32);
%             t.setTag('SampleFormat', Tiff.SampleFormat.IEEEFP);
%             t.setTag('SamplesPerPixel', size(WF_registered,3));
%             % t.setTag('TileWidth', 128);
%             % t.setTag('TileLength', 128);
%             t.setTag('Compression', Tiff.Compression.None);
%             t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
%             t.setTag('Software', 'MATLAB');
%             t.write(WF_registered);
%             t.close();
       
%             save_tiff_modi(WFt_registered, ['RegisteredFOV',name,num2str(FOV),'_Trans.tif']);
%             save_tiff_modi(WF_registered, ['RegisteredFOV',name,num2str(FOV),'_Cy3.tif']);
            
        if ~isempty(listing_T)
            save_tiff_modi(WFt_registered, strcat('Registered_',WidefieldT_file));
        end
            save_tiff_modi(WF_registered, strcat('Registered_',Widefield_file));
            
        if ~isempty(listing_2)
            save_tiff_modi(WF2_registered, strcat('Registered_',Widefield2_file));
        end
            
        end
        
        cd ..
    end
end


%% other functions

function save_tiff_modi(img_tiff,outputName,new_path)
% later saved as matrix2tif_2D in the code management
if nargin < 3
    new_path = cd;
end

old_path = cd;

if (exist(new_path,'dir')==0)
    mkdir(new_path);
end
cd(new_path);

t = Tiff(outputName, 'w');
%
tagstruct.ImageLength = size(img_tiff, 1);
tagstruct.ImageWidth = size(img_tiff, 1);
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 32;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
t.setTag(tagstruct);
t.write(img_tiff);
t.close();
cd(old_path);
% B=imread(outputName)

end