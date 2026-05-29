function [XYAligned, XYZAligned] = ProcessZStackTimeCourse(SourceDirectory, TargetDirectory,channel_num,n_channels)

if isfolder(SourceDirectory)
    animal = regexp(SourceDirectory, '[A-Z]{2,3}0*\d+', 'match'); animal = animal{1};

    allTCfiles = fastdir(SourceDirectory);
elseif isfile(SourceDirectory)
    fparts = regexp(SourceDirectory, filesep, 'split');
    SourceDirectory = fullfile(fparts{1:end-1});
    allTCfiles{1} = fparts{end};
end

if(~isdir(TargetDirectory))
    mkdir(TargetDirectory);
end

cd(TargetDirectory)

accumulated_slice_count = 1;
leftover = false;

for i = 1:length(allTCfiles)
    targetfn = [allTCfiles{i}(1:end-4),'.tif'];
    info = imfinfo(fullfile(SourceDirectory, targetfn));

    if(info(1).ImageDescription(1)=='f')
        %ScanImage 5
        SI = assignments2StructOrObj(info(1).Software);

        loggingFramesPerFile = SI.hScan2D.logFramesPerFile;
        channelsSave = SI.hChannels.channelSave;
        FramesPerSlice = SI.hStackManager.framesPerSlice;
        stackNumSlices = SI.hStackManager.numSlices;
        fastZEnable = SI.hFastZ.enable;
        if(fastZEnable)
            FramesPerSlice = SI.hFastZ.numVolumes;
            stackNumSlices  = SI.hFastZ.numFramesPerVolume;
        end % SI5 can do slow zstack with fastz, but this script is not compatible.
        
        clear SI
    else
        %ScanImage 4
        tmp = assignments2StructOrObj(info(1).ImageDescription);

        loggingFramesPerFile = tmp.SI4.loggingFramesPerFile;
        channelsSave = tmp.SI4.channelsSave;
        FramesPerSlice = tmp.SI4.acqNumFrames;
        stackNumSlices = tmp.SI4.stackNumSlices;
        fastZEnable = tmp.SI4.fastZEnable;
        if(fastZEnable)
            FramesPerSlice = tmp.SI4.fastZNumVolumes;
        end
        clear tmp;
    end
    
   

    numChannels = numel(channelsSave);
    totalFrames = FramesPerSlice * stackNumSlices;
    numFiles = ceil(totalFrames / loggingFramesPerFile);
    
    im{i} = read_tiff(fullfile(SourceDirectory, targetfn),channel_num,n_channels);
    
    numFrames = size(im{i},3);
    if leftover
        z = [zeros(1,remaining_frames_in_slice),floor((0:(numFrames-remaining_frames_in_slice-1))/FramesPerSlice)+1]+1;
    else
        z = floor((0:numFrames-1)/FramesPerSlice)+1;
    end
    
    for slice = 1:max(z)
        fprintf('%d/%d\n',accumulated_slice_count,stackNumSlices);
        im_slice = im{i}(:,:,z == slice);
        im_slice(isnan(im_slice(:)))= 0;
        if leftover
            im_slice = cat(3,im_slice,leftover_frames);
        end
        if ~leftover && sum(z==slice)<FramesPerSlice && i < length(allTCfiles)
            leftover = true;
            leftover_frames = im_slice;
            remaining_frames_in_slice = FramesPerSlice-size(leftover_frames,3);
            break
        else
            leftover = false;
        end
        [StableSliceAverage(:,:,accumulated_slice_count), t] = make_stable_slice_average(im_slice(:,:,:));
%         for ch=1:numChannels
%             if(ch == ch_align)
%                 continue;
%             end
            im_ch = zeros(size(im_slice));
            for j=1:size(im_slice,3)
                im_ch(:,:,j)=BilinearImageRegistrator.shift(im_slice(:,:,j),t(j,:));
            end
            correctedZStack(:,:,accumulated_slice_count) = mean(im_ch,3);
%             figure; subplot(1,2,1); imagesc(StableSliceAverage(:,:,slice)); subplot(1,2,2); imagesc(correctedZStack(:,:,slice))
%         end
        accumulated_slice_count = accumulated_slice_count+1;
    end
end

correctedZStack(isnan(correctedZStack(:))) = 0;
disp('Aligning slices');
center_z_index = floor((stackNumSlices+1)/2);
cumulative_t = cell(stackNumSlices,1);
aligned = zeros(size(correctedZStack));
aligned(:,:,center_z_index ) = correctedZStack(:,:,center_z_index);
for direction = [-1 1]
    if(direction>0)
        max_dz=stackNumSlices-center_z_index;
    else
        max_dz=center_z_index-1;
    end
    cumulative_t{center_z_index}=[0 0];
    for dz = 1:max_dz
        tz = cvMotionCorrect...
            (correctedZStack(:,:,center_z_index+direction*dz),...
            correctedZStack(:,:,center_z_index+direction*(dz-1)));
        cumulative_t{center_z_index+direction*dz}...
            = cumulative_t{center_z_index+direction*(dz-1)} + tz(:,1:2);
        ind = center_z_index+direction*dz;
        %             for ch = 1:numChannels
        aligned(:,:,ind) = BilinearImageRegistrator.shift(correctedZStack(:,:,ind),cumulative_t{ind});
        %             end
    end
end
%     OriginalImages = int16(im{i});
    
StableAverage = int16(StableSliceAverage);
XYAligned = int16(correctedZStack);
XYZAligned = int16(aligned);

write_tiff(fullfile(TargetDirectory,'StableAverage.tif'),StableAverage);
write_tiff(fullfile(TargetDirectory,'XYAligned.tif'),XYAligned);
write_tiff(fullfile(TargetDirectory,'XYZAligned.tif'), XYZAligned);
