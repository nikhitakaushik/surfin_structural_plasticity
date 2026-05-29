mouse_id = "JL122";
date_list = ["250403, 250404, 250405, 250406, 250407, 250408, 250409, 250410, 250411, 250412, 250413, 250414, 250415"];
fov_list = ["FOV1","FOV2","FOV3","FOV4","FOV5","FOV6","FOV7","FOV8","FOV9","FOV10"];
align_ch = [1,1,1,1,1,3,3,3,2,1];
n_channels = 2;
%%
mouse_id = "JL123";
date_list = ["250403, 250404, 250405, 250406, 250407, 250408, 250409, 250410, 250411, 250412, 250413, 250414, 250415"];
fov_list = ["FOV1","FOV2","FOV3","FOV4","FOV5","FOV6","FOV7","FOV8"];
align_ch = [2,2,3,3,1,1,1,3];
n_channels = 2;

%%
%mouse_id = "JL117";
%date_list = ["250224","250225"];
%fov_list = ["FOV1","FOV2","FOV3","FOV4","FOV5","FOV6","FOV7"];
%align_ch = [1,1,2,2,1,1,2];
%n_channels = 2;

%%

for date_id = date_list
    for fov_i = 1:length(fov_list)
        fov_id = fov_list(fov_i);
    
        if align_ch(fov_i) == 1
            color_name = ["GreenCh"];
        elseif align_ch(fov_i) == 2
            color_name = ["RedCh"];
        elseif align_ch(fov_i) == 3
            color_name = ["GreenCh","RedCh"];
        end

        channel_id.GreenCh = 1;
        channel_id.RedCh = 2;

        for ch_name = color_name
            try
                if contains(fov_id,"_")
                    fov_id = strjoin(strsplit(fov_id,"_"),"\");
                end
                ProcessZStackTimeCourse(fullfile('Z:\Data\ImagingRig3\',date_id,mouse_id,date_id,fov_id),fullfile('Z:\People\JenLi\Imaging_Data\L23_L5_structural_imaging\motion_corrected\',mouse_id,date_id,fov_id,ch_name,'\'),channel_id.(ch_name),n_channels);
                
            catch
                disp([fov_id,'cant find'])
            end
        end
    end
end


% ProcessZStackTimeCourse(fullfile('Z:\Data\Bruker\',date_id,mouse_id,date_id,fov_id),fullfile('Z:\People\JenLi\Imaging_Data\',mouse_id,date_id,fov_id,color_name,'\'),align_ch,n_channels);