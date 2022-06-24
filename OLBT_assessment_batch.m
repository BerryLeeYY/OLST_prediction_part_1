clear
addpath 'D:\files\Thesis_first_part_data'
addpath 'C:\Users\a1003\OneDrive\桌面\Thesis'

filenames = dir('D:\files\Thesis_first_part_data');

for n = 1:length(filenames)
    try
        subjfile = filenames(n).name;

        filename = load(filenames(n).name);
        name = filename.(subsref(fieldnames(filename),substruct('{}',{1})));
        label = filename.(subsref(fieldnames(filename),substruct('{}',{1}))).Trajectories.Labeled.Labels;
        path = filename.(subsref(fieldnames(filename),substruct('{}',{1}))).Trajectories.Labeled.Labels;
        split_filename = split(subjfile, ["_", "-"]);
        if ismember(["eyecl"], split_filename) || ismember(["eyescl"], split_filename)
            trial = "EC";
        else
            trial = "EP";
        end
        sub_ID = "sub" + "_" + subjfile(4:5) + "_" + trial + "_" + "W2";
        
        FP1_COP_data = name.Force(1).COP;
        FP2_COP_data = name.Force(2).COP;
        FP3_COP_data = name.Force(3).COP;
        FP4_COP_data = name.Force(4).COP;
        FP5_COP_data = name.Force(5).COP;
        FP6_COP_data = name.Force(6).COP;
        FP7_COP_data = name.Force(7).COP;

        FP1_Force_data = name.Force(1).Force;
        FP2_Force_data = name.Force(2).Force;
        FP3_Force_data = name.Force(3).Force;
        FP4_Force_data = name.Force(4).Force;
        FP5_Force_data = name.Force(5).Force;
        FP6_Force_data = name.Force(6).Force;
        FP7_Force_data = name.Force(7).Force;

        data_len = length(FP1_COP_data);
        frq = length(FP1_COP_data) / length(name.Trajectories.Labeled.Data(26,1,:));
        time = data_len / frq;
        %% 1.	Lifting forefoot or heel 
        %%%	variability of the signal %%%
        %%% mean of the signal %%%
        %%% maximum of the signal %%%

        %%%% R(L)DM1 position in z axis %%%%
        if sum(strcmp( path, 'RDM1')) == 0
            RDM1_position = find(strcmp( path, 'RDMT1'));
        else
            RDM1_position = find(strcmp( path, 'RDM1'));
        end
        RDM1_height = name.Trajectories.Labeled.Data(RDM1_position,3,:);
        RDM1_height = reshape(RDM1_height, [1, time])/10;

        if sum(strcmp( path, 'LDM1')) == 0
            LDM1_position = find(strcmp( path, 'LDMT1'));
        else
            LDM1_position = find(strcmp( path, 'LDM1'));
        end
        LDM1_height = name.Trajectories.Labeled.Data(LDM1_position,3,:);
        LDM1_height = reshape(LDM1_height, [1, time])/10;

        %%%% R(L)CAL1 position in z axis %%%%
        RCAL1_position = find(strcmp( path, 'RCAL1'));
        RCAL1_height = name.Trajectories.Labeled.Data(RCAL1_position,3,:);
        RCAL1_height = reshape(RCAL1_height, [1, time])/10;
        LCAL1_position = find(strcmp( path, 'LCAL1'));
        LCAL1_height = name.Trajectories.Labeled.Data(LCAL1_position,3,:);
        LCAL1_height = reshape(LCAL1_height, [1, time])/10;

        %{
        figure
        plot(RDM1_height)
        hold on
        plot(RCAL1_height)
        plot(LDM1_height)
        plot(LCAL1_height)
        title('Lifting forefoot or heel')
        legend('RDM1', 'RCAL1', 'LDM1', 'LCAL1')
        %}

        %%%	variability of the signal %%%
        %%% mean of the signal %%%
        %%% maximum of the signal %%%
        if std(RDM1_height) > std(LDM1_height)
            sorted_hei_var = [nanstd(LDM1_height), nanstd(LCAL1_height), nanstd(RDM1_height), nanstd(RCAL1_height)];
            sorted_hei_mean = [nanmean(LDM1_height), nanmean(LCAL1_height), nanmean(RDM1_height), nanmean(RCAL1_height)];
            sorted_hei_max = [nanmax(LDM1_height), nanmax(LCAL1_height), nanmax(RDM1_height), nanmax(RCAL1_height)];
        else
            sorted_hei_var = [nanstd(RDM1_height), nanstd(RCAL1_height), nanstd(LDM1_height), nanstd(LCAL1_height)];
            sorted_hei_mean = [nanmean(RDM1_height), nanmean(RCAL1_height), nanmean(LDM1_height), nanmean(LCAL1_height)];
            sorted_hei_max = [nanmax(RDM1_height), nanmax(RCAL1_height), nanmax(LDM1_height), nanmax(LCAL1_height)];
        end
        features_hei = [sorted_hei_var, sorted_hei_mean, sorted_hei_max];
        hei_col = ["SL_fore_var_height", "SL_heel_var_height", "OL_fore_var_height", "OL_heel_var_height",... 
            "SL_fore_mean_height", "SL_heel_mean_height", "OL_fore_mean_height", "OL_heel_mean_height",...
            "SL_fore_max_height", "SL_heel_max_height", "OL_fore_max_height", "OL_heel_max_height"];
        hei_df = [hei_col; features_hei];

        foot_height_count = 0;
        window = 5;
        if abs(nanmean(RDM1_height)) > abs(nanmean(LDM1_height)) %%% left foot on the ground
            landing_foot_height = LDM1_height;
            heel_height = LCAL1_height;
            forefoot_threshold = nanmean(landing_foot_height) + 4.3*nanstd(landing_foot_height);
            heel_threshold = nanmean(heel_height) + 5*nanstd(heel_height);
            for i = 1:(time-1)
                try
                    if landing_foot_height(i) < forefoot_threshold && landing_foot_height(i+1) > forefoot_threshold && landing_foot_height(i+window) > forefoot_threshold || (heel_height(i) < heel_threshold && heel_height(i+1) > heel_threshold && heel_height(i+window) > heel_threshold)
                        foot_height_count = foot_height_count + 1;
                    end
                catch
                    continue
                end
            end
        elseif abs(nanmean(RDM1_height)) < abs(nanmean(LDM1_height))%%% right foot on the ground
            landing_foot_height = RDM1_height;
            heel_height = RCAL1_height;
            forefoot_threshold = nanmean(landing_foot_height) + 4.3*nanstd(landing_foot_height);
            heel_threshold = nanmean(heel_height) + 5*nanstd(heel_height);
            for i = 1:(time-1)
                try
                    if landing_foot_height(i) < forefoot_threshold && landing_foot_height(i+1) > forefoot_threshold && landing_foot_height(i+window) > forefoot_threshold || (heel_height(i) < heel_threshold && heel_height(i+1) > heel_threshold && heel_height(i+window) > heel_threshold)
                        foot_height_count = foot_height_count + 1;
                    end
                catch
                    continue
                end
            end
        end

        %% 2.1 (felxion)Moving hip into more than 30 degrees of flexion or adduction/abduction (reference)
        %%%	Angle in sagittal plane between RFTC-RTH vector and STRN-MidASI %%%
        %%% Angle in sagittal plane between LFTC-LTH vector and STRN-MidASI %%%
        %%%	Angle in frontal plane between RFTC-RTH vector and STRN-MidASI %%%
        %%% Angle in frontal plane between LFTC-LTH vector and STRN-MidASI %%%
        %%%	variability of the signal %%%
        %%% mean of the signal %%%
        %%% maximum of the signal %%%
        %%% Duration of the degree more than 30 degrees %%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%% 2.1 flexion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% R(L)TH position in xyz axis %%%%  =====> replace R(L)TH to R(L)FLE
        RFLE_position = find(strcmp( path, 'RFLE'));
        RFLE_xyz_position = name.Trajectories.Labeled.Data(RFLE_position,1:3,:);
        RFLE_xyz_position = reshape(RFLE_xyz_position, [3, time])/10;
        LFLE_position = find(strcmp( path, 'LFLE'));
        LFLE_xyz_position = name.Trajectories.Labeled.Data(LFLE_position,1:3,:);
        LFLE_xyz_position = reshape(LFLE_xyz_position, [3, time])/10;

        %%%% R(L)FTC position in xyz axis %%%%
        RFTC_position = find(strcmp( path, 'RFTC'));
        RFTC_xyz_position = name.Trajectories.Labeled.Data(RFTC_position,1:3,:);
        RFTC_xyz_position = reshape(RFTC_xyz_position, [3, time])/10;
        LFTC_position = find(strcmp( path, 'LFTC'));
        LFTC_xyz_position = name.Trajectories.Labeled.Data(LFTC_position,1:3,:);
        LFTC_xyz_position = reshape(LFTC_xyz_position, [3, time])/10;

        %%%% R(L)ASI position in xyz axis %%%%
        RASI_position = find(strcmp( path, 'RASI'));
        RASI_xyz_position = name.Trajectories.Labeled.Data(RASI_position,1:3,:);
        RASI_xyz_position = reshape(RASI_xyz_position, [3, time])/10;
        LASI_position = find(strcmp( path, 'LASI'));
        LASI_xyz_position = name.Trajectories.Labeled.Data(LASI_position,1:3,:);
        LASI_xyz_position = reshape(LASI_xyz_position, [3, time])/10;
        MidASI_xyz_position = (RASI_xyz_position + LASI_xyz_position)/2;
        %%%% STRN position in xyz axis %%%%
        if sum(strcmp( path, 'STRN')) == 0
            STRN_position = find(strcmp( path, 'STERNSUP'));
        else
            STRN_position = find(strcmp( path, 'STRN'));
        end
        STRN_xyz_position = name.Trajectories.Labeled.Data(STRN_position,1:3,:);
        STRN_xyz_position = reshape(STRN_xyz_position, [3, time])/10;

        %%%% Vectors RFTC-RTH, LFTC-LTH, STRN-MidASI%%%%
        V_RFTC_RFLE_all = RFLE_xyz_position(:, :) - RFTC_xyz_position(:, :);
        V_LFTC_LFLE_all = LFLE_xyz_position(:, :) - LFTC_xyz_position(:, :);
        V_STRN_MidASI_all = MidASI_xyz_position(:, :) - STRN_xyz_position(:, :);

        %%%% angle between vectors => cos(angle) = dot(A,B) / (norm(A).*norm(B))
        r_a = V_RFTC_RFLE_all;
        r_b = V_STRN_MidASI_all;
        for i = 1:length(r_a)
            R_cos(i)= dot(r_a(:,i),r_b(:,i)) / (norm(r_a(:,i), 2)*norm(r_b(:,i), 2));
        end
        R_angle_all = acosd(R_cos);

        l_a = V_LFTC_LFLE_all;
        l_b = V_STRN_MidASI_all;
        for i = 1:length(l_a)
            L_cos(i)= dot(l_a(:,i),l_b(:,i)) / (norm(l_a(:,i), 2)*norm(l_b(:,i), 2));
        end
        L_angle_all = acosd(L_cos);

        %%%% Vectors RFTC-RTH, LFTC-LTH, STRN-MidASI%%%%
        V_RFTC_RFLE = RFLE_xyz_position(1:2:3, :) - RFTC_xyz_position(1:2:3, :);
        V_LFTC_LFLE = LFLE_xyz_position(1:2:3, :) - LFTC_xyz_position(1:2:3, :);
        V_STRN_MidASI = MidASI_xyz_position(1:2:3, :) - STRN_xyz_position(1:2:3, :);

        %%%% angle between vectors => cos(angle) = dot(A,B) / (norm(A).*norm(B))
        r_a = sqrt(V_RFTC_RFLE(1,:).*V_RFTC_RFLE(1,:) + V_RFTC_RFLE(2,:).*V_RFTC_RFLE(2,:));
        r_b = sqrt(V_STRN_MidASI(1,:).*V_STRN_MidASI(1,:) + V_STRN_MidASI(2,:).*V_STRN_MidASI(2,:));
        r = abs(r_a).*abs(r_b);
        dot_r = V_RFTC_RFLE(1,:).*V_STRN_MidASI(1,:) + V_RFTC_RFLE(2,:).*V_STRN_MidASI(2,:);
        R_cos = dot_r ./ r;
        R_angle = acosd(R_cos);
        l_a = sqrt(V_LFTC_LFLE(1,:).*V_LFTC_LFLE(1,:) + V_LFTC_LFLE(2,:).*V_LFTC_LFLE(2,:));
        l_b = sqrt(V_STRN_MidASI(1,:).*V_STRN_MidASI(1,:) + V_STRN_MidASI(2,:).*V_STRN_MidASI(2,:));
        l = abs(l_a).*abs(l_b);
        dot_l = V_LFTC_LFLE(1,:).*V_STRN_MidASI(1,:) + V_LFTC_LFLE(2,:).*V_STRN_MidASI(2,:);
        L_cos = dot_l ./ l;
        L_angle = acosd(L_cos);

        %{
        figure
        plot(R_angle)
        hold on
        plot(L_angle)
        title('Hip flexion')
        legend('right', 'left')
        %}

        %%%	variability of the signal %%%
        %%% mean of the signal %%%
        %%% maximum of the signal %%%
        %%% Duration of the degree more than 30 degrees %%%
        R_flex_duration = 0;
        for i = 1:length(R_angle)
            if R_angle(1,i) >= 30
                R_flex_duration = R_flex_duration + 1;
            end
        end
        L_flex_duration = 0;
        for i = 1:length(L_angle)
            if L_angle(1,i) >= 30
                L_flex_duration = L_flex_duration + 1;
            end
        end

        if mean(RDM1_height) > mean(LDM1_height) %%%%% left limb is supporting limb
            features_flex_angle = [nanstd(L_angle), nanstd(R_angle), nanmean(L_angle), nanmean(R_angle), nanmax(L_angle), nanmax(R_angle), L_flex_duration, R_flex_duration];
        else
            features_flex_angle = [nanstd(R_angle), nanstd(L_angle), nanmean(R_angle), nanmean(L_angle), nanmax(R_angle), nanmax(L_angle), R_flex_duration, L_flex_duration];
        end

        flex_ang_col = ["SL_var_flex", "OL_var_flex", "SL_mean_flex", "OL_mean_flex", ... 
            "SL_max_flex", "OL_max_flex", "SL_duration_flex", "OL_duration_flex"];

        flex_ang_df = [flex_ang_col; features_flex_angle];

        %% 2.2 (adduction/abduction) Moving hip into more than 30 degrees of flexion or adduction/abduction
        %%%	Angle in sagittal plane between RFTC-RTH vector and STRN-MidASI %%%
        %%% Angle in sagittal plane between LFTC-LTH vector and STRN-MidASI %%%
        %%%	Angle in frontal plane between RFTC-RTH vector and STRN-MidASI %%%
        %%% Angle in frontal plane between LFTC-LTH vector and STRN-MidASI %%%
        %%%	variability of the signal %%%
        %%% mean of the signal %%%
        %%% maximum of the signal %%%
        %%% Duration of the degree more than 30 degrees %%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%% 2.1 flexion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% R(L)TH position in xyz axis %%%%  =====> replace R(L)TH to R(L)FLE
        RFLE_position = find(strcmp( path, 'RFLE'));
        RFLE_xyz_position = name.Trajectories.Labeled.Data(RFLE_position,1:3,:);
        RFLE_xyz_position = reshape(RFLE_xyz_position, [3, time])/10;
        LFLE_position = find(strcmp( path, 'LFLE'));
        LFLE_xyz_position = name.Trajectories.Labeled.Data(LFLE_position,1:3,:);
        LFLE_xyz_position = reshape(LFLE_xyz_position, [3, time])/10;

        %%%% R(L)FTC position in xyz axis %%%%
        RFTC_position = find(strcmp( path, 'RFTC'));
        RFTC_xyz_position = name.Trajectories.Labeled.Data(RFTC_position,1:3,:);
        RFTC_xyz_position = reshape(RFTC_xyz_position, [3, time])/10;
        LFTC_position = find(strcmp( path, 'LFTC'));
        LFTC_xyz_position = name.Trajectories.Labeled.Data(LFTC_position,1:3,:);
        LFTC_xyz_position = reshape(LFTC_xyz_position, [3, time])/10;

        %%%% R(L)ASI position in xyz axis %%%%
        RASI_position = find(strcmp( path, 'RASI'));
        RASI_xyz_position = name.Trajectories.Labeled.Data(RASI_position,1:3,:);
        RASI_xyz_position = reshape(RASI_xyz_position, [3, time])/10;
        LASI_position = find(strcmp( path, 'LASI'));
        LASI_xyz_position = name.Trajectories.Labeled.Data(LASI_position,1:3,:);
        LASI_xyz_position = reshape(LASI_xyz_position, [3, time])/10;
        MidASI_xyz_position = (RASI_xyz_position + LASI_xyz_position)/2;
        %%%% STRN position in xyz axis %%%%
        if sum(strcmp( path, 'STRN')) == 0
            STRN_position = find(strcmp( path, 'STERNSUP'));
        else
            STRN_position = find(strcmp( path, 'STRN'));
        end
        STRN_xyz_position = name.Trajectories.Labeled.Data(STRN_position,1:3,:);
        STRN_xyz_position = reshape(STRN_xyz_position, [3, time])/10;

        %%%% Vectors RFTC-RTH, LFTC-LTH, STRN-MidASI%%%%
        V_RFTC_RFLE = RFLE_xyz_position(2:3, :) - RFTC_xyz_position(2:3, :);
        V_LFTC_LFLE = LFLE_xyz_position(2:3, :) - LFTC_xyz_position(2:3, :);
        V_STRN_MidASI = MidASI_xyz_position(2:3, :) - STRN_xyz_position(2:3, :);

        %%%% angle between vectors => cos(angle) = dot(A,B) / (norm(A).*norm(B))
        r_a = sqrt(V_RFTC_RFLE(1,:).*V_RFTC_RFLE(1,:) + V_RFTC_RFLE(2,:).*V_RFTC_RFLE(2,:));
        r_b = sqrt(V_STRN_MidASI(1,:).*V_STRN_MidASI(1,:) + V_STRN_MidASI(2,:).*V_STRN_MidASI(2,:));
        r = abs(r_a).*abs(r_b);
        dot_r = V_RFTC_RFLE(1,:).*V_STRN_MidASI(1,:) + V_RFTC_RFLE(2,:).*V_STRN_MidASI(2,:);
        R_cos = dot_r ./ r;
        R_abd_angle = acosd(R_cos);
        l_a = sqrt(V_LFTC_LFLE(1,:).*V_LFTC_LFLE(1,:) + V_LFTC_LFLE(2,:).*V_LFTC_LFLE(2,:));
        l_b = sqrt(V_STRN_MidASI(1,:).*V_STRN_MidASI(1,:) + V_STRN_MidASI(2,:).*V_STRN_MidASI(2,:));
        l = abs(l_a).*abs(l_b);
        dot_l = V_LFTC_LFLE(1,:).*V_STRN_MidASI(1,:) + V_LFTC_LFLE(2,:).*V_STRN_MidASI(2,:);
        L_cos = dot_l ./ l;
        L_abd_angle = acosd(L_cos);

        %{
        figure
        plot(R_abd_angle)
        hold on
        plot(L_abd_angle)
        title('Hip adduction/abduction')
        legend('right', 'left')
        %}

        %%%	variability of the signal %%%
        %%% mean of the signal %%%
        %%% maximum of the signal %%%
        %%% Duration of the degree more than 30 degrees %%%
        R_abd_duration = 0;
        for i = 1:length(R_abd_angle)
            if R_abd_angle(1,i) >= 30
                R_abd_duration = R_abd_duration + 1;
            end
        end
        L_abd_duration = 0;
        for i = 1:length(L_abd_angle)
            if L_abd_angle(1,i) >= 30
                L_abd_duration = L_abd_duration + 1;
            end
        end

        if mean(RDM1_height) > mean(LDM1_height) %%%%% left limb is supporting limb
            features_abd_angle = [nanstd(L_abd_angle), nanstd(R_abd_angle), nanmean(L_abd_angle), nanmean(R_abd_angle), nanmax(L_abd_angle), nanmax(R_abd_angle), L_abd_duration, R_abd_duration];
        else
            features_abd_angle = [nanstd(R_abd_angle), nanstd(L_abd_angle), nanmean(R_abd_angle), nanmean(L_abd_angle), nanmax(R_abd_angle), nanmax(L_abd_angle), R_abd_duration, L_abd_duration];
        end

        abd_ang_col = ["SL_var_abd", "OL_var_abd", "SL_mean_abd", "OL_mean_abd",...
            "SL_max_abd", "OL_max_abd", "SL_duration_abd", "OL_duration_abd"];

        abd_ang_df = [abd_ang_col; features_abd_angle];

        angle_count = 0;
        if abs(nanmean(R_angle_all)) > abs(nanmean(L_angle_all)) %%% right leg main moving leg
            for i = 1:(time-1)
                if abs(R_angle_all(i)- abs(nanmean(R_angle_all(1:1300)))) < 12 &&  abs(R_angle_all(i+1) - abs(nanmean(R_angle_all(1:1300)))) > 12
                    angle_count = angle_count + 1;
                end
            end
        elseif abs(nanmean(R_angle_all)) < abs(nanmean(L_angle_all)) %%% left leg main moving leg
            for i = 1:(time-1)
                if abs(L_angle_all(i)- abs(nanmean(L_angle_all(1:1300)))) < 12 &&  abs(L_angle_all(i+1)- abs(nanmean(L_angle_all(1:1300)))) > 12
                    angle_count = angle_count + 1;
                end
            end
        end



        %% 3.	Stepping, stumbling, or falling 
        %%% stepping frequency 
        %%% the total duration of stepping
        %%% Force platform signal %%%

        force_mean_all = [nanmean(0-FP1_Force_data(3,:)); nanmean(0-FP2_Force_data(3,:)); nanmean(0-FP3_Force_data(3,:)); nanmean(0-FP4_Force_data(3,:)); nanmean(0-FP5_Force_data(3,:)); nanmean(0-FP6_Force_data(3,:)); nanmean(0-FP7_Force_data(3,:))];
        force_std_all = [nanstd(0-FP1_Force_data(3,:)); nanstd(0-FP2_Force_data(3,:)); nanstd(0-FP3_Force_data(3,:)); nanstd(0-FP4_Force_data(3,:)); nanstd(0-FP5_Force_data(3,:)); nanstd(0-FP6_Force_data(3,:)); nanstd(0-FP7_Force_data(3,:))];

        force_all = [FP1_Force_data(3,:); FP2_Force_data(3,:); FP3_Force_data(3,:); FP4_Force_data(3,:); FP5_Force_data(3,:); FP6_Force_data(3,:); FP7_Force_data(3,:)];
        [pks, loc] = findpeaks(force_std_all, "MinPeakHeight", 1);

        if length(loc) == 2
            if nanmean(0-force_all(loc(1),:)) > nanmean(0-force_all(loc(2),:))
                supporting_FP = force_all(loc(1),:);
                swinging_FP = force_all(loc(2),:);
            elseif nanmean(0-force_all(loc(1),:)) < nanmean(0-force_all(loc(2),:))
                supporting_FP = force_all(loc(2),:);
                swinging_FP = force_all(loc(1),:);
            end
        elseif length(loc) == 1
            supporting_FP = zeros(size(FP1_Force_data(3,:),1), size(FP1_Force_data(3,:),2));
            swinging_FP = zeros(size(FP1_Force_data(3,:),1), size(FP1_Force_data(3,:),2));
        end


        %{ 
        figure
        plot(supporting_FP)
        hold on
        plot(swinging_FP)
        title('force plat signal')
        findpeaks(0-swinging_FP, 'MinPeakProminence',100, 'MinPeakHeight',100)
        %}

        %%% stepping frequency 
        if nanmean(swinging_FP) < 0
            swinging_FP = 0- swinging_FP;
        end
        stepping_freq = length(findpeaks(swinging_FP, 'MinPeakProminence',100, 'MinPeakHeight',100));
        %%% the total duration of stepping
        stepping_duration = 0;
        for i  = 1:length(swinging_FP)
            if abs(swinging_FP(1,i)) >= 100
                stepping_duration = stepping_duration + 1;
            end
        end

        features_stepping = [stepping_freq, stepping_duration];
        stepping_col = ["stepping_freq", "stepping_duration"];

        stepping_df = [stepping_col; features_stepping];

        for i =  1:length(swinging_FP)
            if abs(swinging_FP(i)) > 100
                filtered_swinging_FP(i) = swinging_FP(i);
            elseif abs(swinging_FP(i)) < 100
                filtered_swinging_FP(i) = 0;
            end
        end

        stepping_count = 0;
        stepping_index(1) = nan;
        [count_pks, count_loc] = findpeaks(swinging_FP, "MinPeakHeight", 20, "MinPeakProminence", 40, "MinPeakDistance", 1000);
        if length(loc) == 2
            stepping_count = length(count_pks);
            if stepping_count > 10
                [count_pks, count_loc] = findpeaks(swinging_FP, "MinPeakHeight", 20, "MinPeakProminence", 40, "MinPeakDistance", 1000);
                stepping_count = length(count_pks);
                if stepping_count > 10
                    [count_pks, count_loc] = findpeaks(swinging_FP, "MinPeakHeight", 100, "MinPeakProminence", 40, "MinPeakDistance", 5000);
                    stepping_count = length(count_pks);
                end
            end
        else
            stepping_count = 0;
        end




        %% 4.	Lifting hands off iliac crests 
        %%% Distance between R(L)WRL and R(L)ASI %%% ====> RWRL and LWRL
        %%% variability 
        %%% max - min

        %%%% R(L)WRL position in xyz axis %%%%
        if sum(strcmp( path, 'RWRL')) == 0
            RWRL_position = find(strcmp( path, 'RWRR'));
        else
            RWRL_position = find(strcmp( path, 'RWRL'));
        end
        RWRL_xyz_position = name.Trajectories.Labeled.Data(RWRL_position,1:3,:);
        RWRL_xyz_position = reshape(RWRL_xyz_position, [3, time])/10;

        if sum(strcmp( path, 'LWRL')) == 0
            LWRL_position = find(strcmp( path, 'LWRR'));
        else
            LWRL_position = find(strcmp( path, 'LWRL'));
        end
        LWRL_xyz_position = name.Trajectories.Labeled.Data(LWRL_position,1:3,:);
        LWRL_xyz_position = reshape(LWRL_xyz_position, [3, time])/10;

        %{
        %%%% R(L)ASI position in xyz axis %%%%
        RASI_position = find(strcmp( path, 'RASI'));
        RASI_xyz_position = name.Trajectories.Labeled.Data(RASI_position,1:3,:);
        RASI_xyz_position = reshape(RASI_xyz_position, [3, time])/10;
        LASI_position = find(strcmp( path, 'LASI'));
        LASI_xyz_position = name.Trajectories.Labeled.Data(LASI_position,1:3,:);
        LASI_xyz_position = reshape(LASI_xyz_position, [3, time])/10;


        %%%% Distance calculation %%%%
        R_xyz_distance = abs(RWRL_xyz_position - RASI_xyz_position);
        R_x_sqr = R_xyz_distance(1,:) .* R_xyz_distance(1,:);
        R_y_sqr = R_xyz_distance(1,:) .* R_xyz_distance(2,:);
        R_z_sqr = R_xyz_distance(1,:) .* R_xyz_distance(3,:);
        R_distance = sqrt(R_x_sqr + R_y_sqr + R_z_sqr);

        L_xyz_distance = abs(LWRL_xyz_position - LASI_xyz_position);
        L_x_sqr = L_xyz_distance(1,:) .* L_xyz_distance(1,:);
        L_y_sqr = L_xyz_distance(1,:) .* L_xyz_distance(2,:);
        L_z_sqr = L_xyz_distance(1,:) .* L_xyz_distance(3,:);
        L_distance = sqrt(R_x_sqr + R_y_sqr + R_z_sqr);

        %}

        R_L_xyz_disatnce = abs(LWRL_xyz_position - RWRL_xyz_position);
        R_L_x_sqr = R_L_xyz_disatnce(1,:) .* R_L_xyz_disatnce(1,:);
        R_L_y_sqr = R_L_xyz_disatnce(2,:) .* R_L_xyz_disatnce(2,:);
        R_L_z_sqr = R_L_xyz_disatnce(3,:) .* R_L_xyz_disatnce(3,:);
        R_L_distance = sqrt(R_L_x_sqr + R_L_y_sqr + R_L_z_sqr);

        %{
        figure

        plot(R_L_distance)
        hold on
        title('Lifting hands off iliac crests')
        %}

        features_RL_dis = [nanstd(R_L_distance), nanmax(R_L_distance) - nanmin(R_L_distance)];
        RL_dis_col = ["var_RL_dis", "Max_min_RL_dis"];
        RL_dis_df = [RL_dis_col;features_RL_dis];

        hand_iliac_count = 0;
        for i = 1:(time-1)
            try
                if abs(R_L_distance(i)) < (nanmean(R_L_distance)+nanstd(R_L_distance)*3.5) &&  abs(R_L_distance(i+1)) > (nanmean(R_L_distance)+nanstd(R_L_distance)*3.5) %%% mean + st*(optimal value)
                    hand_iliac_count = hand_iliac_count + 1;
                end
            catch
                continue
            end
        end

        %% 5. Remaining out of the test position for more than 5s

        out_position_count = 0;
        for i = 1:(time-1)
            try
                if (abs(landing_foot_height(i)) > 10 && abs(landing_foot_height(i+5000)) > 10) || (abs(main_moving_leg(i)) > 30 && abs(main_moving_leg(i+5000)) > 30) || (abs(swinging_FP(i)) > 100 && abs(swinging_FP(i+5000)) > 100) %%% stance (instead of supporting) and swing 
                    out_position_count = out_position_count + 1;
                end
            catch
                continue
            end
        end


        %% Summary of BESS
        total_score = foot_height_count + angle_count + stepping_count + hand_iliac_count + out_position_count;
        BESS_value = [sub_ID, foot_height_count, angle_count, stepping_count, hand_iliac_count, out_position_count, total_score];
        BESS_col = ["sub_ID","foot_height_count", "angle_count", "stepping_count", "hand_iliac_count", "out_position_count", "totel"];
        BESS_summary = [BESS_col; BESS_value];


        %% save file
        final_df = [hei_df, flex_ang_df, abd_ang_df, stepping_df, RL_dis_df, BESS_summary];
        saved_path_name = "C:\Users\a1003\OneDrive\桌面\Thesis\data\OLBT\sub1to31\" + sub_ID + ".csv";
        writematrix(final_df, saved_path_name)
    catch
        subjfile
        %clear
        addpath 'D:\files\Thesis_first_part_data'
        addpath 'C:\Users\a1003\OneDrive\桌面\Thesis'
        filenames = dir('D:\files\Thesis_first_part_data');
        continue
    end
end



