clear all;
FD_lst = importdata('./fdd.lst');
% FD_lst=FD_lst_old(~ismember(FD_lst_old,"CQ.HCB"));
for index = 1:1:length(FD_lst)
    %for index=1:1:1
    clear DATAALL;
    clear S;
    clear DATAALL_SORT_RMSE;
    clear sacfilename;
    StaCode = char(FD_lst(index));
    %StaCode='CQ.ROC'
    sacfilename = dir(['./' StaCode '/*.srf']);
    nums = size(sacfilename);
    plotfig = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1:1:nums(1)
        sacnm = [sacfilename(ii).folder, '/', sacfilename(ii).name];
        S = readsac(sacnm);
        TimeSerie = (S.B):S.DELTA:(S.E + S.DELTA);
        DATAALL(ii, :) = S.DATA1';

        if (plotfig == 1)
            plot(TimeSerie, S.DATA1, 'k', 'Linewidth', 2);
            xlim([S.B S.E]);
            hold on;
        end

    end

    sumtrace = sum(DATAALL) / nums(1);

    if (plotfig == 1)
        plot(TimeSerie, sumtrace, 'r', 'Linewidth', 2);
    end

    for ii = 1:1:nums(1)
        rmse = sqrt(mean((DATAALL(ii, 1:1001) - sumtrace) .^ 2));
        DATAALL(ii, 1002) = rmse;
        DATAALL(ii, 1003) = ii;
    end

    if (plotfig == 1)
        figure
        subplot(121);
        imagesc(DATAALL(:, 1:1001));
        subplot(122);
        plot(DATAALL(:, 1002), DATAALL(:, 1003));
        ylim([0 nums(1)]);
        set(gca, 'Ydir', 'reverse')
    end

    if (plotfig == 1)
        figure;
    end

    DATAALL_SORT_RMSE = sortrows(DATAALL, 1002);

    if (plotfig == 1)
        subplot(121);
        imagesc(DATAALL_SORT_RMSE(:, 1:1001));
        subplot(122);
    end

    for ii = 1:1:nums(1)
        DATAALL_SORT_RMSE(ii, 1004) = ii;
    end

    if (plotfig == 1)
        plot(DATAALL_SORT_RMSE(:, 1002), DATAALL_SORT_RMSE(:, 1004));
        ylim([0 nums(1)]);
        set(gca, 'Ydir', 'reverse')
    end

    if (plotfig == 1)
        figure;

        for ii = 1:1:floor(nums(1) * 0.25)
            plot(TimeSerie, DATAALL_SORT_RMSE(ii, 1:1001), 'k', 'Linewidth', 2);
            xlim([S.B S.E]);
            hold on;
        end

    end

    sumtrace_rmse_0_05 = sum(DATAALL_SORT_RMSE(1:floor(nums(1) * 0.05), 1:1001)) / floor(nums(1) * 0.05);
    sumtrace_rmse_0_10 = sum(DATAALL_SORT_RMSE(1:floor(nums(1) * 0.10), 1:1001)) / floor(nums(1) * 0.10);
    sumtrace_rmse_0_15 = sum(DATAALL_SORT_RMSE(1:floor(nums(1) * 0.15), 1:1001)) / floor(nums(1) * 0.15);
    sumtrace_rmse_0_20 = sum(DATAALL_SORT_RMSE(1:floor(nums(1) * 0.20), 1:1001)) / floor(nums(1) * 0.20);
    sumtrace_rmse_0_25 = sum(DATAALL_SORT_RMSE(1:floor(nums(1) * 0.25), 1:1001)) / floor(nums(1) * 0.25);
    sumtrace_rmse_0_30 = sum(DATAALL_SORT_RMSE(1:floor(nums(1) * 0.30), 1:1001)) / floor(nums(1) * 0.30);
    sumtrace_rmse_0_40 = sum(DATAALL_SORT_RMSE(1:floor(nums(1) * 0.40), 1:1001)) / floor(nums(1) * 0.40);
    sumtrace_rmse_0_50 = sum(DATAALL_SORT_RMSE(1:floor(nums(1) * 0.50), 1:1001)) / floor(nums(1) * 0.50);
    sumtrace_rmse_0_60 = sum(DATAALL_SORT_RMSE(1:floor(nums(1) * 0.60), 1:1001)) / floor(nums(1) * 0.60);
    sumtrace_rmse_0_70 = sum(DATAALL_SORT_RMSE(1:floor(nums(1) * 0.70), 1:1001)) / floor(nums(1) * 0.70);
    sumtrace_rmse_0_80 = sum(DATAALL_SORT_RMSE(1:floor(nums(1) * 0.80), 1:1001)) / floor(nums(1) * 0.80);
    sumtrace_rmse_0_90 = sum(DATAALL_SORT_RMSE(1:floor(nums(1) * 0.90), 1:1001)) / floor(nums(1) * 0.90);

    if (plotfig == 1)
        plot(TimeSerie, sumtrace_rmse_0_20, 'g', 'Linewidth', 2);
        hold on;
        plot(TimeSerie, sumtrace, 'r', 'Linewidth', 2);
    end

    if (plotfig == 1)
        figure;
        plot(TimeSerie, sumtrace_rmse_0_05, 'r+');
        hold on;
        plot(TimeSerie, sumtrace_rmse_0_10, 'r-');
        hold on;
        plot(TimeSerie, sumtrace_rmse_0_15, 'ro');
        hold on;
        plot(TimeSerie, sumtrace_rmse_0_20, 'r*');
        hold on;
        plot(TimeSerie, sumtrace_rmse_0_25, 'rx');
        hold on;
        plot(TimeSerie, sumtrace_rmse_0_30, 'k+');
        hold on;
        plot(TimeSerie, sumtrace_rmse_0_40, 'k-');
        hold on;
        plot(TimeSerie, sumtrace_rmse_0_50, 'ko');
        hold on;
        plot(TimeSerie, sumtrace_rmse_0_60, 'k*');
        hold on;
        plot(TimeSerie, sumtrace_rmse_0_70, 'kx');
        hold on;
        plot(TimeSerie, sumtrace_rmse_0_80, 'b+');
        hold on;
        plot(TimeSerie, sumtrace_rmse_0_90, 'b-');

        hold on;
        plot(TimeSerie, sumtrace, 'g', 'Linewidth', 2);
        xlim([S.B S.E]);
        legend('0.05', '0.10', '0.15', '0.20', '0.25', '0.30', '0.40', '0.50', '0.60', '0.70', '0.80', '0.90', '1.00');
    end

    Name_Out = ['./' StaCode '/' StaCode '.rmse10.lst'];
    Name_DatOut = ['./' StaCode '/' StaCode '.rmse10.rmse.dat'];
    fid = fopen(Name_Out, 'w');
    fid_dat = fopen(Name_DatOut, 'w');

    for ii = 1:1:floor(nums(1) * 0.10)
        name_index = DATAALL_SORT_RMSE(ii, 1003);
        fprintf(fid, 'cp %s %s.rmse10\n', sacfilename(name_index).name, sacfilename(name_index).name);
        fprintf(fid_dat, 'sac<<eof\nr %s.rmse10\n ch user1 %f\nwh\nq\neof\n', sacfilename(name_index).name, DATAALL_SORT_RMSE(ii, 1002));
    end

    Name_Out1 = ['./' StaCode '/' StaCode '.rmse20.lst'];
    Name_DatOut1 = ['./' StaCode '/' StaCode '.rmse20.rmse.dat'];
    fid_dat1 = fopen(Name_DatOut1, 'w');
    fid1 = fopen(Name_Out1, 'w');

    for ii = 1:1:floor(nums(1) * 0.20)
        name_index = DATAALL_SORT_RMSE(ii, 1003);
        fprintf(fid1, 'cp %s %s.rmse20\n', sacfilename(name_index).name, sacfilename(name_index).name);
        fprintf(fid_dat1, 'sac<<eof\nr %s.rmse20\n ch user1 %f\nwh\nq\neof\n', sacfilename(name_index).name, DATAALL_SORT_RMSE(ii, 1002));
    end

    Name_Out2 = ['./' StaCode '/' StaCode '.rmse30.lst'];
    Name_DatOut2 = ['./' StaCode '/' StaCode '.rmse30.rmse.dat'];
    fid_dat2 = fopen(Name_DatOut2, 'w');
    fid2 = fopen(Name_Out2, 'w');

    for ii = 1:1:floor(nums(1) * 0.30)
        name_index = DATAALL_SORT_RMSE(ii, 1003);
        fprintf(fid2, 'cp %s %s.rmse30\n', sacfilename(name_index).name, sacfilename(name_index).name);
        fprintf(fid_dat2, 'sac<<eof\nr %s.rmse30\n ch user1 %f\nwh\nq\neof\n', sacfilename(name_index).name, DATAALL_SORT_RMSE(ii, 1002));
    end

    Name_DatOutall = ['./' StaCode '/' StaCode '.all.rmse.dat'];
    fid_datall = fopen(Name_DatOutall, 'w');

    for ii = 1:1:nums(1)
        name_index = DATAALL_SORT_RMSE(ii, 1003);
        fprintf(fid_datall, 'sac<<eof\nr %s\n ch user1 %f\nwh\nq\neof\n', sacfilename(name_index).name, DATAALL_SORT_RMSE(ii, 1002));
    end

    fclose(fid);
    fclose(fid1);
    fclose(fid2);
    fclose(fid_dat);
    fclose(fid_dat1);
    fclose(fid_dat2);
    fclose(fid_datall);

    grdname_before = ['./' StaCode '/' StaCode 'SRF_before.xyz'];
    rmse_before = ['./' StaCode '/' StaCode 'rmse_before.xy'];
    grdname_after = ['./' StaCode '/' StaCode 'SRF_after.xyz'];
    rmse_after = ['./' StaCode '/' StaCode 'rmse_after.xy'];

    fid_grdBefore = fopen(grdname_before, 'w');
    fid_grdAfter = fopen(grdname_after, 'w');
    fid_rmseB = fopen(rmse_before, 'w');
    fid_rmseA = fopen(rmse_after, 'w');

    for ii = 1:1:nums(1)
        fprintf(fid_rmseB, '%f %f\n', ii, DATAALL(ii, 1002));
        fprintf(fid_rmseA, '%f %f\n', ii, DATAALL_SORT_RMSE(ii, 1002));

        for jj = 1:1:1001
            time_index = TimeSerie(jj);
            name_index = DATAALL_SORT_RMSE(ii, 1003);
            fprintf(fid_grdBefore, '%f %f %f\n', time_index, ii, DATAALL(ii, jj));
            fprintf(fid_grdAfter, '%f %f %f\n', time_index, ii, DATAALL_SORT_RMSE(ii, jj));
        end

    end

    fclose(fid_grdBefore);
    fclose(fid_grdAfter);
    fclose(fid_rmseB);
    fclose(fid_rmseA);

    Ref_Trace = ['./' StaCode '/' StaCode '.ref_sum.xy'];
    fid_reft = fopen(Ref_Trace, 'w');

    for jj = 1:1:1001
        time_index = TimeSerie(jj);
        fprintf(fid_reft, '%.1f %f %f %f %f %f %f %f %f %f %f %f %f %f\n', ...
            time_index, sumtrace(jj), sumtrace_rmse_0_05(jj), sumtrace_rmse_0_10(jj), ...
            sumtrace_rmse_0_15(jj), sumtrace_rmse_0_20(jj), sumtrace_rmse_0_25(jj), ...
            sumtrace_rmse_0_30(jj), sumtrace_rmse_0_40(jj), sumtrace_rmse_0_50(jj), ...
            sumtrace_rmse_0_60(jj), sumtrace_rmse_0_70(jj), sumtrace_rmse_0_80(jj), ...
            sumtrace_rmse_0_90(jj));
    end

    fclose(fid_reft);
    disp([StaCode ' DONE ~']);
end
