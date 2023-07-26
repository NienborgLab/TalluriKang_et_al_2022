
function alignTriggers(sesID,trigDir,tsDir,outDir)

% load triggers
load(fullfile(trigDir,sesID,'both_trigTimes'),'StM','StNMM','evM','evNMM');

% load time stamps
load(fullfile(tsDir,sesID,'preprocessed_data'),'D');
preprocessed_time_stamps = D.time_stamps;

fprintf('processing... %s \n',sesID)
%% align triggers to time-stamps
trigIdxS    = cell(size(StM));
trigIdxSNM  = cell(size(StM));
trigIdxM    = cell(size(StM));
trigIdxNM   = cell(size(StM));
trigIdxM_u    = cell(size(StM));
trigIdxNM_u   = cell(size(StM));

for h = 1:size(StM,1)
    for s = 1:size(StM,2)
        for n =1:length(StM{h,s})
            [val,minIdx] = min(abs(preprocessed_time_stamps-StM{h,s}(n)));
            if val<0.02  % only include time-stamps that are within one video frame
                trigIdxS{h,s} = [trigIdxS{h,s},minIdx];
            else
                StM{h,s}(n) = NaN;
            end
        end
        for n = 1:length(StNMM{h,s})
            [val,minIdx] = min(abs(preprocessed_time_stamps-StNMM{h,s}(n)));
            if val<0.02 % align to within 20ms (frame duration is 16ms)
                trigIdxSNM{h,s} = [trigIdxSNM{h,s} minIdx];
            else 
                StNMM{h,s}(n) = NaN;
            end
        end
        for n =1:length(evM{h,s})
            [val,minIdx] = min(abs(preprocessed_time_stamps-evM{h,s}(n)));
            if val<0.02
                trigIdxM{h,s} = [trigIdxM{h,s} minIdx];
            else
                evM{h,s}(n) = NaN;
            end
        end
        for n = 1:length(evNMM{h,s})
            [val,minIdx] = min(abs(preprocessed_time_stamps-evNMM{h,s}(n)));
            if val <0.02
                trigIdxNM{h,s} = [trigIdxNM{h,s} minIdx];
            else
                evNMM{h,s}(n) = NaN;
            end
        end

    end
end


if contains(sesID,'le_')
    dirNm = strrep(sesID,'le_','M1_');
else
    dirNm = strrep(sesID,'ma_','M2_');
end

if ~exist(fullfile(outDir,dirNm),'dir')
    mkdir(fullfile(outDir,dirNm))
end
cd(fullfile(outDir,dirNm))
save('both_trigIdx','trigIdxS','trigIdxSNM','trigIdxM','trigIdxNM')

cd(fullfile(trigDir,dirNm))
save('both_trigTimes_clean','StM','StNMM','evM','evNMM')



