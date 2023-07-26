function newSessionID = convertSessionID(sID,rootDir)
% function newSessionID = convertSessionID(sID,rootDir)
% convert sessionID from our lab's raw one to abbreviated one and back

newSessionID = '';

% identify whether current format is raw or abbreviated
rawSflag = false;
if contains(sID,'.')
    rawSflag = true;
end

if rawSflag  % convert from raw to abbreviated
    us = strfind(sID,'_');
    bDir = sID(1:us(2)-1);
    cd(rootDir);
    dn=dir([bDir '*']);
    if length(dn)>1
        dots = strfind(sID,'.');
        bDir = sID(1:dots(1)-1);
        dn = dir([bDir '*']);
    end
    if isempty(dn)
        fprintf('%s not available\n',bDir);
        return
    else
        newSessionID = dn(1).name;
    end
    
else  % convert from abbreviated to raw
    sID = strrep(sID,'all1','all.');
    cd(rootDir)
    dn = dir([sID '*.CO']);
    
    if isempty(dn)
        dn = dir([sID '*rds*DX2']);
        if isempty(dn)
            fprintf('%s not available\n',sID);
        return
        else
           newSessionID = dn(end).name;
        end
    else
        newSessionID = dn(end).name;
    end
end
