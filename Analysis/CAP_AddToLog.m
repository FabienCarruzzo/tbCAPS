% Displays a message specifying what command is currently being run, and
% saves the same command and relevant information into a text file
function [Log] = CAP_AddToLog(Log,new_s,Params,Params_name)

    n = length(Log);
    Log{n+1}{1} = [' '];
    %Log{n+1}{1} = [date,'/',num2str(round(hour(now))),'/',num2str(round(minute(now))),'/',num2str(round(second(now)))];
    Log{n+1}{2} = new_s;
    
    disp(new_s);
    
    % If we also want to save parameters, we do so
    if nargin > 2
        for n2 = 1:length(Params)
            Log{n+1}{n2+2} = [Params_name{n2},': ',num2str(Params{n2})];
            disp(Log{n+1}{n2+2});
        end
    end
end