%% Finds the type of data that has been entered, using the file type of the
% matlab variable
function DataType = CAP_FindDataType(Data)

    if iscell(Data)
        DataType = 'Data';
    elseif isfloat(Data)
        DataType = 'Motion';
    elseif islogical(Data)
        DataType = 'Mask';
    elseif isstruct(Data)
        DataType = 'Info';
    else
        DataType = 'Unknown';
    end

end