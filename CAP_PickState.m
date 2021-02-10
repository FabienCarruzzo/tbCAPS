function [State2] = CAP_PickState(State,TM)

    % Picks a random number between 0 and 1
    T = randi([0,100000])/100000;

    % Our probabilities of interest (knowing the start state)
    POI = TM(State,:);
    CP = cumsum(POI);
    
    if T <= CP(1)
        State2 = 1;
    elseif T <= CP(2)
        State2 = 2;
    elseif T <= CP(3)
        State2 = 3;
    elseif T <= CP(4)
        State2 = 4;
    elseif T <= CP(5)
        State2 = 5;
    else
        errordlg('WTH');
    end
end