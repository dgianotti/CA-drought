clear;
clc;

% Determine which stations are in CA:
all_IDs = load_goodIDs;

CA_IDs = {};
bad_states = {};

n = 0;
m = 0;

for i = 1:774
    id = all_IDs{i};
    stn = load_stn_data(id,'stn');
    state = stn.human_location{3};
    
    if strcmp(state,'CA')
        % it's in CA!
        n = n+1;
        CA_IDs{n} = id;
    else
        % Not in CA, but just to check
        m = m+1;
        bad_states{m} = state;
    end
end
      
save('CA_ids.mat','CA_IDs');