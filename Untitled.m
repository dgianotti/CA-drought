
clear;
clc;

load('CA_ids.mat')

dat = nan(125,13);

for i = 1:13
    id = good_CA_IDs{i};
    filename = ['GHCN-Daily/',id,'.dly'];
    
    
    % Read in fixed-width data:
    fprintf('Re-formatting data for station %s...\n',id);
    [datenums,precip] = format_GHCN_precip_data(filename,'PadFirstLastYears','ExcludeLeapDays');

    %     % Set last day to to:
    %     feb28_2014 = datenum('2014-02-28','yyyy-mm-dd');
    %     datenums = [datenums(:); ((datenums(end)+1):feb28_2014)'];
    %     precip = [precip(:); nan([length(datenums)-length(precip),1])];
        
    new_precip = precip(datenums>= datenum('2010-01-01') & datenums<= today);
    new_datenums = datenums(datenums>= datenum('2010-01-01') & datenums<= today);

    
    % Load old precip data:
    ImpStn = load_stn_data(id,'ImpStn');

    % Determine likelihood of data given model:
    
    nan_padded_new_data = nan(365*5,1);
    nan_padded_new_data(1:length(new_precip)) = new_precip;
    
    data = [ImpStn.intensity_data; reshape(nan_padded_new_data,[365,5])'];
    
    data(end,:) = [];
    
    dat( (end-size(data,1) +1):end,i) = nanmean(data,2);
    
end

%%

    figure;
    for i = 1:13
        plot(1888:2013,365*dat(:,i));
        hold on;
        
    end
    
    
    print(gcf,'-dpng','annual_precip.png');
    
    