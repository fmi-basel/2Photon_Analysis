%% Generate Experimental settings file

function Generate_ExpDetails(selpath)

    % Get Folder
    dims = [1 50];
    formatOut = 'dd/mm/yy';

    % Get the main Information
    title = 'Experimental Details';
    prompt = {'Date', 'Animal number', 'Number of planes', 'Spacing [um]', 'Number of channels', 'Imaging Rate (Hz)', ...
        'td-Tomato', 'Indicator 1', 'Promoter 1', 'AAV-Serotype 1', 'Indicator 2',  'Promoter 2', 'AAV-Serotype 2', ...
        'Lens location', 'Surgical Implantation Side (DV, AP, MM)', 'Age of animal (weeks)', 'Sex', ...
        'Habituation to head fixation (Ses.)', 'Gain', 'Laser Power', 'Protocol' 'Notes'};

    definput = {datestr(now,formatOut), '888888', '3', '80', '1', '30', '0', 'GCaMP6f', ...
        'CamKII', '9', 'flex-tdTom', 'cag', '1', 'Left', '4.2, 1.5, 3.3' '12', 'male', '3', '600', '300', 'US',''};
    answer = inputdlg(prompt,title,dims,definput);

    % Prepare Output and save in Excel sheet
    savefast([selpath '\Experimental_Details.mat'], 'answer')

end





