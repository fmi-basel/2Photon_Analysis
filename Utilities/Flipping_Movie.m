%% Fixing Andras's Setting change
Session = { ...
    '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Imaging_2p\US\2021_08\A_1040736\Session_18\processed_data_Folder\', ...
    '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Imaging_2p\US\2021_08\A_1040736\Session_19\processed_data_Folder\', ...
    '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Imaging_2p\US\2021_08\A_1040738\Session_17\processed_data_Folder\', ...
    '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Imaging_2p\US\2021_08\A_1040738\Session_18\processed_data_Folder\', ...
    '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Imaging_2p\US\2021_08\A_1040738\Session_19\processed_data_Folder\', ...
    '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Imaging_2p\US\2021_08\A_1040743\Session_17\processed_data_Folder\', ...
    '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Imaging_2p\US\2021_08\A_1040743\Session_18\processed_data_Folder\', ...
    '\\tungsten-nas.fmi.ch\tungsten\scratch\gluthi\hinzjuli\Data\Imaging_2p\US\2021_08\A_1040743\Session_19\processed_data_Folder\'};

for c = 1:size(Session, 2)
    
    % CorrImages
    load([Session{c} 'Corr_Imgs.mat'])
    obj.Cor_Images_1.Corr_Img = flipdim(obj.Cor_Images_1.Corr_Img, 2);
    obj.Cor_Images_2.Corr_Img = flipdim(obj.Cor_Images_2.Corr_Img, 2);
    obj.Cor_Images_3.Corr_Img = flipdim(obj.Cor_Images_3.Corr_Img, 2);
    savefast([Session{c} 'Corr_Imgs.mat'], 'obj')
    clear obj

    % DS_MAT
    load([Session{c} 'DS_Dat.mat'])
    DS_Dat = flipdim(DS_Dat, 2);
    savefast([Session{c} 'DS_Dat.mat'], 'DS_Dat')
    clear DS_Dat

    % DS_MAT_vis
    load([Session{c} 'DS_Dat_vis.mat'])
    DS_Dat_vis = flipdim(DS_Dat_vis, 2);
    savefast([Session{c} 'DS_Dat_vis.mat'], 'DS_Dat_vis')
    clear DS_Dat_vis

    % MC_r
    load([Session{c} 'MC_r.mat'])
    Y = flipdim(Y, 2);
    savefast([Session{c} 'MC_r.mat'], 'Y')
    clear Y

    % Motion_correction_templates
    load([Session{c} 'Motion_correction_templates.mat'])
    Templates_rg{1} = flipdim(Templates_rg{1}, 2);
    Templates_rg{2} = flipdim(Templates_rg{2}, 2);
    Templates_rg{3} = flipdim(Templates_rg{3}, 2);
    savefast([Session{c} 'Motion_correction_templates.mat'], 'Templates_rg')
    clear Templates_rg

    % Raw
    load([Session{c} 'Raw.mat'])
    Y = flipdim(Y, 2);
    savefast([Session{c} 'Raw.mat'], 'Y')
    clear Y
    
    % Delete related products
    try
        delete([Session{c} 'Master_structure.mat']); 
        delete([Session{c} 'PospP_CNMF_SVM.mat']); 
        delete([Session{c} 'Results_CNMF.mat']); 
    catch
       disp('No Post-Processing, yet!') 
    end
end