clc;
clear;

ResultsFolder = 'D:\xuxiaoyu\DMRI_network_development\Normative_model\results_ABCD\';
SAaxis = cifti_read('D:\xuxiaoyu\Atlas\S_A_axis\S-A_ArchetypalAxis\FSLRVertex\SensorimotorAssociation_Axis.dscalar.nii');

%% Yeo17 vertex
Yeo17_Data_Mat =cifti_read('D:\xuxiaoyu\Atlas\Yeo2011_17Networks_N1000.dlabel.nii');
Yeo17_Group_Label_L = cifti_struct_dense_extract_surface_data(Yeo17_Data_Mat,'CORTEX_LEFT');
Yeo17_Group_Label_R = cifti_struct_dense_extract_surface_data(Yeo17_Data_Mat,'CORTEX_RIGHT');
Yeo17_Group_Label = [Yeo17_Group_Label_L;Yeo17_Group_Label_R];
tabulate(Yeo17_Group_Label)

SAaxis_Label_L = cifti_struct_dense_extract_surface_data(SAaxis,'CORTEX_LEFT');
SAaxis_Label_R = cifti_struct_dense_extract_surface_data(SAaxis,'CORTEX_RIGHT');
SAaxis_Label = [SAaxis_Label_L;SAaxis_Label_R];

MeanSA = groupsummary(SAaxis_Label, Yeo17_Group_Label, 'mean');
MedianSA = groupsummary(SAaxis_Label, Yeo17_Group_Label, 'median');

NetworkLabel = ["None", "VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", ...
    "LimbicA", "LimbicB",  "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB"]';
Yeo17_SAaxis = table(MeanSA, MedianSA, NetworkLabel);
writetable(Yeo17_SAaxis, [ResultsFolder, '\Yeo17_SArank.csv']);

Yeo17_SAaxis_vertex = table(Yeo17_Group_Label, SAaxis_Label);
writetable(Yeo17_SAaxis_vertex, [ResultsFolder, '\Yeo17_SArank_vertex.csv']);

%% Schaefer 400 Yeo 17 vertex
Yeo17S400_Data_Mat =cifti_read('D:\xuxiaoyu\Atlas\qsirecon_AtlasPack_loc\tpl-fsLR_atlas-4S456Parcels_den-91k_dseg.dlabel.nii');
Yeo17S400_Group_Label_L = cifti_struct_dense_extract_surface_data(Yeo17S400_Data_Mat,'CORTEX_LEFT');
Yeo17S400_Group_Label_R = cifti_struct_dense_extract_surface_data(Yeo17S400_Data_Mat,'CORTEX_RIGHT');
Yeo17S400_Group_Label = [Yeo17S400_Group_Label_L;Yeo17S400_Group_Label_R];
Yeo17S400_Group_Label_NetLabel = zeros(length(Yeo17S400_Group_Label), 1);
S400tab = readtable('D:\xuxiaoyu\Atlas\qsirecon_AtlasPack_loc\schaefer400_index_SA.csv');
NetworkLabel = ["None", "VisCent", "VisPeri", "SomMotA", "SomMotB", "DorsAttnA", "DorsAttnB", "SalVentAttnA", "SalVentAttnB", ...
    "LimbicA", "LimbicB",  "ContC", "ContA", "ContB", "TempPar", "DefaultC", "DefaultA", "DefaultB"]';

for i = 2:18
    Nettmp = NetworkLabel(i);
    regionidx = S400tab.index(find(strcmp(S400tab.network_label_17network, Nettmp)));
    idx = ismember(Yeo17S400_Group_Label, regionidx);
    Yeo17S400_Group_Label_NetLabel(idx) = i-1;
end
tabulate(Yeo17S400_Group_Label_NetLabel)

MeanSA = groupsummary(SAaxis_Label, Yeo17S400_Group_Label_NetLabel, 'mean');
MedianSA = groupsummary(SAaxis_Label, Yeo17S400_Group_Label_NetLabel, 'median');

Yeo17S400_SAaxis = table(MeanSA, MedianSA, NetworkLabel);
writetable(Yeo17S400_SAaxis, [ResultsFolder, '\Yeo17S400_SArank.csv']);

Yeo17S400_SAaxis_vertex = table(Yeo17S400_Group_Label_NetLabel, SAaxis_Label);
writetable(Yeo17S400_SAaxis_vertex, [ResultsFolder, '\Yeo17S400_SArank_vertex.csv']);

Yeo17_Data_Mat.cdata = Yeo17S400_Group_Label_NetLabel;

cifti_write(Yeo17_Data_Mat, 'D:\xuxiaoyu\Atlas\qsirecon_AtlasPack_loc\tpl-fsLR_atlas-4SYeo17Parcels_den-91k_dseg.dlabel.nii');
