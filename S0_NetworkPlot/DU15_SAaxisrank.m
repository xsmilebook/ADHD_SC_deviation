clc;
clear;

ResultsFolder = 'D:\xuxiaoyu\DMRI_network_development\Normative_model\results\';
SAaxis = cifti_read('D:\xuxiaoyu\Atlas\S_A_axis\S-A_ArchetypalAxis\FSLRVertex\SensorimotorAssociation_Axis.dscalar.nii');
Du15_Data_Mat =cifti_read('D:\xuxiaoyu\Atlas\DU15NET\DU15NET\HCP\fsLR_32k\DU15NET_Consensus_fsLR_32k.dlabel.nii');
Du15_Group_Label_L = cifti_struct_dense_extract_surface_data(Du15_Data_Mat,'CORTEX_LEFT');
Du15_Group_Label_R = cifti_struct_dense_extract_surface_data(Du15_Data_Mat,'CORTEX_RIGHT');
Du15_Group_Label = [Du15_Group_Label_L;Du15_Group_Label_R];
tabulate(Du15_Group_Label)

SAaxis_Label_L = cifti_struct_dense_extract_surface_data(SAaxis,'CORTEX_LEFT');
SAaxis_Label_R = cifti_struct_dense_extract_surface_data(SAaxis,'CORTEX_RIGHT');
SAaxis_Label = [SAaxis_Label_L;SAaxis_Label_R];

MeanSA = groupsummary(SAaxis_Label, Du15_Group_Label, 'mean');
MedianSA = groupsummary(SAaxis_Label, Du15_Group_Label, 'median');

NetworkLabel = ["NONE", "VIS-P", "CG-OP", "DN-B", "SMOT-B", "AUD", "PM-PPr", "dATN-B", "SMOT-A", "LANG", "FPN-B", "FPN-A", "dATN-A", "VIS-C", "SAL/PMN", "DN-A", "UNCERTAIN"]';
Du15_SAaxis = table(MeanSA, MedianSA, NetworkLabel);
writetable(Du15_SAaxis, [ResultsFolder, '\Du15_SArank.csv']);

Du15_SAaxis_vertex = table(Du15_Group_Label, SAaxis_Label);
writetable(Du15_SAaxis_vertex, [ResultsFolder, '\Du15_SArank_vertex.csv']);


