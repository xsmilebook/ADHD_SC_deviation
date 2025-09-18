% convert Yeo17 52 regions to 17 network
clc; clear;

atlaspath='D:\xuxiaoyu\Atlas\CBIG-stable_projects-brain_parcellation-Yeo2011_fcMRI_clustering\1000subjects_reference\Yeo_JNeurophysiol11_SplitLabels';
Yeo17MNI1mm=niftiread([atlaspath, '\MNI152\Yeo2011_17Networks_N1000.split_components.FSL_MNI152_1mm.nii.gz']);
MNI1mminfo = niftiinfo([atlaspath, '\MNI152\Yeo2011_17Networks_N1000.split_components.FSL_MNI152_1mm.nii.gz']);
indextable=readtable([atlaspath, '\Yeo2011_17networks_N1000.split_components.glossary.csv']);
indextable.index=(0:114)';

% 1. do not dinguish left and right
networklabel=unique(indextable.NetworkName);
networklabel=networklabel([1,2,13,16,17,9,10,14,15,11,12,5,3,4,18,8,6,7]);

indextable.networkindex = categorical(indextable.NetworkName, networklabel');
indextable.networkindex = double(indextable.networkindex);
indextable.networkindex = indextable.networkindex - 1;

Yeo17MNI1mm_Network = zeros(182, 218, 182);
for i = 0:17
    idxvect = indextable.index(find(indextable.networkindex==i));
    idx = ismember(Yeo17MNI1mm, idxvect);
    Yeo17MNI1mm_Network(idx) = i;
end
Yeo17MNI1mm_Network = single(Yeo17MNI1mm_Network);

niftiwrite(Yeo17MNI1mm_Network, [atlaspath, '\MNI152\Yeo2011_17Networks_FSL_MNI152_1mmNetworkIndex'], MNI1mminfo);

% 2. dinguish left and right
indextable.NetworkNameRL = indextable.NetworkName;
indextable.NetworkNameRL(2:58) = cellfun(@(x) [x '_LH'], indextable.NetworkNameRL(2:58), 'UniformOutput', false);
indextable.NetworkNameRL(59:115) = cellfun(@(x) [x '_RH'], indextable.NetworkNameRL(59:115), 'UniformOutput', false);

networklabel_L=unique(indextable.NetworkName);
networklabel_L=networklabel_L([1,2,13,16,17,9,10,14,15,11,12,5,3,4,18,8,6,7]);
networklabel_L(2:18)=cellfun(@(x) [x '_LH'], networklabel_L(2:18), 'UniformOutput', false);

networklabel_R=unique(indextable.NetworkName);
networklabel_R=networklabel_R([1,2,13,16,17,9,10,14,15,11,12,5,3,4,18,8,6,7]);
networklabel_R=cellfun(@(x) [x '_RH'], networklabel_R(2:18), 'UniformOutput', false);

networklabel_LR = [networklabel_L; networklabel_R];
indextable.networkindexLR = categorical(indextable.NetworkNameRL, networklabel_LR');
indextable.networkindexLR = double(indextable.networkindexLR);
indextable.networkindexLR = indextable.networkindexLR - 1;

Yeo17MNI1mm_Network_LR = zeros(182, 218, 182);
for i = 0:34
    idxvect = indextable.index(find(indextable.networkindexLR==i));
    idx = ismember(Yeo17MNI1mm, idxvect);
    Yeo17MNI1mm_Network_LR(idx) = i;
end
Yeo17MNI1mm_Network_LR = single(Yeo17MNI1mm_Network_LR);

niftiwrite(Yeo17MNI1mm_Network_LR, [atlaspath, '\MNI152\Yeo2011_17Networks_FSL_MNI152_1mmNetworkIndex_LR'], MNI1mminfo);



