%%%%%%% tomo_cplxgenerate - This file will simulate protein complexes from template library, and keep track of the minimum bounding sphere radius and copy numbers of each complex bbchange

addpath(genpath(pwd)); % I don't think this does anything.

tic; clear %bbchange

tomogram_size=[100 100 100]; % bbchange - was 500 500 200, and was in tomo_densgenerate part I THINK THIS IS NM
pym=importdata('Alber_model_1d.txt'); % bb INPUT FILE - Contains locations and rotations of all particles. CHANGE THIS as needed, then change vol_dim (line 87) to match, and change output filename in tom_mrcwrite (line 233)

pym=pym.data;

rots=[];
for i=1:size(pym,1)/2
        rots(i,1:3)=pym(2*i-1,1:3);
end
disp(rots)

disp('Initialization Start');

work_dir=pwd;
ws = struct();
ws.map.map_resolution_s = '40'; %bb This (cubed) is the size of all of the individual protein maps. They will all be the same size. Situs generates maps of different sizes, so they need to be rescaled
ws.map.map_resolution = str2double(ws.map.map_resolution_s);
ws.map.rescale_individual_maps = false;
ws.map.protein_names = {'1BXR'}%{'1A1S','1EQR','1GYT','1KYI','1VPX','1W6T','2AWB','2BYU','2GLS','2IDB','3DY4','1BXR','1F1B','1KP8','1QO1','1VRG','1YG6','2BO9','2GHO','2H12','2REC'};

sys_opt.data_dir = strcat(work_dir,'/PDB'); 
ws.sys.data_dir = sys_opt.data_dir;
ws.map.volumes =  GenerateSimulationMap.load_maps(ws);

[vols_large, max_siz] = VolumeUtil.resize_vols(ws.map.volumes, 1.0);

% definition of reconstruction model parameters
% set default values for SNR, missing wedge and ctf
ws.reconstruction_param.model = struct();
ws.reconstruction_param.model.SNR = 0.05;
ws.reconstruction_param.model.missing_wedge_angle = 30;
ws.reconstruction_param.model.ctf = GenerateSimulationMap.get_ctf_param(ws.map.map_resolution);
ws.reconstruction_param.model.ctf.voltage=300;

% set the contour level for each macromolecular complex
% calculate the minimum bounding spheres for each macromolecular complexes
contour_threshold = 0.2;
for i=1:numel(vols_large)
        density_threshold(i)=max(vols_large{i}(:))*contour_threshold;
        standard_sphere(i)=amprs_wthr(vols_large{i},density_threshold(i),max_siz(1));
        complex_volume(i)=sum(vols_large{i}(:)>density_threshold(i));
end

sphere_radius=zeros(numel(vols_large),1);
for i=1:numel(vols_large)
        sphere_radius(i)=standard_sphere(i).r;
end

freq_t = rand(numel(vols_large),1);

% here is the frequncy used in our study
%freq_t = [    
%    0.8147
%    0.9058
%    0.1270
%    0.9134
%    0.6324
%    0.0975
%    0.2785
%    0.5469
%    0.9575
%    0.9649
%    0.1576
%    0.9706
%    0.9572
%    0.4854
%    0.8003
%    0.1419
%    0.4218
%    0.9157
%    0.7922
%    0.9595
%    0.6557
%];


% vol_dim indicates the total number of particles in the simulation
% for example the first round we set 10*10*20=2000 particles and then
% the second round we set 10*10*50=5000 particles, as long as the 
% product of dim is consistant, the total number of particles for simulation
% will be the same, i.e., vol_dim=[10,10,20] is equal to vol_dim=[20,20,5] 
vol_dim=[1,1,4;10,10,50;10,10,80]; %bb get rid of 2nd and 3rd. The number of particles here must agree with the number of particles in the INPUT FILE.

disp('Initialization finished')

models=struct();
models=cell(3,1);
gss_opt=struct();
gss_opt=cell(3,1);
protein_complexes=struct();
protein_complexes=cell(3,1);
copy_number=zeros(numel(vols_large),3);

for j=1%:3 bbchange
	models{j}=struct();
	models{j}.vol_dim=vol_dim(j,:);
	models{j}.freq = freq_t / sum(freq_t);
	models{j}.class_lbls = MappingFrequencyMass_bb.generate_labels_via_frequency(prod(models{j}.vol_dim), models{j}.freq); %bb I changed MappingFrequencyMass to read the color values from the input file as euler angles.

	vols_grid_params = MappingFrequencyMass_bb.generate_model(models{j}, rots);
	gss_opt{j}.vols = vols_large;
	gss_opt{j}.grid_params = vols_grid_params;
	gss_opt{j}.reconstruction_param = ws.reconstruction_param;

	for i=1:numel(vols_large)
		copy_number(i,j)=length(find(models{j}.class_lbls==i));
	end

	% begins the simulation of each individual particles 
	tic
	disp('msb for protein complexes start');
	protein_complexes{j} = struct();
	protein_complexes{j} = cell(numel(gss_opt{j}.grid_params),1);
	for i=1:numel(gss_opt{j}.grid_params)
		protein_complexes{j}{i}=struct();
		disp('protein complex');	
		disp(i);
		protein_complexes{j}{i}.vols=VolumeUtil.rotate_vol_pad0(gss_opt{j}.vols{models{j}.class_lbls(i)},gss_opt{j}.grid_params{i}.rm,standard_sphere(models{j}.class_lbls(i)).c);
		protein_complexes{j}{i}.msb=standard_sphere(models{j}.class_lbls(i));
		protein_complexes{j}{i}.class_lbls=models{j}.class_lbls(i);
	end
	disp('msb for protein complexes finish');
	toc
end

% the radius of all particles will be IMP input
%dlmwrite('imp_radius2000',sphere_radius(models{1}.class_lbls)); bbchange
%dlmwrite('imp_radius5000',sphere_radius(models{2}.class_lbls)); bbchange
%dlmwrite('imp_radius8000',sphere_radius(models{3}.class_lbls)); bbchange

cd ../workspace

% bbchange save -v7.3 tomo_crowd_30.mat;







%%%%%%% tomo_densgenerate - simulate whole cell tomogram based on the location info generated by IMP and apply back projection to realistically simulate tomograms with factors like SNR, missing wedge and ctf
%%% needs: tomo_crowd.mat workspace file created by tomo_cplxgenerate.m - I think this includes the density maps of all of the proteins; 2k_ps_bb.pym workspace file created by tomo_locsgenerate.py - I think this can be edited manually - HOWEVER, I think it must agree with the tomo_crowd.mat file

% bbchange cd ../workspace
% bbchange load tomo_crowd_30; %bbchange - was tomo_crowd

% tomogram_size is the 3D volume to represent whole cell tomogram region
% in our study we set a 500*500*200 cubic 3D volume
% load the location file of spheres (particles) from the .pym file


fact=[];
for i=1:size(pym,1)/2
	fact(i,:)=pym(2*i,:);
end
cd ../code

% select different crowding levels, this whole part could be written in a loop form to
% set whole tomogram simulation for different crowding levels. here j=1 means the first crowding level (2000 particles)
j=1;

% to avoid edge effect of cutting a particle, we process the boundary
% by adding a sphere distance of the maximum sphere radius to the whole cell 3D boundaries 
fact(:,1:3)=fact(:,1:3)+ceil(max(sphere_radius));
cubic_size=tomogram_size+2*ceil(max(sphere_radius));

% this part is record each particle's category
for i=1:numel(vols_large)
	fact_sg{i}=fact(abs(fact(:,4)-sphere_radius(i))<1e-4,:);
end

fact_cb=[];
for i=1:numel(vols_large)
	fact_cb=[fact_cb;fact_sg{i}];
end

% permute all the particles by their category
orders=struct();
orders=[];
for i=1:numel(vols_large)
        category_index = find(models{j}.class_lbls==i);
        orders=[orders,category_index'];
end

% add each particle to its corresponding location to generate simulated tomogram
tic
vol_den=struct();
vol_den=zeros(cubic_size);
for i=1:size(fact_cb,1)
        disp('adding particle')
	disp(i)
	rotation_center = protein_complexes{j}{orders(i)}.msb.c;
        shifting_to_location = fact_cb(i,1:3);
        disp(gss_opt{j}.grid_params{orders(i)}.rm)
        vol_t_mut=VolumeUtil.rotate_vol_pad0(protein_complexes{j}{orders(i)}.vols, gss_opt{j}.grid_params{orders(i)}.rm,rotation_center,shifting_to_location,cubic_size,'cubic');
        vol_den=vol_den+vol_t_mut;
end
toc

% apply back projection to realistically simulate the tomograms with different SNR
vol_den_bp=struct();
vol_den_bp=cell(6,1);
SNR=[1] % bbchange, original: SNR=[100,50,20,10,5,1] %0.5,0.1,0.05,0.01];
for i=1:numel(SNR)
	disp('in the stage, where SNR is ')
	disp(SNR(i))
        vol_den_bp{i}=GenerateSimulationMap.backprojection_reconstruction(gss_opt{j}.reconstruction_param, vol_den, SNR(i));
end

clear vol_t_mut;
clear category_index;
clear rotation_center;
clear shifting_to_location;
clear i;

cd ../workspace

%save -v7.3 tomo_imp_2k_bb.mat; % bbchange was tomo_imp_2k.mat bbhelp - what is this file for? I guess it's just for saving the workspace

cd ../tomograms

% bbchange
% tom_mrcwrite(vol_den_bp{2},'name','tomo_2k_snr50.mrc','style','fei');
% tom_mrcwrite(vol_den_bp{3},'name','tomo_2k_snr20.mrc','style','fei');
% tom_mrcwrite(vol_den_bp{4},'name','tomo_2k_snr10.mrc','style','fei');
% tom_mrcwrite(vol_den_bp{6},'name','tomo_2k_snr1.mrc','style','fei');
tom_mrcwrite(vol_den_bp{1},'name','Alber_model_1d.mrc','style','fei');

cd ..; toc %bbchange

