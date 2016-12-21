% MappingFrequencyMass.m
% functions for handling frequency assignment and instances numbers for generating mixture of macromplecular complexes

classdef MappingFrequencyMass_bb

    methods (Static = true)
        
	%The grid_params represents each individual instance of particles (macromolecular complex)
        function grid_params = generate_model(model, rots)            
            grid_params = cell(model.vol_dim);
            for vol_i = 1 : numel(model.class_lbls)
                grid_params{vol_i} = struct();
                grid_params{vol_i}.class_lbls = model.class_lbls(vol_i);
		%class_lbls will indicate the type of macromolecular complex the instance (starting from type 1 to 21 in our study)
                disp(rots(vol_i,:))
                grid_params{vol_i}.rm = AngLocUtil.rotation_matrix_zyz(rots(vol_i,:) * pi); %bbchange - was AngLocUtil.rotation_matrix_zyz(rand(1, 3) * pi)
		%rm is the rotation matrix. different instances of particles from the same type of macromolecular complex could illustrate different orientations by applying the corresponding rotational matrix
            end            
        end     
                        
        function lbls = generate_labels_via_frequency(element_num, freq)
            assert(abs( sum(freq) - 1 ) < 1e-5);            
	    freq_cdf = zeros(size(freq));      
	    freq_cdf(1) = freq(1);
	    %calculate the cumulated density function for all macromolecular complexes
            for lbl = 2 : numel(freq)
                freq_cdf(lbl) = freq_cdf(lbl-1) + freq(lbl);
            end
            freq_cdf(numel(freq_cdf)) = 1;      
	    % just to avoid numerical imprecision
            freq_cdf_num = freq_cdf * element_num;
            
            lbls = zeros(element_num, 1);    
	    lbl_current = 1;
            for element_i = 1 : element_num
                while element_i > freq_cdf_num(lbl_current)
                    lbl_current = lbl_current + 1;
                end
                
                lbls(element_i) = lbl_current;
            end
            [r_s, r_i] = sort(rand(numel(lbls),1));      
	    lbls = lbls(r_i);         % shuffle
            if false
                for lbl = 1 : numel(freq)
                    disp(abs(sum(lbls == lbl) / numel(lbls)  - freq(lbl)));
                end
            end            
        end     % generate_labels_via_frequency()                        
    end     % methods (Static = true)
end     

