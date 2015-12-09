function outlier = removeOutlier( file_outlier, file_fasta )
% Remove outlier
% > removeOutlier('dat/hm_h3n2_noram_1968_2013/hm_h3n2_noram_outlier_sets.csv', 'dat/hm_h3n2_noram_1968_2013/hm_h3n2_noram_cds.fas');

    fid = fopen(file_outlier);
    c = textscan(fid, '%s', 'Delimiter', ',');
    outlier = c{1};
    
    [Header Sequences] = fastaread(file_fasta);
    rm_id = []; %list of ID that should be removed 
    for i=1:length(Header)
        for j=1:length(outlier)
            %TF = strcmp(outlier, Header(i));
            pos = cell2mat(regexp(Header(i),outlier(j)))
            if ~isempty(pos)
                rm_id(end+1) = i;
                disp 'remove outlier';
            end
        end
    end
    rm_id;
    Header(rm_id) = [];
    Sequences(rm_id) = [];
    
    % Save fasta
    pos_s = regexp(file_fasta,'\.');
    fname = [file_fasta(1:pos_s-1) '_std'];
    fastawrite([fname '.fas'], Header, Sequences); 
end
