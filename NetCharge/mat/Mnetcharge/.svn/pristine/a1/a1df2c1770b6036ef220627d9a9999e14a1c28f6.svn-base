%return a gene2protein table
function [gb] = create_geneprot_table(filename, gbacc, merged_dat)
    if ~exist('filename')
        geneprot(:,1) = textread('dat/h1n1/hm_h1n1_flu_noram_gene_access.txt','%s');
        geneprot(:,2) = textread('dat/h1n1/hm_h1n1_flu_noram_prot_access.txt','%s');
    else
        geneprot(:,1) = textread(['dat/' filename '_gene_access.txt'],'%s');
        geneprot(:,2) = textread(['dat/' filename '_prot_access.txt'],'%s');
    end
    
    gbacc_all = cellstr(gbacc);
    for i=1:length(gbacc_all)  
        %i = 403
        geneacc = geneprot(strcmp(geneprot(:,2),char(gbacc_all(i))),1);
        if isempty(geneacc)
            %gbacc_all(i,2);
            gbacc_all(i,2) = {''};
        else
            gbacc_all(i,2)= geneprot(strcmp(geneprot(:,2),char(gbacc_all(i))),1);
        end
    end
    
    dat = load(['dat/' merged_dat]);
    Viruses = dat.Viruses;
    gbacc = dat.gbacc;
    ages = dat. ages;
    iso_date_num = dat.iso_date_num;
    
    save(['dat/' merged_dat], 'Viruses', 'gbacc', 'ages', 'iso_date_num', 'gbacc_all');
end