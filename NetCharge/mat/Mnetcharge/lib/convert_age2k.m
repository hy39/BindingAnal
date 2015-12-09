%return a gene2protein table
function [] = convert_age2k(merged_dat)
     
    age_f = [0 1.1 4.4 8.1 12.3 16.9 22.3 28.5 35.7 44.2 54.4 66.6 81.4 99.3];

    dat = load(['dat/' merged_dat]);
    Viruses = dat.Viruses;
    gbacc = dat.gbacc;
    ages = dat.ages;
    iso_date_num = dat.iso_date_num;
    gbacc_all = dat.gbacc_all;
    
    infect_k = 0;
    for i1=1:length(Viruses(:,1))
        if Viruses(i1,1)==0
             ages(i1,2)=-1;
             Viruses(i1,6)=-1;
             continue;
        end
        for i2=1:length(age_f)-1
            if Viruses(i1,1)>age_f(i2) & Viruses(i1,1)<age_f(i2+1)
                infect_k = i2-1;
                ages(i1,2)=infect_k;
                Viruses(i1,6)=infect_k;
            end
        end
    end  
    
    save(['dat/' merged_dat], 'Viruses', 'gbacc', 'ages', 'iso_date_num', 'gbacc_all');
end