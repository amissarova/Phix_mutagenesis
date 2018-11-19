function AA = getAAfromCodon(Codon , AA_table)
%%
idx = find(strcmp(AA_table.Codon , Codon));
AA = AA_table.AA{idx};

end