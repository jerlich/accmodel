function err=create_temp_table(table,field,type,data)





sqlstr=['create temporary table `' table '` (`' field '` ' type ') engine=myisam'];


bdata(sqlstr);

if isnumeric(data)
    for dx=1:numel(data)
        bdata('insert into "{S}" ("{S}") values ("{S}")',table, field, data(dx));
    end
elseif iscell(data)
    
    for dx=1:numel(data)
        bdata('insert into "{S}" ("{S}") values ("{S}")',table, field, data{dx});
    end
end
