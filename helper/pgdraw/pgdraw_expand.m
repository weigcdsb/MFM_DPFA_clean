function pgSamp = pgdraw_expand(b,c)

b = ceil(b);
n = length(c);

c_expand = repelem(c, b);
grp_idx = repelem([1:n]', b);
pg_expand = pgdraw(c_expand);

pgSamp = accumarray(grp_idx,pg_expand);

end