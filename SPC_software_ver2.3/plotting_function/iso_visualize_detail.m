function [p] = iso_visualize_detail(Z,iso,xx,yy,zz);

[d1 d2 d3] = size(Z);
[x, y, z] = meshgrid (xx, yy, zz);

view (120, 30);
p = patch(isosurface (x, y, z, Z, iso));
isonormals(x,y,z, Z, p)
set(p, 'FaceColor', 'green', 'EdgeColor', 'none');
daspect([1 1 1]); axis tight; 
colormap(prism(28))
camlight; lighting phong

