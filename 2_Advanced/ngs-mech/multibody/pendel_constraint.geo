#
## a pendel
#
algebraic3d

solid cube = orthobrick (0.06, -0.05, -0.05; 1, 0.05, 0.05);

# cut cylinder by planes:

solid lager = cylinder ( 0, -0.08, 0; 0, 0.08, 0; 0.05 ) -bc=1;

solid fincyl1 = cylinder ( 0, -0.08, 0; 0, 0.08, 0; 0.1 )
	and plane (0, -0.08, 0; 0, -1, 0)
	and plane (0, 0.08, 0; 0, 1, 0)
	and not lager;
	
solid pendel= cube or fincyl1 -bc=2;

tlo pendel;

