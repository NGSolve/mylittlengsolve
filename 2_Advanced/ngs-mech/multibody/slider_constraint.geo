#
## Kurbelwelle + Pleuel
#
algebraic3d

# cut cylinder by planes:

solid p1 = plane (-5, 0, 0; -1, 0, 0) -bc=1;
solid p2 = plane (35, 0, 0; 1, 0, 0) -bc=2;
solid p3 = cylinder (9, 0, 20; 21, 0, 20; 3.5) -bc=3;

solid shaft1 = cylinder ( -5, 0, 0; 1, 0, 0; 5)
	and p1
	and plane (1, 0, 0; 10, 0, 0);

solid shaft2 = 
	    plane (0, 0, 0; -1, 0, 0)
	and plane (10, 0, 0; 1, 0, 0)
	and (cylinder (0, 0., 0.; 10, 0., 0.; 11)
	     or cylinder (0, 0, 20; 10, 0, 20; 6)
	     or (    plane (0, 5, 0; 0, 1, 0)
	    	 and plane (0, -5, 0; 0, -1, 0)
		 and plane (0, 0, 0; 0, 0, -1)
		 and plane (0, 0, 20; 0, 0, 1)));  

solid shaft3 = p3
	and plane (9, 0, 0; -1, 0, 0)
	and plane (21, 0, 0; 1, 0, 0);

solid shaft4 = 
	    plane (20, 0, 0; -1, 0, 0)
	and plane (30, 0, 0; 1, 0, 0)
	and (cylinder (20, 0., 0.; 30, 0., 0.; 11)
	     or cylinder (20, 0, 20; 30, 0, 20; 6)
	     or (    plane (0, 5, 0; 0, 1, 0)
	    	 and plane (0, -5, 0; 0, -1, 0)
		 and plane (0, 0, 0; 0, 0, -1)
		 and plane (0, 0, 20; 0, 0, 1)));  

solid shaft5 = cylinder ( 29, 0, 0; 35, 0, 0; 5)
	and plane (29, 0, 0; -1, 0, 0)
	and p2;
	
solid shaft = shaft1 or shaft2 or shaft3 or shaft4 or shaft5 -bc=11  -maxh=50;


solid p4 = cylinder (12, 0, 20; 18, 0, 20; 3.6) -bc=4;
solid p5 = cylinder (12, 0, 80; 18, 0, 80; 3.6) -bc=5;
solid p9 = plane (0, -4, 24; 0, -1, 0) -bc=9;
solid p10 = plane (0,  4, 24; 0,  1, 0) -bc=10;

solid pl1 = cylinder (10.1, 0, 20; 19.9, 0, 20; 7)
	and not p4
	and plane (10.1, 0, 0; -1, 0, 0)
	and plane (19.9, 0, 0;  1, 0, 0);

solid pl2 = cylinder (9.1, 0, 80; 20.9, 0, 80; 7)
	and not p5
	and plane (9.1, 0, 0; -1, 0, 0)
	and plane (20.9, 0, 0;  1, 0, 0);

solid pl1x = plane (13.5, 0, 0; -1, 0, 0) -maxh=4;
solid pl2x = plane (16.5, 0, 0;  1, 0, 0) -maxh=4;

solid pl3 = p9
	and p10
	and plane (0,  4, 24.5; 0,  0,-1)
	and plane (0,  4, 75.5; 0,  0, 1)
	and pl1x
	and pl2x;
	
solid pleuel = (pl1 or pl2 or pl3) -bc=11;


solid p6 = cylinder (-1, 0, 80; 31, 0, 80; 3.5) -bc=6;
solid p7 = plane (0., -15, 70; 0, -1, 0) -bc=7;
solid p8 = plane (30., 15, 95; 0, 1, 0) -bc=8;
solid p12 = plane (30., 15, 95;  0, 0, 1) -bc=12;

solid cub1 = 
	    plane (0., -15, 70; -1, 0, 0)
	and p7
	and plane (0., -15, 70; 0, 0, -1)
	and plane (30., 15, 95;  1, 0, 0)
	and p8
	and p12;
	
solid cub2 = 
	    plane (9, -12, 60; -1, 0, 0)
	and plane (9, -12, 60; 0, -1, 0)
	and plane (9, -12, 60; 0, 0, -1)
	and plane (21, 12, 90;  1, 0, 0)
	and plane (21, 12, 90;  0, 1, 0)
	and plane (21, 12, 90;  0, 0, 1);
	
solid cyl3 = p6
	and plane (1, 0, 80; -1, 0, 0)
	and plane (29, 0, 80; 1, 0, 0);

solid cyl = cub1 and not cub2 or cyl3 -bc=11 -maxh=25;


tlo shaft -transparent;
tlo pleuel -col=[0,1,0];
tlo cyl -col=[1,0,0];

