# An octave(matlab) program to solve Einstein Field Equations in 1 dimension (spherical reduction)

* Working parameters: 
    * inB = 1.5, outB = 2.5;
    * N = 15, meshNum = 1;
    * paragamma0 = 10;
    * bdry\_type = FREEZING.

* Want to make working: 
    * outB >= 5;
    * bdry_type = FREEZING or DIRICHLET;
    * paragamma0?

* How to combine freezing incoming characteristic fields condition with filter?
    * If we apply freezing condition first, then applying the filter will introduce instabilities(known from tests);
    * Adjusting U(vmapO) will actually affect the deriU in outermost cell, even though all the nodes' value(except vmapO) are the same;
    *   
