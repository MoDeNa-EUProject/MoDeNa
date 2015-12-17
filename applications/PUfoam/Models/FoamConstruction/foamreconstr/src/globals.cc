namespace globals {
    bool createNodes; //create struts at cell vertices
    bool createEdges; //create struts at cell edges
    bool openCell; //make open cell foam
    bool save_dat; //save file in old dx style
    bool save_vtk; //save file in new paraview style
    double dstrut; //parameter influencing size of struts in cell vertices
    double dedge; //parameter influencing size of struts in cell edges
    double strutPorosity; //desired porosity of only struts
    double RANDOM; //number between 0.0 and 1.0
    //how much are positions of cell centers perturbated from even position
    int grid; //grid type for cell centers, 1=cubic grid,
    //2=hexagonal lattice, ABAB...,3=hexagonal lattice, ABCABC...
    int nx; //domain size in X
    int ny; //domain size in Y
    int nz; //domain size in Z
    int sx; //size of cells in X
    int sy; //size of cells in Y
    int sz; //size of cells in Z
    bool save_voro_diag1; //gnuplot Voronoi diagram
    bool save_voro_diag2; //alternative gnuplot Voronoi diagram
    bool import_vtk; //import morphology from vtk
    bool progress_report; //show detailed progress report
}
