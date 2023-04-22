#include "Stencil.h"



// ------------------------------- Stencil -------------------------------
void Stencil::SetCoefficientCount()
{ nCoefficients = ((nOrder + 1) / 2) * ((nOrder + 1) / 2); }
void Stencil::AllocateBuffers()
{
    w.resize(nDir);
    cx.resize(nDir);
    cy.resize(nDir);
    cz.resize(nDir);
    theta.resize(nDir);
    phi.resize(nDir);
    neighbour0.resize(nDir * 2*nDir);
    neighbour1.resize(nDir * 2*nDir);
    neighbour2.resize(nDir * 2*nDir);
}
void Stencil::SortDirections()
{
    int index[nDir];
    for(int i = 0; i < nDir; i++)
        index[i] = i;

    // sort index array based on custom compare function
    std::sort(index, index + nDir, [this](int i, int j)
    {
        // if(phi[i] != phi[j])
            // return phi[i] < phi[j];
        // else
            // return theta[i] < theta[j];
            
        if(theta[i] != theta[j])
            return theta[i] < theta[j];
        else
            return phi[i] < phi[j];
    });

    // create temporary arrays to hold sorted values
    double sorted_w[nDir];
    double sorted_cx[nDir];
    double sorted_cy[nDir];
    double sorted_cz[nDir];
    double sorted_theta[nDir];
    double sorted_phi[nDir];

    // copy values from original arrays to temporary arrays
    for(int i = 0; i < nDir; i++)
    {
        sorted_w[i] = w[index[i]];
        sorted_cx[i] = cx[index[i]];
        sorted_cy[i] = cy[index[i]];
        sorted_cz[i] = cz[index[i]];
        sorted_theta[i] = theta[index[i]];
        sorted_phi[i] = phi[index[i]];
    }

    // overwrite original arrays with sorted values
    for(int i = 0; i < nDir; i++)
    {
        w[i] = sorted_w[i];
        cx[i] = sorted_cx[i];
        cy[i] = sorted_cy[i];
        cz[i] = sorted_cz[i];
        theta[i] = sorted_theta[i];
        phi[i] = sorted_phi[i];
    }
}
void Stencil::InitializeConnectedTriangles()
{
    // Setup vertices vector for convex hull triangulation:
    std::vector<Vector3> vertices;
    vertices.reserve(nDir);
    for(int d=0; d<nDir; d++)
    {
        Vector3 v(Cx(d),Cy(d),Cz(d));
        vertices.push_back(v);
    }

    // Triangulate vertices:
    ConvexHull convexHull(vertices);
    convexHull.OriginalOrdering(vertices);

    // Extract connectedTriangles:
    std::vector<Vector3Int> allTriangles = convexHull.GetTriangles();
    for(int i=0; i<nDir; i++)
    {
        std::vector<Vector3Int> triangles;
        for(int j=0; j<allTriangles.size(); j++)
        {
            Vector3Int triangle = allTriangles[j];
            if(triangle[0] == i
            || triangle[1] == i
            || triangle[2] == i)
                triangles.push_back(triangle);
        }
        connectedTriangles.AddRow(triangles);
    }
}
void Stencil::InitializeConnectedVertices()
{
    // First order of connected vertices:
    for(int d=0; d<nDir; d++)
    {
        // Get unique indices of connected triangles:
        std::set<size_t> indices;
        for(int k=connectedTriangles.Start(d); k<connectedTriangles.End(d); k++)
        {
            Vector3Int triangle = connectedTriangles[k];
            indices.insert(triangle[0]);
            indices.insert(triangle[1]);
            indices.insert(triangle[2]);
        }

        // Extract corresponding vertices:
        std::vector<size_t> vertices;
        for(auto it : indices)
            vertices.push_back(it);

        // Add list of vertices to connected vertices:
        connectedVerticesOrder1.AddRow(vertices);
    }

    // Second order of connected vertices:
    for(int d=0; d<nDir; d++)
    {
        // Get unique indices of connected triangles:
        std::set<size_t> indices;
        for(int k=connectedVerticesOrder1.Start(d); k<connectedVerticesOrder1.End(d); k++)
        {
            size_t index = connectedVerticesOrder1[k];
            for(int p=connectedTriangles.Start(index); p<connectedTriangles.End(index); p++)
            {
                Vector3Int triangle = connectedTriangles[p];
                indices.insert(triangle[0]);
                indices.insert(triangle[1]);
                indices.insert(triangle[2]);
            }
        }

        // Extract corresponding vertices:
        std::vector<size_t> vertices;
        for(auto it : indices)
            vertices.push_back(it);

        // Add list of vertices to connected vertices:
        connectedVerticesOrder2.AddRow(vertices);
    }
}

double Stencil::W(size_t d) const
{ return w[d]; }
double Stencil::Theta(size_t d) const
{ return theta[d]; }
double Stencil::Phi(size_t d) const
{ return phi[d]; }
double Stencil::Cx(size_t d) const
{ return cx[d]; }
double Stencil::Cy(size_t d) const
{ return cy[d]; }
double Stencil::Cz(size_t d) const
{ return cz[d]; }
Tensor3 Stencil::C(size_t d) const
{ return Tensor3(Cx(d), Cy(d), Cz(d)); }

void Stencil::InitializeNearestNeighbourGrid()
{
    sphereGrid = SphereGrid(nDir, 2*nDir);

    for(size_t j=0; j<sphereGrid.nPh; j++)
    for(size_t i=0; i<sphereGrid.nTh; i++)
    {
        Tensor3 c = sphereGrid.C(i,j);
        double dot0 = -1;
        double dot1 = -1;
        double dot2 = -1;
        int index0 = -1;
        int index1 = -1;
        int index2 = -1;
        for(size_t d=0; d<nDir; d++)
        {
            double dot = Tensor3::Dot(c,C(d));
            if(dot > dot0)
            {
                dot2 = dot1;
                dot1 = dot0;
                dot0 = dot;
                index2 = index1;
                index1 = index0;
                index0 = d;
            }
            else if(dot > dot1)
            {
                dot2 = dot1;
                dot1 = dot;
                index2 = index1;
                index1 = d;
            }
            else if(dot > dot2)
            {
                dot2 = dot;
                index2 = d;
            }
        }
        neighbour0[sphereGrid.Index(i,j)] = index0;
        neighbour1[sphereGrid.Index(i,j)] = index1;
        neighbour2[sphereGrid.Index(i,j)] = index2;
    }
}
size_t Stencil::GetNeighbourIndex0(const Tensor3& p)
{
    size_t i = std::round(sphereGrid.i(p.Theta()));
    size_t j = ((int)std::round(sphereGrid.j(p.Phi())) + sphereGrid.nPh) % sphereGrid.nPh;
    return neighbour0[sphereGrid.Index(i,j)];
}
size_t Stencil::GetNeighbourIndex1(const Tensor3& p)
{
    size_t i = std::round(sphereGrid.i(p.Theta()));
    size_t j = ((int)std::round(sphereGrid.j(p.Phi())) + sphereGrid.nPh) % sphereGrid.nPh;
    return neighbour1[sphereGrid.Index(i,j)];
}
size_t Stencil::GetNeighbourIndex2(const Tensor3& p)
{
    size_t i = std::round(sphereGrid.i(p.Theta()));
    size_t j = ((int)std::round(sphereGrid.j(p.Phi())) + sphereGrid.nPh) % sphereGrid.nPh;
    return neighbour2[sphereGrid.Index(i,j)];
}
Tensor3 Stencil::GetNeighbour0(const Tensor3& p)
{
    size_t i = std::round(sphereGrid.i(p.Theta()));
    size_t j = ((int)std::round(sphereGrid.j(p.Phi())) + sphereGrid.nPh) % sphereGrid.nPh;
    size_t d = neighbour0[sphereGrid.Index(i,j)];
    return C(d);
}
Tensor3 Stencil::GetNeighbour1(const Tensor3& p)
{
    size_t i = std::round(sphereGrid.i(p.Theta()));
    size_t j = ((int)std::round(sphereGrid.j(p.Phi())) + sphereGrid.nPh) % sphereGrid.nPh;
    size_t d = neighbour1[sphereGrid.Index(i,j)];
    return C(d);
}
Tensor3 Stencil::GetNeighbour2(const Tensor3& p)
{
    size_t i = std::round(sphereGrid.i(p.Theta()));
    size_t j = ((int)std::round(sphereGrid.j(p.Phi())) + sphereGrid.nPh) % sphereGrid.nPh;
    size_t d = neighbour2[sphereGrid.Index(i,j)];
    return C(d);
}

void Stencil::Print() const
{
    std::cout << "        d,\t        w,\t    theta,\t      phi,\t       cx,\t       cy,\t       cz\n";
    for(int d=0; d<nDir; d++)
    {
        std::cout << Format(d) << ",\t" << Format(W(d)) << ",\t" << Format(Theta(d)) << ",\t" << Format(Phi(d)) << ",\t";
        std::cout << Format(Cx(d)) << ",\t" << Format(Cy(d)) << ",\t" << Format(Cz(d)) << "\n";
    }
    std::cout << std::endl;
}
// -----------------------------------------------------------------------



// -------------------------- LebedevStencil ---------------------------
LebedevStencil::LebedevStencil(size_t nOrder)
{
    name = "Lebedev" + std::to_string(nOrder);
    this->nOrder = nOrder;
    SetCoefficientCount();

    switch (nOrder)
    {
        case  3:
                 #include "../stencils/LebedevStencil/LebedevStencil3" 
                 break;
        case  5:
                 #include "../stencils/LebedevStencil/LebedevStencil5" 
                 break;
        case  7:
                 #include "../stencils/LebedevStencil/LebedevStencil7" 
                 break;
        case  9:
                 #include "../stencils/LebedevStencil/LebedevStencil9" 
                 break;
        case 11:
                 #include "../stencils/LebedevStencil/LebedevStencil11"
                 break;
        case 13:
                 #include "../stencils/LebedevStencil/LebedevStencil13"
                 break;
        case 15:
                 #include "../stencils/LebedevStencil/LebedevStencil15"
                 break;
        case 17:
                 #include "../stencils/LebedevStencil/LebedevStencil17"
                 break;
        case 19:
                 #include "../stencils/LebedevStencil/LebedevStencil19"
                 break;
        case 21:
                 #include "../stencils/LebedevStencil/LebedevStencil21"
                 break;
        case 23:
                 #include "../stencils/LebedevStencil/LebedevStencil23"
                 break;
        case 25:
                 #include "../stencils/LebedevStencil/LebedevStencil25"
                 break;
        case 27:
                 #include "../stencils/LebedevStencil/LebedevStencil27"
                 break;
        case 29:
                 #include "../stencils/LebedevStencil/LebedevStencil29"
                 break;
        case 31:
                 #include "../stencils/LebedevStencil/LebedevStencil31"
                 break;
        case 35:
                 #include "../stencils/LebedevStencil/LebedevStencil35"
                 break;
        case 41:
                 #include "../stencils/LebedevStencil/LebedevStencil41"
                 break;
        case 47:
                 #include "../stencils/LebedevStencil/LebedevStencil47"
                 break;
        default:
            ExitOnError("Invalid LebedevStencil nOrder. See stencil/LebedevStencil for valid orders.");
            break;
    }
    
    SortDirections();
    InitializeConnectedTriangles();
    InitializeConnectedVertices();
    InitializeNearestNeighbourGrid();
}
// ---------------------------------------------------------------------



// ----------------------- GaussLegendreStencil ------------------------
GaussLegendreStencil::GaussLegendreStencil(size_t nOrder)
{
    name = "GaussLegendre" + std::to_string(nOrder);
    this->nOrder = nOrder;
    SetCoefficientCount();

    switch (nOrder)
    {
        case  3:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil3" 
                 break;
        case  5:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil5" 
                 break;
        case  7:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil7" 
                 break;
        case  9:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil9" 
                 break;
        case 11:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil11"
                 break;
        case 13:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil13"
                 break;
        case 15:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil15"
                 break;
        case 17:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil17"
                 break;
        case 19:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil19"
                 break;
        case 21:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil21"
                 break;
        case 23:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil23"
                 break;
        case 25:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil25"
                 break;
        case 27:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil27"
                 break;
        case 29:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil29"
                 break;
        case 31:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil31"
                 break;
        case 33:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil33"
                 break;
        case 35:
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil35"
                 break;
        default:
            ExitOnError("Invalid GaussLegendreStencil nOrder. Must be odd and smaller 36.");
            break;
    }
    
    SortDirections();
    InitializeConnectedTriangles();
    InitializeConnectedVertices();
    InitializeNearestNeighbourGrid();
}
// ---------------------------------------------------------------------