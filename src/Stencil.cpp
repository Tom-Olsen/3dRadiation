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
    double sorted_phi[nDir];
    double sorted_theta[nDir];

    // copy values from original arrays to temporary arrays
    for(int i = 0; i < nDir; i++)
    {
        sorted_w[i] = w[index[i]];
        sorted_cx[i] = cx[index[i]];
        sorted_cy[i] = cy[index[i]];
        sorted_cz[i] = cz[index[i]];
        sorted_phi[i] = phi[index[i]];
        sorted_theta[i] = theta[index[i]];
    }

    // overwrite original arrays with sorted values
    for(int i = 0; i < nDir; i++)
    {
        w[i] = sorted_w[i];
        cx[i] = sorted_cx[i];
        cy[i] = sorted_cy[i];
        cz[i] = sorted_cz[i];
        phi[i] = sorted_phi[i];
        theta[i] = sorted_theta[i];
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
    for(int d=0; d<nDir; d++)
    {
        // Get unique indices of connected triangles:
        std::set<size_t> indices;
        int start = connectedTriangles.Start(d);
        int end = connectedTriangles.End(d);
        for(int k=start; k<end; k++)
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
        connectedVertices.AddRow(vertices);
    }
}
void Stencil::InitializeMinMaxDot()
{
    minMaxDot = 1;
    int n = 500;
    for(int i=0; i<n; i++)
    for(int j=0; j<2*n; j++)
    {
        double theta = M_PI * (i + 0.5) / (double)n;
        double phi = 2.0 * M_PI * j / (double)n;
        Tensor3 c(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
        double maxDot = -1;
        for(int d=0; d<nDir; d++)
        {
            double dot = Tensor3::Dot(C(d), c);
            maxDot = std::max(maxDot,dot);
        }
        minMaxDot = std::min(minMaxDot,maxDot);
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

void Stencil::Print() const
{
    std::cout << "        d\t        w\t    theta\t      phi\t       cx\t       cy\t       cz\n";
    for(int d=0; d<nDir; d++)
    {
        std::cout << Format(d) << "\t" << Format(W(d)) << "\t" << Format(Theta(d)) << "\t" << Format(Phi(d)) << "\t";
        std::cout << Format(Cx(d)) << "\t" << Format(Cy(d)) << "\t" << Format(Cz(d)) << "\n";
    }
    std::cout << std::endl;
}
// -----------------------------------------------------------------------



// ------------------------------ MyStencil ------------------------------
MyStencil::MyStencil(size_t nOrder)
{
    nTh = nOrder;
    nPh = 2 * nOrder;
    this->nOrder = nOrder;
    this->nDir = nTh * nPh;
    SetCoefficientCount();

    switch(nOrder)
    {
        case  3: dTheta = 0.88607712379261370; break;
        case  5: dTheta = 0.56708068878468730; break;
        case  7: dTheta = 0.41679707883494480; break;
        case  9: dTheta = 0.32944347754574150; break;
        case 11: dTheta = 0.27234940787623110; break;
        case 13: dTheta = 0.23211697811143664; break;
        case 15: dTheta = 0.20223903140287983; break;
        case 17: dTheta = 0.17917455393695525; break;
        case 19: dTheta = 0.16083171772074056; break;
        default:
            exit_on_error("MyStencil, invalid nOrder.");
    }
    
    InitializeConnectedTriangles();
    InitializeConnectedVertices();
    InitializeMinMaxDot();

    // Initialize weights:
    w0 = 0;
    for(size_t d1=0; d1<nPh; d1++)
    for(size_t d0=0; d0<nTh; d0++)
        w0 += MySin(Theta(d0,d1));
    w0 = 4.0 * M_PI / w0;
}
// Overrides:
double MyStencil::W(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return W(d0,d1);
}
double MyStencil::Theta(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return Theta(d0,d1);
}
double MyStencil::Phi(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return Phi(d0,d1);
}
double MyStencil::Cx(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return Cx(d0,d1);
}
double MyStencil::Cy(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return Cy(d0,d1);
}
double MyStencil::Cz(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return Cz(d0,d1);
}
Tensor3 MyStencil::C(size_t d) const
{
    int d0 = d % nTh;
    int d1 = d / nTh;
    return Tensor3(Cx(d0,d1), Cy(d0,d1), Cz(d0,d1));
}

// Indexing:
size_t MyStencil::Index(size_t d0, size_t d1) const
{ return d0 + d1 * nTh; }
double MyStencil::d0(double theta) const
{ return nTh / 2 + (theta - M_PI_2) / dTheta; }
double MyStencil::d1(double phi) const
{ return phi * nPh / (2.0 * M_PI) - 0.5; }

// Acces with two indices:
double MyStencil::W(size_t d0, size_t d1) const
{ return w0 * MySin(Theta(d0,d1)); }

double MyStencil::Theta(double d0, double d1) const
{ return M_PI_2 + (d0 - nTh / 2) * dTheta; }
double MyStencil::Phi(double d0, double d1) const
{ return 2.0 * M_PI * (d1 + 0.5) / nPh; }

double MyStencil::Cx(double d0, double d1) const
{ return MySin(Theta(d0,d1)) * MyCos(Phi(d0,d1)); }
double MyStencil::Cy(double d0, double d1) const
{ return MySin(Theta(d0,d1)) * MySin(Phi(d0,d1)); }
double MyStencil::Cz(double d0, double d1) const
{ return MyCos(Theta(d0,d1)); }
Tensor3 MyStencil::C(double d0, double d1) const
{ return Tensor3(Cx(d0,d1), Cy(d0,d1), Cz(d0,d1)); }
// ---------------------------------------------------------------------



// -------------------------- LebedevStencil ---------------------------
LebedevStencil::LebedevStencil(size_t nOrder)
{
    this->nOrder = nOrder;
    SetCoefficientCount();

    switch (nOrder)
    {
        case  3: nDir =   6; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil3" 
                 break;
        case  5: nDir =  14; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil5" 
                 break;
        case  7: nDir =  26; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil7" 
                 break;
        case  9: nDir =  38; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil9" 
                 break;
        case 11: nDir =  50; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil11"
                 break;
        case 13: nDir =  74; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil13"
                 break;
        case 15: nDir =  86; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil15"
                 break;
        case 17: nDir = 110; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil17"
                 break;
        case 19: nDir = 146; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil19"
                 break;
        case 21: nDir = 170; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil21"
                 break;
        case 23: nDir = 194; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil23"
                 break;
        case 25: nDir = 230; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil25"
                 break;
        case 27: nDir = 266; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil27"
                 break;
        case 29: nDir = 302; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil29"
                 break;
        case 31: nDir = 350; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil31"
                 break;
        case 35: nDir = 434; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil35"
                 break;
        case 41: nDir = 590; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil41"
                 break;
        case 47: nDir = 770; AllocateBuffers();
                 #include "../stencils/LebedevStencil/LebedevStencil47"
                 break;
        default:
            exit_on_error("Invalid LebedevStencil nOrder. See stencil/LebedevStencil for valid orders.");
            break;
    }
    
    SortDirections();
    InitializeConnectedTriangles();
    InitializeConnectedVertices();
    InitializeMinMaxDot();
}
// ---------------------------------------------------------------------



// ----------------------- GaussLegendreStencil ------------------------
GaussLegendreStencil::GaussLegendreStencil(size_t nOrder)
{
    this->nOrder = nOrder;
    this->nDir = nOrder * (2 * nOrder);
    SetCoefficientCount();

    switch (nOrder)
    {
        case  3: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil3" 
                 break;
        case  5: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil5" 
                 break;
        case  7: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil7" 
                 break;
        case  9: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil9" 
                 break;
        case 11: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil11"
                 break;
        case 13: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil13"
                 break;
        case 15: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil15"
                 break;
        case 17: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil17"
                 break;
        case 19: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil19"
                 break;
        case 21: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil21"
                 break;
        case 23: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil23"
                 break;
        case 25: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil25"
                 break;
        case 27: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil27"
                 break;
        case 29: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil29"
                 break;
        case 31: AllocateBuffers();
                 #include "../stencils/GaussLegendreStencil/GaussLegendreStencil31"
                 break;
        default:
            exit_on_error("Invalid GaussLegendreStencil nOrder. Must be odd and smaller 32.");
            break;
    }
    
    SortDirections();
    InitializeConnectedTriangles();
    InitializeConnectedVertices();
    InitializeMinMaxDot();
}
// ---------------------------------------------------------------------