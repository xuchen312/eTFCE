#include <Rcpp.h>
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <random>
#include <chrono>

//------------------------- (0) PREPARATIONS FOR 2D INPUTS -------------------------//

// Find valid neighbours of a pixel
std::vector<int> findNeighbours2D(Rcpp::IntegerVector& MASK,   // 2D mask of original orders (1:m)
                                  Rcpp::IntegerVector& DIMS,   // image dimensions
                                  int index,                   // 0-based pixel index of interest
                                  int conn)                    // connectivity criterion
{
    // compute [x y] coordinates based on pixel index
    int x = index % DIMS[0] + 1;
    int y = index / DIMS[0] + 1;
    
    // xy coordinate adjustment vectors
    int DX[8] = {1,-1,0, 0,  1,-1, 1,-1};
    int DY[8] = {0, 0,1,-1,  1, 1,-1,-1};
    
    // find all valid neighbours of a pixel
    std::vector<int> IDS;
    IDS.reserve(conn);
    for (int i = 0; i < conn; i++)
    {
        int nx = x + DX[i];
        int ny = y + DY[i];
        int id = (ny-1)*DIMS[0] + (nx-1);
        
        if (nx >= 1 && nx <= DIMS[0] &&
            ny >= 1 && ny <= DIMS[1] &&
            MASK[id] != 0)
        {
            IDS.push_back(MASK[id]);
        }
    }
    
    return IDS;
}

// Find the adjacency list for all in-mask voxels
std::vector<std::vector<int>> findAdjList2D(Rcpp::IntegerVector& MASK,     // 2D mask of original orders (1:m)
                                            Rcpp::IntegerVector& DIMS,     // image dimensions
                                            std::vector<int>&    INDEXP,   // voxel indices of unsorted p-values
                                            int m,                         // number of in-mask voxels
                                            int conn)                      // connectivity criterion
{
    std::vector<std::vector<int>> ADJ;
    ADJ.reserve(m);
    for (int i = 0; i < m; i++)
    {
        ADJ.emplace_back(findNeighbours2D(MASK, DIMS, INDEXP[i], conn));
    }
    
    return ADJ;
}


//------------------------- (1) PREPARATIONS FOR 3D INPUTS -------------------------//

// Find valid neighbours of a voxel
std::vector<int> findNeighbours(Rcpp::IntegerVector& MASK,   // 3D mask of original orders (1:m)
                                Rcpp::IntegerVector& DIMS,   // image dimensions
                                int index,                   // voxel index of interest
                                int conn)                    // connectivity criterion
{
    // compute [x y z] coordinates based on voxel index
    int x = index % DIMS[0] + 1;
    int y = (index/DIMS[0]) % DIMS[1] + 1;
    int z = index/(DIMS[0]*DIMS[1]) + 1;
    
    // xyz coordinate adjustment vectors
    int DX[26] = {1,-1,0, 0,0, 0,  1,-1, 1,-1,1,-1, 1,-1,0, 0, 0, 0,  1,-1, 1,-1, 1,-1, 1,-1};
    int DY[26] = {0, 0,1,-1,0, 0,  1, 1,-1,-1,0, 0, 0, 0,1,-1, 1,-1,  1, 1,-1,-1, 1, 1,-1,-1};
    int DZ[26] = {0, 0,0, 0,1,-1,  0, 0, 0, 0,1, 1,-1,-1,1, 1,-1,-1,  1, 1, 1, 1,-1,-1,-1,-1};
    
    // find all valid neighbours of a voxel
    std::vector<int> IDS;
    IDS.reserve(conn);
    for (int i = 0; i < conn; i++)
    {
        int nx = x + DX[i];
        int ny = y + DY[i];
        int nz = z + DZ[i];
        int id = (nz-1)*DIMS[0]*DIMS[1] + (ny-1)*DIMS[0] + (nx-1);
        
        if (nx >= 1 && nx <= DIMS[0] &&
            ny >= 1 && ny <= DIMS[1] &&
            nz >= 1 && nz <= DIMS[2] &&
            MASK[id] != 0)
        {
            IDS.push_back(MASK[id]);
        }
    }
    
    return IDS;
}

// Find the adjacency list for all in-mask voxels
std::vector<std::vector<int>> findAdjList(Rcpp::IntegerVector& MASK,     // 3D mask of original orders (1:m)
                                          Rcpp::IntegerVector& DIMS,     // image dimensions
                                          std::vector<int>&    INDEXP,   // voxel indices of unsorted p-values
                                          int m,                         // number of in-mask voxels
                                          int conn)                      // connectivity criterion
{
    std::vector<std::vector<int>> ADJ;
    ADJ.reserve(m);
    for (int i = 0; i < m; i++)
    {
        ADJ.emplace_back(findNeighbours(MASK, DIMS, INDEXP[i], conn));
    }
    
    return ADJ;
}


//------------------------- (2) FIND ALL STCS (USING SORTING RANKS) -------------------------//

// Find function for disjoint set data structure
// (more stable iterative version)
int Find(int x,
         std::vector<int>& PARENT)
{
    while (PARENT[x] != x)
    {
        PARENT[x] = PARENT[PARENT[x]];
        x         = PARENT[x];
    }
    
    return x;
}

// Union function of a disjoint-set data structure based on the "union by size" technique.
// The disjoint sets will represent the components of a forest, and the data structure is
// augmented to keep track of the forest root of each component. UnionBySize(i,j) merges
// the sets S_i and S_j (where S_k denotes the set containing k), and assigns the forestroot
// of S_i to be the forestroot of the resulting union.
void UnionBySize(int i,
                 int j,
                 std::vector<int>&    PARENT,
                 std::vector<int>&    FORESTROOT,
                 std::vector<int>&    SIZEchange,
                 std::vector<int>&    SIZE,
                 std::vector<double>& MASS)
{
    int irep = Find(i, PARENT);  // root in union-find data structure
    int jrep = Find(j, PARENT);
    
    // if i and j are already in the same set
    if (irep == jrep) return;
    
    // i and j are not in same set, so we merge
    int iroot = FORESTROOT[irep];  // real root in directed rooted forest!
    int jroot = FORESTROOT[jrep];
    if (SIZE[iroot] < SIZE[jroot])
    {
        PARENT[irep] = jrep;
        FORESTROOT[jrep] = iroot;
    }
    else
    {
        PARENT[jrep] = irep;
    }
    SIZE[iroot] += SIZE[jroot];
    MASS[iroot] += MASS[jroot];
    
    // record the parent of jroot
    if (SIZEchange[jroot] == jroot) SIZEchange[jroot] = iroot;
}

// Compute all supra-threshold clusters (STCs)
struct ClusterResults
{
    std::vector<int> SIZEchange;
    std::vector<int> SIZE;
    std::vector<double> MASS;
};
ClusterResults findClusters(int m,                                 // number of nodes
                            std::vector<std::vector<int>>& ADJ,    // a list of neighbours for all nodes (unsorted!)
                            std::vector<double>&           STATS,  // unsorted statistics
                            std::vector<int>&              ORD,    // sorted orders for non-ascending stats
                            std::vector<int>&              RANK)   // sorting ranks
{
    // initialize output (1): a vector of each voxel's parents in the tree
    std::vector<int> SIZEchange(m);
    // initialize output (2): a vector of subtree size
    std::vector<int> SIZE(m, 1);
    // initialize output (3): a vector of subtree mass
    std::vector<double> MASS = STATS;
    
    // initialize forest roots
    std::vector<int> ROOT(m, -1);
    
    // prepare disjoint set data structure
    std::vector<int> PARENT(m);
    std::vector<int> FORESTROOT(m);
    for (int i = 0; i < m; i++)
    {
        PARENT[i] = i;
        FORESTROOT[i] = i;
        SIZEchange[i] = i;
    }
    
    // loop through all nodes in the non-ascending order of stats
    // // for C++11 and 0-based ORD, the use of v is not needed
    // for(int i : ORD)
    for (int i = 0; i < m; i++)
    {
        int v = ORD[i];
        
        // find neighbours for node with the ith largest stat
        const std::vector<int>& IDS = ADJ[v];
        
        // loop through all its neighbours
        for (int j = 0; j < IDS.size(); j++)
        {
            int nb = IDS[j] - 1;
            
            if (RANK[nb] < i+1)  // check if the neighbour has a smaller rank
            {
                int jrep = Find(nb, PARENT);  // representative of the tree
                int    w = FORESTROOT[jrep];        // forest root of the tree
                
                if (v != w)
                {
                    // Merge S_v and S_w=S_{jrep}
                    UnionBySize(v, jrep, PARENT, FORESTROOT, SIZEchange, SIZE, MASS);
                }
            }
        }
    }
    
    //    // flag forest roots
    //    for (int i = 0; i < m; i++)
    //    {
    //        if (PARENT[i] == i) ROOT[FORESTROOT[i]] = FORESTROOT[i];
    //    }
    
    return ClusterResults{std::move(SIZEchange), std::move(SIZE), std::move(MASS)};
}

// Compute TFCE scores for all voxels
struct TFCEResults
{
    std::vector<double> TFCE;
    double maxTFCE;
};
TFCEResults allTFCE(int m,
                    double E,
                    double H,
                    std::vector<double>& stats,   // sorted in non-ascending order
                    std::vector<int>&    sizes,
                    std::vector<int>&    ords,
                    std::vector<int>&    parents)
{
    // initialize output:
    double maxTFCE = 0.0;
    std::vector<double> TFCEscore(m, 0.0);
    
    double invH = 1.0 / (H + 1.0);
    // ------ (1) For special case (default settings)
    if (E == 0.5 && H == 2.0)
    {
        for (int i = m-1; i >= 0; i--)
        {
            int v = ords[i];
            int p = parents[v];
            
            double sv = stats[v];
            if (sv <= 0.0) continue;
            
            double sizeE  = std::sqrt((double)sizes[v]);
            double statH  = sv * sv * sv;
            double statHp = (stats[p] > 0.0) ? stats[p]*stats[p]*stats[p] : 0.0;
            
            double val;
            if (p != v)  // check if root has been reached
            {
                val = TFCEscore[p] + sizeE * (statH - statHp) * invH;
            }
            else
            {
                val = sizeE * statH * invH;
            }
            TFCEscore[v] = val;
            
            // compute max(TFCE)
            if (val > maxTFCE) maxTFCE = val;
        }
        
        return TFCEResults{std::move(TFCEscore), maxTFCE};
    }
    
    // ------ (2) For general case
    
    const double Hplus = H + 1.0;
    
    // pre-calculation
    std::vector<double> statsH(m, 0.0);
    std::vector<double> sizesE(m);
    for (int i = 0; i < m; i++)
    {
        double s = stats[i];

        // stats^H（one-sided）
        if (s > 0.0)
        {
            if (H == 2.0)
                statsH[i] = s * s * s;
            else if (H == 1.0)
                statsH[i] = s * s;
            else if (H == 3.0)
                statsH[i] = s * s * s * s;
            else
                statsH[i] = std::exp(Hplus * std::log(s));
            //statsH[i] = (H == 2.0) ? (s * s * s) : std::pow(s, Hplus);
        }

        // sizes^E
        if (E == 0.5)
            sizesE[i] = std::sqrt((double)sizes[i]);
        else
            sizesE[i] = std::exp(E * std::log((double)sizes[i]));
        //sizesE[i] = (E == 0.5) ? std::sqrt((double)sizes[i]) : std::pow((double)sizes[i], E);
    }
    
    for (int i = m-1; i >= 0; i--)
    {
        int v = ords[i];
        int p = parents[v];
        
        if (stats[v] <= 0) continue;
        
        double val;
        if (p != v)
        {
            val = TFCEscore[p] + sizesE[v] * (statsH[v] - statsH[p]) * invH;
        }
        else
        {
            val = sizesE[v] * statsH[v] * invH;
        }
        TFCEscore[v] = val;
        
        // compute max(TFCE)
        if (val > maxTFCE) maxTFCE = val;
    }
    
    return TFCEResults{std::move(TFCEscore), maxTFCE};
}

// Compute max(TFCE)
double maxTFCE(int m,
               double E,
               double H,
               std::vector<double>& stats,
               std::vector<int>&    sizes,
               std::vector<int>&    ords,
               std::vector<int>&    parents)
{
    // initialize output:
    double maxTFCE = 0.0;
    
    double invH = 1.0 / (H + 1.0);
    std::vector<double> TFCEscore(m, 0.0);
    
    // ------ (1) For special case (default settings)
    if (E == 0.5 && H == 2.0)
    {
        for (int i = m-1; i >= 0; i--)
        {
            int v = ords[i];
            int p = parents[v];
            
            double sv = stats[v];
            if (sv <= 0) continue;
            
            double sizeE  = std::sqrt((double)sizes[v]);
            double statH  = sv * sv * sv;
            double statHp = (stats[p] > 0.0) ? stats[p]*stats[p]*stats[p] : 0.0;
            
            double val;
            if (p != v)  // check if root has been reached
            {
                val = TFCEscore[p] + sizeE * (statH - statHp) * invH;
            }
            else
            {
                val = sizeE * statH * invH;
            }
            TFCEscore[v] = val;
            
            // find max(TFCE)
            if (val > maxTFCE) maxTFCE = val;
        }
        return maxTFCE;
    }
    
    // ------ (2) For general case
    
    const double Hplus = H + 1.0;
    
    // pre-calculation
    std::vector<double> statsH(m, 0.0);
    std::vector<double> sizesE(m);
    for (int i = 0; i < m; i++)
    {
        double s = stats[i];
        
        // stats^H（one-sided）
        if (H == 2.0)
            statsH[i] = s * s * s;
        else if (H == 1.0)
            statsH[i] = s * s;
        else if (H == 3.0)
            statsH[i] = s * s * s * s;
        else
            statsH[i] = std::exp(Hplus * std::log(s));
        //statsH[i] = (s > 0.0) ? std::pow(s, H + 1.0) : 0.0;
        
        // sizes^E
        if (E == 0.5)
            sizesE[i] = std::sqrt((double)sizes[i]);
        else
            sizesE[i] = std::exp(E * std::log((double)sizes[i]));
    }
    
    for (int i = m-1; i >= 0; i--)
    {
        int v = ords[i];
        int p = parents[v];
        
        if (stats[v] <= 0) continue;
        
        double val;
        if (p != v)
        {
            val = TFCEscore[p] + sizesE[v] * (statsH[v] - statsH[p]) * invH;
        }
        else
        {
            val = sizesE[v] * statsH[v] * invH;
        }
        TFCEscore[v] = val;
        
        // max(TFCE)
        if (val > maxTFCE) maxTFCE = val;
    }
    
    return maxTFCE;
}

// find the minimum (rightmost) voxel height h>=t0
int findminVox(double t0,
               int m,
               std::vector<double>& stats,
               std::vector<int>&    ords)
{
    int l = 0;
    int r = m-1;
    int left = 0;
    int right = m-1;
    
    if (stats[ords[r]] >= t0) return r;
    if (stats[ords[l]] < t0) return -1;
    
    while (l < r)
    {
        int mid = (l+r+1)/2;  // (1) binary search
        if (stats[ords[mid]] >= t0)
        {
            l = mid;
        }
        else
        {
            r = mid - 1;
        }
        
        right--;              // (2) linear search part
        if (stats[ords[right]] >= t0) return right;
        left++;
        if (stats[ords[left]] < t0) return left-1;
    }
    
    return l;
}

// Compute max cluster stats (extent & mass) for CDT=t0
std::vector<double> maxCS(double t0,
                          int m,
                          std::vector<double>& stats,
                          std::vector<int>&    sizes,
                          std::vector<double>& masses,
                          std::vector<int>&    ords,
                          std::vector<int>&    ranks,
                          std::vector<int>&    parents)
{
    int i = findminVox(t0, m, stats, ords);
    if (i < 0) return std::vector<double> {0.0, 0.0};
    
    // find maxCE & maxCM
    int maxCE = 0;
    double maxCM = 0.0;
    for (int j = 0; j <= i; j++)
    {
        int v = ords[j];
        int p = parents[v];
        
        if (p == v || ranks[p] > i+1)
        {
            int    s = sizes[v];
            double m = masses[v];
            if (s > maxCE) maxCE = s;
            if (m > maxCM) maxCM = m;
        }
    }
    
    return std::vector<double> {(double)maxCE, maxCM};
}

// Compute original cluster stats (extent & mass) for CDT=t0
struct OrgCSResults
{
    std::vector<int> orgCE;
    std::vector<double> orgCM;
};
OrgCSResults orgCS(double t0,
                   int m,
                   std::vector<double>& stats,
                   std::vector<int>&    sizes,
                   std::vector<double>& masses,
                   std::vector<int>&    ords,
                   std::vector<int>&    ranks,
                   std::vector<int>&    parents)
{
    int i = findminVox(t0, m, stats, ords);
    if (i < 0) return {};
    
    std::vector<int> orgCE;
    std::vector<double> orgCM;
    orgCE.reserve(i+1);
    orgCM.reserve(i+1);
    for (int j = 0; j <= i; j++)
    {
        int v = ords[j];
        int p = parents[v];
        
        if (p == v || ranks[p] > i+1)
        {
            orgCE.push_back(sizes[v]);
            orgCM.push_back(masses[v]);
        }
    }
    
    return OrgCSResults{std::move(orgCE), std::move(orgCM)};
}


//------------------------- (3) TFCE/CS/CM implementations -------------------------//

struct GraphPre
{
    std::vector<int> idt;
    std::vector<std::vector<int>> adj;
};
GraphPre precompute_graph(Rcpp::IntegerVector& mask,
                          Rcpp::IntegerVector& dims,
                          int conn = 26)
{
    int vox = mask.size();
    
    std::vector<int> idt;
    for (int i = 0; i < vox; i++)
    {
        if (mask[i] != 0)
        {
            idt.push_back(i);
            mask[i] = idt.size();  // update mask: 0 -> 0, 1 -> 1,2,...
        }
    }
    int m = idt.size();
    
    std::vector<std::vector<int>> adj = findAdjList(mask, dims, idt, m, conn);
    
    return GraphPre{std::move(idt), std::move(adj)};
}

// eTFCE computation (for permuted data)
struct ETFCEPermResults
{
    double mTFCE;
    double mCZ;
    double mCM;

    double time_sort;
    double time_tree;
    double time_tfce;
    double time_clus;
};
ETFCEPermResults eTFCE(std::vector<double>&           tval,
                       std::vector<std::vector<int>>& adj,
                       std::vector<int>&              idt,
                       int m,
                       double E,
                       double H,
                       double t0)
{
    auto t_start = std::chrono::steady_clock::now();
    
    // ---- sort ----
    std::vector<int> ord(m);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return tval[a] > tval[b]; });
    
    // ---- rank ----
    std::vector<int> rank(m);
    for (int i = 0; i < m; i++)
    {
        rank[ord[i]] = i + 1;
    }
    
    auto t_tree0 = std::chrono::steady_clock::now();
    
    // ---- find all clusters ----
    ClusterResults cl = findClusters(m, adj, tval, ord, rank);
    
    auto t_tree1 = std::chrono::steady_clock::now();
    
    // ---- TFCE ----
    double mTFCE = maxTFCE(m, E, H, tval, cl.SIZE, ord, cl.SIZEchange);
    
    auto t_tfce1 = std::chrono::steady_clock::now();
    
    // ---- cluster size & mass ----
    std::vector<double> cs = maxCS(t0, m, tval, cl.SIZE, cl.MASS, ord, rank, cl.SIZEchange);
    
    auto t_clus1 = std::chrono::steady_clock::now();
    
    return ETFCEPermResults{mTFCE, cs[0], cs[1],
        std::chrono::duration<double>(t_tree0 - t_start).count(),   // sorting
        std::chrono::duration<double>(t_tree1 - t_tree0).count(),   // tree construction
        std::chrono::duration<double>(t_tfce1 - t_tree1).count(),   // TFCE
        std::chrono::duration<double>(t_clus1 - t_tfce1).count()};  // cluster size & mass
}

// eTFCE computation (for original data)
struct ETFCEResults
{
    double maxTFCE;
    std::vector<double> TFCE;
    std::vector<int> CE;
    std::vector<double> CM;
    
    double time_sort;
    double time_tree;
    double time_tfce;
    double time_clus;
};
ETFCEResults eTFCE_org(std::vector<double>&           tval,
                       std::vector<std::vector<int>>& adj,
                       std::vector<int>&              idt,
                       int m,
                       double E,
                       double H,
                       double t0)
{
    auto t_start = std::chrono::steady_clock::now();
    
    // ---- sort ----
    std::vector<int> ord(m);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return tval[a] > tval[b]; });
    
    // ---- rank ----
    std::vector<int> rank(m);
    for (int i = 0; i < m; i++)
    {
        rank[ord[i]] = i + 1;
    }
    
    auto t_tree0 = std::chrono::steady_clock::now();
    
    // ---- find all clusters ----
    ClusterResults cl = findClusters(m, adj, tval, ord, rank);
    
    auto t_tree1 = std::chrono::steady_clock::now();
    
    // ---- TFCE ----
    TFCEResults res = allTFCE(m, E, H, tval, cl.SIZE, ord, cl.SIZEchange);
    
    auto t_tfce1 = std::chrono::steady_clock::now();
    
    // ---- cluster size & mass ----
    OrgCSResults out = orgCS(t0, m, tval, cl.SIZE, cl.MASS, ord, rank, cl.SIZEchange);
    
    auto t_clus1 = std::chrono::steady_clock::now();
    
    return ETFCEResults{res.maxTFCE,std::move(res.TFCE), std::move(out.orgCE), std::move(out.orgCM),
        std::chrono::duration<double>(t_tree0 - t_start).count(),   // sorting
        std::chrono::duration<double>(t_tree1 - t_tree0).count(),   // tree construction
        std::chrono::duration<double>(t_tfce1 - t_tree1).count(),   // TFCE
        std::chrono::duration<double>(t_clus1 - t_tfce1).count()};  // cluster size & mass
}

// [[Rcpp::export]]
Rcpp::List eTFCE_perms(Rcpp::NumericVector& conts,  // dim1*dim2*dim3 * n
                       Rcpp::IntegerVector& mask,
                       Rcpp::IntegerVector& dims,
                       int n,
                       int n_perm = 5000,
                       double E = 0.5,
                       double H = 2.0,
                       double t0 = 3.1,
                       int conn = 26,
                       int seed = 12345,  // optional seed specification
                       Rcpp::Nullable<Rcpp::NumericVector> sign_matrix = R_NilValue,  // optional permutation matrix: n_perm * n
                       bool auto_signflip = true)
{
    auto t_start = std::chrono::steady_clock::now();
    
    int vox = dims[0] * dims[1] * dims[2];
    
    // ---- graph (ONCE) ----
    GraphPre pre = precompute_graph(mask, dims, conn);
    std::vector<std::vector<int>>& adj = pre.adj;
    std::vector<int>& idt = pre.idt;
    int m = idt.size();
    
    auto t1 = std::chrono::steady_clock::now();
    
    // ---- outputs ----
    Rcpp::NumericVector maxTFCE(n_perm);
    Rcpp::NumericVector maxCZ(n_perm);
    Rcpp::NumericVector maxCM(n_perm);
    Rcpp::NumericVector time_sec(11);
    
    // ---- constants ----
    double sqrt_n = std::sqrt((double)n);
    double inv_n  = 1.0 / n;
    double inv_n1 = 1.0 / (n - 1);
    
    // =========================
    // (1) ORIGINAL RESULT
    // =========================
    std::vector<double> tval0(m, 0.0);
    std::vector<double> SS(m, 0.0);
    
    for (int i = 0; i < m; i++)
    {
        int v = idt[i];
        
        double sum = 0.0;
        double  ss = 0.0;
        
        for (int j = 0; j < n; j++)
        {
            double val = conts[v + vox*j];
            sum += val;
            ss  += val * val;
        }
        
        SS[i] = ss;
        
        double mean = sum * inv_n;
        double var  = (ss - n * mean * mean) * inv_n1;
        
        if (var > 0) tval0[i] = sqrt_n * mean / std::sqrt(var);
    }
    
    // observed TFCE
    ETFCEResults org = eTFCE_org(tval0, adj, idt, m, E, H, t0);
    //std::vector<double> orgTFCE(vox, 0.0);
    //for (int i = 0; i < m; i++)
    //{
    //    orgTFCE[idt[i]] = org.TFCE[i];
    //}
    
    auto t2 = std::chrono::steady_clock::now();
    
    // =========================
    // (2) PERMUTATIONS
    // =========================
    std::vector<double> signs(n);
    std::vector<double> tval(m);
    
    double time_sort = 0.0;
    double time_tree = 0.0;
    double time_tfce = 0.0;
    double time_clus = 0.0;
    
    // ---- prepare sign_matrix (ONCE) ----
    Rcpp::NumericVector sm;
    bool use_sign_matrix = false;
    
    if (!auto_signflip && sign_matrix.isNotNull())
    {
        sm = Rcpp::NumericVector(sign_matrix);
        use_sign_matrix = true;
    }
    
    // ---- RNG ----
    std::mt19937 rng(seed);
    std::bernoulli_distribution dist(0.5);
    
    for (int b = 0; b < n_perm; b++)
    {
        // ---- signs ----
        if (use_sign_matrix)
        {
            for (int j = 0; j < n; j++)
            {
                signs[j] = sm[b*n + j];
            }
        }
        else
        {
            for (int j = 0; j < n; j++)
            {
                signs[j] = dist(rng) ? 1.0 : -1.0;
            }
        }
        
        // ---- t-values ----
        std::fill(tval.begin(), tval.end(), 0.0);
        
        for (int i = 0; i < m; i++)
        {
            int v = idt[i];
            
            double sum = 0.0;
            
            for (int j = 0; j < n; j++)
            {
                sum += conts[v + vox*j] * signs[j];
            }
            
            double mean = sum * inv_n;
            double var  = (SS[i] - n * mean * mean) * inv_n1;
            
            if (var > 0) tval[i] = sqrt_n * mean / std::sqrt(var);
        }
        
        // ---- TFCE ----
        ETFCEPermResults res = eTFCE(tval, adj, idt, m, E, H, t0);
        
        maxTFCE[b] = res.mTFCE;
        maxCZ[b]   = res.mCZ;
        maxCM[b]   = res.mCM;
        
        time_sort += res.time_sort;
        time_tree += res.time_tree;
        time_tfce += res.time_tfce;
        time_clus += res.time_clus;
    }
    
    auto t3 = std::chrono::steady_clock::now();
    
    time_sec[0] = std::chrono::duration<double>(t1 - t_start).count();  // graph
    time_sec[1] = std::chrono::duration<double>(t2 - t1).count();       // original
    time_sec[2] = std::chrono::duration<double>(t3 - t2).count();       // permutation
    time_sec[3] = org.time_sort;
    time_sec[4] = org.time_tree;
    time_sec[5] = org.time_tfce;
    time_sec[6] = org.time_clus;
    time_sec[7] = time_sort / n_perm;
    time_sec[8] = time_tree / n_perm;
    time_sec[9] = time_tfce / n_perm;
    time_sec[10] = time_clus / n_perm;
    
    
    return Rcpp::List::create(Rcpp::Named("n_perm") = n_perm,
                              Rcpp::Named("time_sec") = time_sec,
                              Rcpp::Named("maxTFCE") = maxTFCE,
                              Rcpp::Named("maxCZ") = maxCZ,
                              Rcpp::Named("maxCM") = maxCM,
                              Rcpp::Named("mTFCE") = org.maxTFCE,
                              Rcpp::Named("orgTFCE") = Rcpp::wrap(org.TFCE),
                              Rcpp::Named("orgCZ") = Rcpp::wrap(org.CE),
                              Rcpp::Named("orgCM") = Rcpp::wrap(org.CM));
}
