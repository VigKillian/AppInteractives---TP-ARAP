#ifndef LAPLACIANWEIGHTS_H
#define LAPLACIANWEIGHTS_H

#include <vector>
#include <map>
#include "Mesh.h"



//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//
//
// This class describes the implementation of the well known cotangent weights
// Here, weight( v1 , v2 ) is half the sum of the cotangent angles of the opposite corners in the triangles (v1,v2,other)
//   Careful ! these weights can be negative !
//
// Additionally, we provide a scheme for non-symmetric (!!!) POSITIVE edge weights
//   Careful ! these weights are non-symetric ! e_ij != e_ji
//
//-------------------------------------------------------------------------------------//
//-------------------------------------------------------------------------------------//





// access to an edge is of complexity O( log(val) ), with val the average valence of the vertices

//---------------------------------   YOU DO NOT NEED TO CHANGE THE FOLLOWING CODE  --------------------------------//
class LaplacianWeights{
private:
    unsigned int n_vertices;
    std::vector< std::map< unsigned int , double > > edge_weights;
    std::vector< double > vertex_weights;

public:
    LaplacianWeights() : n_vertices(0) {}
    void clear()
    {
        n_vertices = 0;
        edge_weights.clear();
        vertex_weights.clear();
    }
    ~LaplacianWeights() { clear(); }

    double sumVertexWeights() const {
        double s = 0.0;
        for(unsigned int v = 0 ; v < n_vertices ; ++v) {
            s += vertex_weights[v];
        }
        return s;
    }

    void resize( unsigned int nVertices )
    {
        clear();
        if( nVertices > 0 )
        {
            n_vertices = nVertices;
            edge_weights.resize(nVertices);
            vertex_weights.resize(nVertices , 0.0);

            for( unsigned int v = 0 ; v < nVertices ; ++v )
            {
                edge_weights[v].clear();
                vertex_weights[v] = 0.0;
            }
        }
    }
    unsigned int get_n_adjacent_edges( unsigned int vertex_index ) const
    {
        return edge_weights[vertex_index].size();
    }
    double get_edge_weight( unsigned int v1 , unsigned int v2 ) const
    {
        std::map< unsigned int , double >::const_iterator it = edge_weights[v1].find(v2);
        if( it == edge_weights[v1].end() ) return 0.0;
        return it->second;
    }
    unsigned int get_n_vertices() const
    {
        return n_vertices;
    }


    std::map< unsigned int , double >::iterator get_weight_of_adjacent_edges_it_begin( unsigned int v1 )
    {
        return edge_weights[v1].begin();
    }
    std::map< unsigned int , double >::iterator get_weight_of_adjacent_edges_it_end( unsigned int v1 )
    {
        return edge_weights[v1].end();
    }

    std::map< unsigned int , double >::const_iterator get_weight_of_adjacent_edges_it_begin( unsigned int v1 ) const
    {
        return edge_weights[v1].begin();
    }
    std::map< unsigned int , double >::const_iterator get_weight_of_adjacent_edges_it_end( unsigned int v1 ) const
    {
        return edge_weights[v1].end();
    }

    double get_vertex_weight( unsigned int v ) const
    {
        return vertex_weights[v];
    }

    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //---------------------------------  CODE TO CHANGE  --------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//

    //--------------------------   Cotangent weights:   ------------------------//
    // Weight of edge eij : i<->j is the sum of cotangent of opposite angles divided by 2
    // wij = 1/2 * (cot(alpha_ij) + cot(beta_ij)) alpha_ij and beta_ij being the two opposite angles of the edge ij
    
    void buildCotangentWeightsOfTriangleMesh( const Mesh& mesh)
    {
        resize(mesh.V.size());

        // For each triangle : 
            // Compute its edge p0, p1, p2
            // Compute opposite angle for each edge
            // Compute cotangent of this angle
            // Add to edge weight 

        for(unsigned int t=0; t<mesh.T.size(); t++){
            unsigned int v0 = mesh.T[t][0];
            unsigned int v1 = mesh.T[t][1];
            unsigned int v2 = mesh.T[t][2];

            Eigen::Vector3d p0(mesh.V[v0].pInit[0], mesh.V[v0].pInit[1], mesh.V[v0].pInit[2]);
            Eigen::Vector3d p1(mesh.V[v1].pInit[0], mesh.V[v1].pInit[1], mesh.V[v1].pInit[2]);
            Eigen::Vector3d p2(mesh.V[v2].pInit[0], mesh.V[v2].pInit[1], mesh.V[v2].pInit[2]);

            Eigen::Vector3d u = p1 - p0;
            Eigen::Vector3d v = p2 - p0;
            double cot0 = u.dot(v) / u.cross(v).norm();
            edge_weights[v1][v2] += 0.5 * cot0;
            edge_weights[v2][v1] += 0.5 * cot0;

            u = p0 - p1;
            v = p2 - p1;
            double cot1 = u.dot(v) / u.cross(v).norm();
            edge_weights[v0][v2] += 0.5 * cot1;
            edge_weights[v2][v0] += 0.5 * cot1;

            u = p0 - p2;
            v = p1 - p2;
            double cot2 = u.dot(v) / u.cross(v).norm();
            edge_weights[v0][v1] += 0.5 * cot2;
            edge_weights[v1][v0] += 0.5 * cot2;
        }
    }

    //---------------------------------   YOU DO NOT NEED TO CHANGE THE FOLLOWING CODE  --------------------------------//

    //------------------------------------   Barycentric weights:  --------------------------------//
    // Weight of vertex i is its barycentric area (sum of areas of connected triangles divided by 3)
    // Weight of directed edge i->j is its barycentric area divided by the weight of vertex i
    // CAREFUL: weight(i->j) is different from weight(j->i) !
    //          EDGE WEIGHTS ARE NOT SYMMETRIC !
    // These weights are all positive

    template< class vertex_t , class triangle_t >
    void buildBarycentricWeightsOfTriangleMesh(
    const std::vector< vertex_t > & vertices,
    const std::vector< triangle_t > & triangles)
{
    // Initialise les structures (resize doit initialiser n_vertices, edge_weights, vertex_weights)
    resize(vertices.size());

    const double eps = 1e-6;
    const double cotan_eps = 1e-6;
    const double cotan_max = std::cos(cotan_eps) / std::sin(cotan_eps);

    auto safe_cotan = [&](const Vec3 &a, const Vec3 &b) -> double {
        // cot(angle between a and b) = (a·b) / |a×b|
        Vec3 cross = Vec3::cross(a, b);
        double denom = cross.norm();
        if (denom <= eps) return 0.0; // angle degenerate => cot ~ 0 (robuste)
        double val = a.dot(b) / denom;
        if (std::isnan(val) || std::isinf(val)) return 0.0;
        // clamp to avoid huge values when angle ~ 0 or ~ pi
        if (val > cotan_max) val = cotan_max;
        if (val < -cotan_max) val = -cotan_max;
        return val;
    };

    // Pour chaque triangle : calculer les cotangentes aux 3 sommets
    for (unsigned int t = 0; t < triangles.size(); ++t)
    {
        unsigned int i0 = triangles[t][0];
        unsigned int i1 = triangles[t][1];
        unsigned int i2 = triangles[t][2];

        Vec3 p0(vertices[i0][0], vertices[i0][1], vertices[i0][2]);
        Vec3 p1(vertices[i1][0], vertices[i1][1], vertices[i1][2]);
        Vec3 p2(vertices[i2][0], vertices[i2][1], vertices[i2][2]);
        

        // vecteurs autour de chaque sommet
        // angle at p0 between (p1-p0) and (p2-p0)
        Vec3 v01 = p1 - p0;
        Vec3 v02 = p2 - p0;
        double cot0 = safe_cotan(v01, v02);

        // angle at p1 between (p2-p1) and (p0-p1)
        Vec3 v12 = p2 - p1;
        Vec3 v10 = p0 - p1;
        double cot1 = safe_cotan(v12, v10);

        // angle at p2 between (p0-p2) and (p1-p2)
        Vec3 v20 = p0 - p2;
        Vec3 v21 = p1 - p2;
        double cot2 = safe_cotan(v20, v21);

        // ajouter chaque cotangente à l'arête opposée
        // cot0 (angle at i0) -> contributes to edge (i1, i2)
        edge_weights[i1][i2] += cot0;
        edge_weights[i2][i1] += cot0;

        // cot1 -> edge (i2, i0)
        edge_weights[i2][i0] += cot1;
        edge_weights[i0][i2] += cot1;

        // cot2 -> edge (i0, i1)
        edge_weights[i0][i1] += cot2;
        edge_weights[i1][i0] += cot2;
    }

    // Maintenant accumuler la somme des poids incident par sommet (utile pour la diagonale de la Laplacienne)
    for (unsigned int v = 0; v < n_vertices; ++v)
    {
        double sum = 0.0;
        for (std::map<unsigned int, double>::const_iterator it = edge_weights[v].begin();
             it != edge_weights[v].end(); ++it)
        {
            sum += it->second;
        }
        vertex_weights[v] = sum;
    }

    // --- Optionnel : normaliser les poids autour de chaque sommet pour qu'ils somment à 1
    // (Décommente si tu veux reproduire exactement le comportement "barycentric" précédent)
    /*
    for (unsigned int v = 0; v < n_vertices; ++v)
    {
        double vsum = vertex_weights[v];
        if (vsum <= 0.0) continue;
        for (std::map<unsigned int,double>::iterator it = edge_weights[v].begin();
             it != edge_weights[v].end(); ++it)
        {
            edge_weights[v][it->first] = it->second / vsum;
        }
        // après normalisation la 'vertex_weights[v]' vaut 1.0 si tu veux.
        vertex_weights[v] = 1.0;
    }
    */
}


};





#endif // LAPLACIANWEIGHTS_H
