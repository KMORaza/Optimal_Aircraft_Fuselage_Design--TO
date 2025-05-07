#include "geometry.h"
#include <cmath>

Geometry::Geometry(double radius, double length, double skin_thickness,
                   int num_frames, int num_stringers)
    : radius_(radius), length_(length), skin_thickness_(skin_thickness),
      num_frames_(num_frames), num_stringers_(num_stringers),
      num_nodes_(0), num_elements_(0) {
    generateMesh();
}

void Geometry::generateMesh() {
    // Generate a cylindrical mesh
    int nodes_per_ring = num_stringers_;
    int num_rings = num_frames_ + 1;
    num_nodes_ = nodes_per_ring * num_rings;

    // Generate nodes
    nodes_.resize(num_nodes_, std::vector<double>(3));
    double dz = length_ / num_frames_;
    for (int i = 0; i < num_rings; ++i) {
        double z = i * dz;
        for (int j = 0; j < nodes_per_ring; ++j) {
            double theta = 2 * M_PI * j / nodes_per_ring;
            int idx = i * nodes_per_ring + j;
            nodes_[idx][0] = radius_ * std::cos(theta);
            nodes_[idx][1] = radius_ * std::sin(theta);
            nodes_[idx][2] = z;
        }
    }

    // Generate quad elements
    num_elements_ = num_frames_ * num_stringers_;
    elements_.resize(num_elements_, std::vector<int>(4));
    for (int i = 0; i < num_frames_; ++i) {
        for (int j = 0; j < num_stringers_; ++j) {
            int idx = i * num_stringers_ + j;
            int n0 = i * nodes_per_ring + j;
            int n1 = i * nodes_per_ring + (j + 1) % nodes_per_ring;
            int n2 = (i + 1) * nodes_per_ring + (j + 1) % nodes_per_ring;
            int n3 = (i + 1) * nodes_per_ring + j;
            elements_[idx] = {n0, n1, n2, n3};
        }
    }
}

double Geometry::getPanelArea() const {
    // Approximate area of a single panel
    return M_PI * radius_ * length_ / num_elements_;
}