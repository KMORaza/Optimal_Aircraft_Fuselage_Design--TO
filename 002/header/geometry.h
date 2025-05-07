#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

class Geometry {
public:
    Geometry(double radius, double length, double skin_thickness,
             int num_frames, int num_stringers);

    int getNumNodes() const { return num_nodes_; }
    int getNumElements() const { return num_elements_; }
    double getRadius() const { return radius_; }
    double getLength() const { return length_; }
    double getSkinThickness() const { return skin_thickness_; }
    int getNumFrames() const { return num_frames_; }
    int getNumStringers() const { return num_stringers_; }
    const std::vector<std::vector<double>>& getNodes() const { return nodes_; }
    const std::vector<std::vector<int>>& getElements() const { return elements_; }
    double getPanelArea() const;
    void updateSkinThickness(double thickness) { skin_thickness_ = thickness; }

private:
    void generateMesh();

    double radius_;
    double length_;
    double skin_thickness_;
    int num_frames_;
    int num_stringers_;
    int num_nodes_;
    int num_elements_;
    std::vector<std::vector<double>> nodes_; // {x, y, z}
    std::vector<std::vector<int>> elements_; // Node indices per element
};

#endif // GEOMETRY_H