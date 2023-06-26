#include <gmsh.h>
#include <vector>
#include <iostream>

using namespace std;
struct Node{
    size_t id;
    double x;
    double y;
    double z;
};

struct Element{
    size_t id;
    size_t type;
    vector<int> nodeid;
};
void readMSHFile(const string& fileName, vector<Node>& nodes, vector<Element>& elements) {
    gmsh::initialize();
    gmsh::open(fileName);

    vector<size_t> nodeTags;
    vector<double> coord, parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord);

    for(size_t i = 0; i < nodeTags.size(); i++) {
        Node node;
        node.id = nodeTags[i];
        node.x = coord[3*i];
        node.y = coord[3*i+1];
        node.z = coord[3*i+2];
        nodes.push_back(node);
    }

    vector<int> elementTypes;
    gmsh::model::mesh::getElementTypes(elementTypes);

    for(auto type : elementTypes) {
        vector<size_t> elementTags, nodeTags;
        gmsh::model::mesh::getElementsByType(type, elementTags, nodeTags);
        size_t numNodes = nodeTags.size() / elementTags.size();
        for(size_t i = 0; i < elementTags.size(); i++) {
            Element element;
            element.id = elementTags[i];
            element.type = type;
            for(size_t j = 0; j < numNodes; j++) {
                element.nodeid.push_back(nodeTags[i*numNodes+j]);
            }
            elements.push_back(element);
        }
    }

    gmsh::finalize();
}


int main() {
	cout<<"check"<<endl;
    vector<Node> nodes;
    vector<Element> elements;
	cout<<"check"<<endl;
    readMSHFile("./msh/mesh.msh", nodes, elements);
    // Now nodes and elements are filled with data from the .msh file
	cout<<nodes.size()<<endl;
	for(auto node:nodes){
		cout<<node.x<<endl;
	}
    return 0;
}
