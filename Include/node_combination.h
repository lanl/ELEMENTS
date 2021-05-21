#ifndef NODE_COMBINATION_H
#define NODE_COMBINATION_H  

#include "matar.h"
#include "utilities.h"

using namespace utils;
class Node_Combination;
bool operator< (const Node_Combination &object1, const Node_Combination &object2);

class Node_Combination {

public:
    
  CArray<size_t> node_set;
  size_t patch_gid;

  //Default Constructor
  Node_Combination(){}

  //Constructor with initialization
  Node_Combination(CArray<size_t> &nodes_init) {
    node_set = nodes_init;
  }

  //Destructor
  ~Node_Combination( ) {}

  //overload = operator
  Node_Combination& operator= (Node_Combination &not_this){
    node_set = not_this.node_set;
    return *this;
  }

  //overload = operator
  bool operator== (Node_Combination &not_this){
    int this_size = this->node_set.size();
    //check if this node combination is identical
    //first check size of the combination
    if(this_size!=not_this.node_set.size())
      return false;

    //check if the nodes in the set are the same; sort them to simplify
    std::sort(this->node_set.get_pointer(),this->node_set.get_pointer()+this->node_set.size());
    std::sort(not_this.node_set.get_pointer(),not_this.node_set.get_pointer()+not_this.node_set.size());\

    //loop through the sorted nodes to check for equivalence
    for(int i = 0; i < this_size; i++)
      if(this->node_set(i)!=not_this.node_set(i)) return false;

    return true;
    
  }

  
    

private:
    int num_dim_;

};

#endif // end STATE_H
