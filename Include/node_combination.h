#ifndef NODE_COMBINATION_H
#define NODE_COMBINATION_H  

#include "matar.h"

using namespace utils;

class Node_Combination {

public:
    
  CArray<size_t> node_set;
  //Constructor
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

  //overload < operator
  bool operator< (Node_Combination &not_this){
    int this_size = this->node_set.size();
    //check if this node combination is identical
    //first check size of the combination; if smaller evaluate to true
    //the set using this is then ordered first according to size of the combinations
    if(this_size<not_this.node_set.size())
      return true;
    
    //This part sorts for segments of the set where combinations have the same size
    //define < using the sort of both combinations. If the first nonequal element of the lhs combination, w.r.t to 
    //the corresponding element of the rhs, is less than the respective element of the rhs < evaluates to true
    std::sort(this->node_set.get_pointer(),this->node_set.get_pointer()+this->node_set.size());
    std::sort(not_this.node_set.get_pointer(),not_this.node_set.get_pointer()+not_this.node_set.size());\

    //loop through the sorted nodes to check for <
    for(int i = 0; i < this_size; i++){
      if(this->node_set(i)<not_this.node_set(i)) return true;
      else if(this->node_set(i)==not_this.node_set(i)) continue;
      else break;
    }

    return false;
    
  }
    

private:
    int num_dim_;

};


#endif // end STATE_H
