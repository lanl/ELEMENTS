#ifndef NDEBUG
#define NDEBUG 0
#endif 

template <typename T>
class DynamicRaggedRightList {
    
private:
    
    size_t *start_index_;
    T * array_;
    
    size_t dim1_;
    size_t dim2_;
    size_t length_;

    size_t index_;
    size_t irow_;
    
public:
    
    // Default constructor
    DynamicRaggedRightList();
    
    // Overload constructor
    DynamicRaggedRightList(size_t dim1, size_t dim2);
    
    // A method to reset the start indices
    void reset();

    // A method to append a data to the array
    void append(T id);

    // A method to set the start_index for new row
    void rowmark();

    // A method to return stride size for row i
    size_t stride(size_t i) const;

    // A method to sort the data for each row
    void sort();

    // Overload operator() to access data as array(i,j),
    // where i=[0:N-1], j=[0:stride(i)]
    T& operator()(size_t i, size_t j) const;

    // Overload operator= to perform copy assignment
    DynamicRaggedRightList& operator=(const DynamicRaggedRightList temp);

    // Destructor
    ~DynamicRaggedRightList();
}; // end of DynamicRaggedRightList

template <typename T>
DynamicRaggedRightList<T>::DynamicRaggedRightList() {}

// Overloaded constructor
template <typename T>
DynamicRaggedRightList<T>::DynamicRaggedRightList(size_t dim1, size_t dim2) {
    // Dimensions of the array
    // Number of rows
    dim1_ = dim1;
    // Max buffer size, i.e., max number of elements in any row
    dim2_ = dim2;
    // Max number of elements in array
    length_ = (dim1_ * dim2_);

    // Create memory on heap for values
    array_ = new T[length_];

    // Create array for holding start indices
    start_index_ = new size_t[(dim1_ + 1)];     // Note the dim1_ + 1

    // Initialize the indices
    index_ = 0;
    irow_  = 0;

    // 1D array starts at 0
    for (size_t i = 0; i < (dim1_ + 1); i++) {
        start_index_[i] = 0;    // Initialize each row's stride to be 0
    }
}   // End constructor

// Method to reset the start indices
template <typename T>
// void DynamicRaggedRightList<T>::reset() const {
void DynamicRaggedRightList<T>::reset() {
    // Reset the indices
    index_ = 0;
    irow_  = 0;

    // 1D array starts at 0
    for (size_t i = 0; i < (dim1_ + 1); i++) {
        start_index_[i] = 0;  // Reset each row's stride to be 0
    }
}   // End reset()

// A method to append an element to the array
template <typename T>
void DynamicRaggedRightList<T>::append(T id) {
    // Make index start at next element
    index_++;

    // Expand the array, if necessary
    if (index_ >= length_) {
        T* tmp = new T[(2 * length_)];
        
        // NOTE: May not work on GPU; look into this
        std::copy(array_, (array_ + length_), tmp);

        // Deallocate old memory
        delete[] array_;

        // Make array_ point to new memory and update length
        array_ = tmp;
        length_ = (2 * length_);
    }

    array_[(index_ - 1)] = id;
}   // End append()

// A method to set the start index for each new row
template <typename T>
void DynamicRaggedRightList<T>::rowmark() {
    irow_++;

    // Expand the number of rows, if necessary
    if (irow_ > dim1_) {
    // if (irow_ >= dim1_) {
        size_t* tmp = new size_t[((2 * dim1_) + 1)];
        // NOTE: May not work on GPU; look into this
        std::copy(start_index_, (start_index_ + (dim1_ + 1)), tmp);
        delete[] start_index_;

        start_index_ = tmp;
        dim1_ = (2 * dim1_);
    }
    start_index_[irow_] = index_;
}   // End rowmark()

// A method to return row i's stride
template <typename T>
size_t DynamicRaggedRightList<T>::stride(size_t i) const {
    // Check that i is a valid row index (die if i >= dim1_)
    assert(i < dim1_ && "i is out of dim1_ bounds in DynamicRaggedRightList");

    return (start_index_[(i + 1)] - start_index_[i]);
}   // End stride()

// A method to sort the data for each row
template <typename T>
void DynamicRaggedRightList<T>::sort() {
    for(size_t i = 0; i < irow_; i++) { 
        size_t ibeg = start_index_[i];
        // size_t iend = start_index_[(i + 1)] - 1;
        size_t iend = start_index_[(i + 1)];
        std::sort((array_ + ibeg), (array_ + iend));
    }
}   // End sort()

// Overload operator() to access data as array(i, j),
// where i = [0:dim1_], j = [0:stride(i)]
template <typename T>
inline T& DynamicRaggedRightList<T>::operator()(size_t i, size_t j) const {
    // Check that i is a valid row index (die if i >= dim1_)
    assert(i < dim1_ && "i is out of dim1_ bounds in DynamicRaggedRightList");
    // Check that j is a valid column index, i.e., within row i's stride
    // (die if j >= stride(i))
    assert(j < stride(i) && "j is out of stride bounds in DynamicRaggedRightList");

    return array_[start_index_[i] + j];
}   // End operator()

// Overload operator= to perform copy assignment
template <typename T>
inline DynamicRaggedRightList<T>& DynamicRaggedRightList<T>::operator=(const DynamicRaggedRightList temp) {
    // Do nothing if assignment is of the form x = x
    if (this != &temp) {
        dim1_   = temp.dim1_;
        dim2_   = temp.dim2_;
        length_ = temp.length_;
        array_  = new T[length_];
        start_index_ = new size_t[(dim1_ + 1)];     // Note the dim1_ + 1
        index_  = temp.index_;
        irow_   = temp.irow_;
        
        // Copy over row strides
        for (size_t i = 0; i < (dim1_ + 1); i++) {
            start_index_[i] = temp.start_index_[i]; 
        }
    }

    return *this;
} 

// Destructor
template <typename T>
DynamicRaggedRightList<T>::~DynamicRaggedRightList() {
    delete[] array_;
    delete[] start_index_;
}   // End destructor
