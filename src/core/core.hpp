#ifndef CORE_HPP
#define CORE_HPP

#include <cstring>

namespace mp {

    class container {
        double* data;
        unsigned long rows;
        unsigned long cols;
        void copy_data(const double* const data, const unsigned long& rows, const unsigned long& cols)
        {
            if (this->data == data) return;  // prevent self-copy
            double* old_data = this->data;
            if (data)
            {
                this->data = new double[rows * cols];
                std::memcpy(this->data, data, rows * cols * sizeof(double)); 
            }
            else this->data = nullptr;
            if (old_data) delete[] old_data;
            this->rows = rows;
            this->cols = cols;
        }
        void move_data(double* const data, const unsigned long& rows, const unsigned long& cols)
        {
            if (this->data == data) return; // prevent self-move
            double* old_data = this->data;
            this->data = data;
            if (old_data) delete[] old_data;
            this->rows = rows;
            this->cols = cols;
        }
    public:
        container() { this->data = nullptr; this->rows = 0; this->cols = 0; }
        container(const container& c) { this->copy_data(c.data, c.rows, c.cols); }
        container(container&& c) { this->move_data(c.data, c.rows, c.cols); }
        
        virtual ~container() { if (data) delete[] data; }
    };
    
    typedef container matrix;
    
    class vector : public container {
    public:
        vector() : container() {}
    };

}

#endif
