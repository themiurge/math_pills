#ifndef CORE_HPP
#define CORE_HPP

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#define IDX(i, j, cols) ((i)*(cols)+(j))

namespace mp {

    class container {
    protected:
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
        void move_data(double*& data, const unsigned long& rows, const unsigned long& cols)
        {
            if (this->data == data) return; // prevent self-move
            double* old_data = this->data;
            this->data = data;
            if (old_data) delete[] old_data;
            data = nullptr;
            this->rows = rows;
            this->cols = cols;
        }
    public:
        container() { this->data = nullptr; this->rows = 0; this->cols = 0; }
        container(const unsigned long& rows, const unsigned long& cols, const double& d) { this->data = ((rows*cols) ? new double[rows*cols] : nullptr); this->rows = rows; this->cols = cols; fill(d); }
        container(const container& c) { this->copy_data(c.data, c.rows, c.cols); }
        container(container&& c) { this->move_data(c.data, c.rows, c.cols); }

        inline void fill(const double& d) { for (unsigned long i = 0; i != rows*cols; ++i) data[i] = d; }

        inline const unsigned long& getRows() const { return rows; }
        inline const unsigned long& getCols() const { return cols; }

        inline bool operator==(const container& c) { if (rows != c.rows || cols != c.cols) return false; for (unsigned long i = 0; i != rows*cols; ++i) if (data[i] != c.data[i]) return false; return true; }
        inline bool operator!=(const container& c) { return !(*this == c); }

        friend void check_dimensions(const container&, const container&);
        
        virtual ~container() { if (data) delete[] data; }
    };

    void check_dimensions(const container&, const container&);
    
    class vector : public container {
    private:
        inline void check_index(const unsigned long& index) const { if (index >= rows) throw std::out_of_range("index"); }
    public:
        vector() : container(0, 1, 0.0) { }
        vector(const unsigned long& dim) : container(dim, 1, 0.0) { }
        vector(const unsigned long& dim, const double& fill_value) : container(dim, 1, fill_value) { }
        vector(const vector& v) { this->copy_data(v.data, v.rows, v.cols); }
        vector(vector&& v) { this->move_data(v.data, v.rows, v.cols); }

        const double& operator[](const unsigned long& index) const { check_index(index); return data[index]; }
        double& operator[](const unsigned long& index) { check_index(index); return data[index]; }

        friend vector operator+(const vector& v, const double& d);
        friend vector operator+(const double& d, const vector& v);
        friend vector operator-(const vector& v, const double& d);
        friend vector operator-(const double& d, const vector& v);
        friend vector operator*(const vector& v, const double& d);
        friend vector operator*(const double& d, const vector& v);
        friend vector operator/(const vector& v, const double& d);
        friend vector operator/(const double& d, const vector& v);

        inline vector& operator=(const vector& v) { this->copy_data(v.data, v.rows, v.cols); return *this; }
        inline vector& operator=(vector&& v) { this->move_data(v.data, v.rows, v.cols); return *this; }
        inline vector& operator+=(const double& d) { for (unsigned long i = 0; i != rows; ++i) data[i] += d; return *this; }
        inline vector& operator-=(const double& d) { for (unsigned long i = 0; i != rows; ++i) data[i] -= d; return *this; }
        inline vector& operator*=(const double& d) { for (unsigned long i = 0; i != rows; ++i) data[i] *= d; return *this; }
        inline vector& operator/=(const double& d) { for (unsigned long i = 0; i != rows; ++i) data[i] /= d; return *this; }

        friend vector operator+(const vector& v1, const vector& v2);
        friend vector operator-(const vector& v1, const vector& v2);
        friend vector operator*(const vector& v1, const vector& v2);
        friend vector operator/(const vector& v1, const vector& v2);
        friend vector operator+(const vector& v);
        friend vector operator-(const vector& v);

        inline double dot(const vector& v) const { check_dimensions(*this, v); double res = 0.0; for (unsigned long i = 0; i != rows; ++i) res += data[i]*v.data[i]; return res; }
        friend double dot(const vector& v1, const vector& v2);

        inline vector& operator+=(const vector& v) { check_dimensions(*this, v); for (unsigned long i = 0; i != rows; ++i) data[i] += v.data[i]; return *this; }
        inline vector& operator-=(const vector& v) { check_dimensions(*this, v); for (unsigned long i = 0; i != rows; ++i) data[i] -= v.data[i]; return *this; }
        inline vector& operator*=(const vector& v) { check_dimensions(*this, v); for (unsigned long i = 0; i != rows; ++i) data[i] *= v.data[i]; return *this; }
        inline vector& operator/=(const vector& v) { check_dimensions(*this, v); for (unsigned long i = 0; i != rows; ++i) data[i] /= v.data[i]; return *this; }

        friend vector abs(const vector& v);
        friend vector exp(const vector& v);
        friend vector exp2(const vector& v);
        friend vector log(const vector& v);
        friend vector log10(const vector& v);
        friend vector log2(const vector& v);
        friend vector pow(const vector& v, const double& d);
        friend vector pow(const vector& v1, const vector& v2);
        friend vector sqrt(const vector& v);
        friend vector sin(const vector& v);
        friend vector cos(const vector& v);
        friend vector tan(const vector& v);
        friend vector asin(const vector& v);
        friend vector acos(const vector& v);
        friend vector atan(const vector& v);

        inline void resize(const unsigned long& new_size)
        {
            if (new_size == rows) return;
            double* old_data = data;
            if (new_size)
            {
                data = new double[new_size];
                if (old_data) std::memcpy(data, old_data, std::min(rows, new_size) * sizeof(double));
                if (new_size > rows) std::memset(data + rows, 0, (new_size - rows) * sizeof(double));
            }
            else data = nullptr;
            if (old_data) delete[] old_data;
            this->rows = new_size;
        }

        inline void append(const double& d) { this->resize(this->rows+1); this->data[this->rows-1] = d; }
    };

    inline vector operator+(const vector& v, const double& d) { vector res(v); for (unsigned long i = 0; i != v.rows; ++i) res.data[i] += d; return res; }
    inline vector operator+(const double& d, const vector& v) { return v+d; }
    inline vector operator*(const vector& v, const double& d) { vector res(v); for (unsigned long i = 0; i != v.rows; ++i) res.data[i] *= d; return res; }
    inline vector operator*(const double& d, const vector& v) { return v*d; }
    inline vector operator-(const vector& v, const double& d) { vector res(v); for (unsigned long i = 0; i != v.rows; ++i) res.data[i] -= d; return res; }
    inline vector operator-(const double& d, const vector& v) { vector res(v); for (unsigned long i = 0; i != v.rows; ++i) res.data[i] = d - res.data[i]; return res; }
    inline vector operator/(const vector& v, const double& d) { vector res(v); for (unsigned long i = 0; i != v.rows; ++i) res.data[i] /= d; return res; }
    inline vector operator/(const double& d, const vector& v) { vector res(v); for (unsigned long i = 0; i != v.rows; ++i) res.data[i] = d / res.data[i]; return res; }
    
    inline void check_dimensions(const container& c1, const container& c2) { if (c1.rows != c2.rows || c1.cols != c2.cols) throw std::invalid_argument("v2"); }

    inline vector operator+(const vector& v1, const vector& v2) { check_dimensions(v1, v2); vector res(v1); for (unsigned long i = 0; i != res.rows; ++i) res.data[i] += v2.data[i]; return res; }
    inline vector operator-(const vector& v1, const vector& v2) { check_dimensions(v1, v2); vector res(v1); for (unsigned long i = 0; i != res.rows; ++i) res.data[i] -= v2.data[i]; return res; }
    inline vector operator*(const vector& v1, const vector& v2) { check_dimensions(v1, v2); vector res(v1); for (unsigned long i = 0; i != res.rows; ++i) res.data[i] *= v2.data[i]; return res; }
    inline vector operator/(const vector& v1, const vector& v2) { check_dimensions(v1, v2); vector res(v1); for (unsigned long i = 0; i != res.rows; ++i) res.data[i] /= v2.data[i]; return res; }
    inline vector operator+(const vector& v) { return vector(v); }
    inline vector operator-(const vector& v) { vector res(v); for (unsigned long i = 0; i != res.rows; ++i) res.data[i] = -res.data[i]; return res; }

    inline double dot(const vector& v1, const vector& v2) { return v1.dot(v2); }

    vector abs(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::abs(v[i]); return w; }
    vector exp(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::exp(v[i]); return w; }
    vector exp2(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::exp2(v[i]); return w; }
    vector log(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::log(v[i]); return w; }
    vector log10(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::log10(v[i]); return w; }
    vector log2(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::log2(v[i]); return w; }
    vector pow(const vector& v, const double& d) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::pow(v[i], d); return w; }
    vector pow(const vector& v1, const vector& v2) { check_dimensions(v1, v2); vector w(v1.rows); for (unsigned long i = 0; i != v1.rows; ++i) w[i] = std::pow(v1[i], v2[i]); return w; }
    vector sqrt(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::sqrt(v[i]); return w; }
    vector sin(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::sin(v[i]); return w; }
    vector cos(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::cos(v[i]); return w; }
    vector tan(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::tan(v[i]); return w; }
    vector asin(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::asin(v[i]); return w; }
    vector acos(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::acos(v[i]); return w; }
    vector atan(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::atan(v[i]); return w; }
    
    /*class matrix : public container {
        inline void resize(const unsigned long& new_rows, const unsigned long& new_cols)
        {
            if (new_cols == cols)
            {
                if (new_rows == rows) return;

                rows = new_rows;
            }
        }
    };*/

}

#endif
