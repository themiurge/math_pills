#ifndef CORE_HPP
#define CORE_HPP

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <functional>
#include <initializer_list>
#include <random>

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
    protected:
        container(const container& c1, const container& c2, std::function<double(const double&, const double&)> op)
        {
            this->rows = c1.rows; this->cols = c1.cols;
            if (rows*cols)
            {
                this->data = new double[rows*cols];
                for (unsigned long i = 0; i != rows*cols; ++i)
                    this->data[i] = op(c1.data[i], c2.data[i]);
            }
            else
                this->data = nullptr;
        }
        container(const container& c, std::function<double(const double&)> op)
        {
            this->rows = c.rows; this->cols = c.cols;
            if (rows*cols)
            {
                this->data = new double[rows*cols];
                for (unsigned long i = 0; i != rows*cols; ++i)
                    this->data[i] = op(c.data[i]);
            }
            else
                this->data = nullptr;
        }
    public:

        typedef double* iterator;
        typedef const double* const_iterator;
        iterator begin() const { return &data[0]; }
        iterator end() const { return &data[rows*cols]; }
        const_iterator cbegin() const { return &data[0]; }
        const_iterator cend() const { return &data[rows*cols]; }

        container() { this->data = nullptr; this->rows = 0; this->cols = 0; }
        container(const unsigned long& rows, const unsigned long& cols, const double& d) { this->data = ((rows*cols) ? new double[rows*cols] : nullptr); this->rows = rows; this->cols = cols; fill(d); }
        container(const container& c) { this->copy_data(c.data, c.rows, c.cols); }
        container(container&& c) { this->move_data(c.data, c.rows, c.cols); }
        container(std::initializer_list<std::initializer_list<double>> l) 
        { 
            rows = l.size(); 
            if (rows && (cols = l.begin()->size()))
            {
                data = new double[rows*cols];
                unsigned long pos = 0;
                for (auto& ll : l)
                {
                    if (ll.size() != cols) 
                    {
                        delete[] data;
                        throw std::invalid_argument("l");
                    }
                    for (auto d : ll)
                        data[pos++] = d;
                }
            }
            else
            {
                cols = 0;
                data = nullptr;
            }
        }

        inline void fill(const double& d) { for (auto& elem : *this) elem = d; }

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
    protected:
        vector(const vector& v1, const vector& v2, std::function<double(const double&, const double&)> op) : container(v1, v2, op) { }
        vector(const vector& v, std::function<double(const double&)> op) : container(v, op) { }
    public:
        vector() : container(0, 1, 0.0) { }
        vector(const unsigned long& dim) : container(dim, 1, 0.0) { }
        vector(const unsigned long& dim, const double& fill_value) : container(dim, 1, fill_value) { }
        vector(const vector& v) { this->copy_data(v.data, v.rows, v.cols); }
        vector(vector&& v) { this->move_data(v.data, v.rows, v.cols); }
        vector(const double* d, const unsigned long& start, const unsigned long& count) 
        { 
            cols = 1; 
            if ((rows = count) > 0) 
            { 
                data = new double[rows]; 
                std::memcpy(data, d + start, count * sizeof(double)); 
            } 
            else 
                data = nullptr; 
        }
        vector(std::initializer_list<double> l)
        {
            rows = l.size();
            cols = 1;
            if (rows)
            {
                data = new double[rows];
                unsigned long pos = 0;
                for (auto d : l)
                    data[pos++] = d;
            }
            else
                data = nullptr;
        }

        static vector from(const vector& v1, const vector& v2, std::function<double(const double&, const double&)> op) { check_dimensions(v1, v2); return vector(v1, v2, op); }
        static vector from(const vector& v, std::function<double(const double&)> op) { return vector(v, op); }

        static vector linspace(const double& start, const double& end, const unsigned long& n) 
        { 
            vector v;
            v.cols = 1; 
            if ((v.rows = n) > 0)
            {
                v.data = new double[n];
                if (n == 1) v.data[0] = (end - start) / 2.0;
                else 
                {
                    double interval = (end - start) / (double)(n-1);
                    for (unsigned long i = 0; i != n; ++i)
                        v.data[i] = start + interval * (double)i;
                }
            }
            else
                v.data = nullptr;
            return v;
        }

        static vector randn(const unsigned long& n, const double& mean = 0.0, const double& sigma = 1.0)
        {
            vector v;
            v.cols = 1;
            if ((v.rows = n) > 0)
            {
                std::random_device rd;
                std::mt19937 gen(rd());

                std::normal_distribution<double> d(mean, sigma);
                v.data = new double[v.rows];
                for (auto& elem : v)
                    elem = d(gen);
            }
            else
                v.data = nullptr;
            return v;
        }

        inline const unsigned long& size() const { return rows; }

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
        inline vector& operator+=(const double& d) { for (auto& elem : *this) elem += d; return *this; }
        inline vector& operator-=(const double& d) { for (auto& elem : *this) elem -= d; return *this; }
        inline vector& operator*=(const double& d) { for (auto& elem : *this) elem *= d; return *this; }
        inline vector& operator/=(const double& d) { for (auto& elem : *this) elem /= d; return *this; }

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

        inline double mean() const { if (!rows) return 0.0; double m = 0.0; for (const auto& d : *this) m += d; return m / (double)rows; }
        inline double var() const { if (!rows) return 0.0; double m = mean(), v = 0.0; for (const auto& d : *this) v += (d-m)*(d-m); return v / (double)rows; }
        inline double stddev() const { return std::sqrt(var()); }
    };

    inline vector operator+(const vector& v, const double& d) { vector res(v); for (auto& elem : res) elem += d; return res; }
    inline vector operator+(const double& d, const vector& v) { return v+d; }
    inline vector operator*(const vector& v, const double& d) { vector res(v); for (auto& elem : res) elem *= d; return res; }
    inline vector operator*(const double& d, const vector& v) { return v*d; }
    inline vector operator-(const vector& v, const double& d) { vector res(v); for (auto& elem : res) elem -= d; return res; }
    inline vector operator-(const double& d, const vector& v) { vector res(v); for (auto& elem : res) elem = d - elem; return res; }
    inline vector operator/(const vector& v, const double& d) { vector res(v); for (auto& elem : res) elem /= d; return res; }
    inline vector operator/(const double& d, const vector& v) { vector res(v); for (auto& elem : res) elem = d / elem; return res; }
    
    inline void check_dimensions(const container& c1, const container& c2) { if (c1.rows != c2.rows || c1.cols != c2.cols) throw std::invalid_argument("v2"); }

    inline vector operator+(const vector& v1, const vector& v2) { return vector::from(v1, v2, [](const double& d1, const double& d2) { return d1 + d2; }); }
    inline vector operator-(const vector& v1, const vector& v2) { return vector::from(v1, v2, [](const double& d1, const double& d2) { return d1 - d2; }); }
    inline vector operator*(const vector& v1, const vector& v2) { return vector::from(v1, v2, [](const double& d1, const double& d2) { return d1 * d2; }); }
    inline vector operator/(const vector& v1, const vector& v2) { return vector::from(v1, v2, [](const double& d1, const double& d2) { return d1 / d2; }); }
    inline vector operator+(const vector& v) { return vector(v); }
    inline vector operator-(const vector& v) { vector res(v); for (auto& elem : res) elem = -elem; return res; }

    inline double dot(const vector& v1, const vector& v2) { return v1.dot(v2); }

    vector abs(const vector& v) { return vector::from(v, [](const double& d) { return std::abs(d); }); }
    vector exp(const vector& v) { return vector::from(v, [](const double& d) { return std::exp(d); }); }
    vector exp2(const vector& v) { return vector::from(v, [](const double& d) { return std::exp2(d); }); }
    vector log(const vector& v) { return vector::from(v, [](const double& d) { return std::log(d); }); }
    vector log10(const vector& v) { return vector::from(v, [](const double& d) { return std::log10(d); }); }
    vector log2(const vector& v) { return vector::from(v, [](const double& d) { return std::log2(d); }); }
    vector pow(const vector& v, const double& d) { return vector::from(v, [&](const double& dd) { return std::pow(dd, d); }); }
    vector pow(const vector& v1, const vector& v2) { return vector::from(v1, v2, [](const double& d1, const double& d2) { return std::pow(d1, d2); }); }
    vector sqrt(const vector& v) { return vector::from(v, [](const double& d) { return std::sqrt(d); }); }
    vector sin(const vector& v) { return vector::from(v, [](const double& d) { return std::sin(d); }); }
    vector cos(const vector& v) { return vector::from(v, [](const double& d) { return std::cos(d); }); }
    vector tan(const vector& v) { return vector::from(v, [](const double& d) { return std::tan(d); }); }
    vector asin(const vector& v) { return vector::from(v, [](const double& d) { return std::asin(d); }); }
    vector acos(const vector& v) { return vector::from(v, [](const double& d) { return std::acos(d); }); }
    vector atan(const vector& v) { return vector::from(v, [](const double& d) { return std::atan(d); }); }

    std::ostream& operator<<(std::ostream& out, const vector& v)
    {
        out << '[';
        for (unsigned long i = 0; i != v.size(); ++i) out << (i ? ", " : "") << v[i];
        out << ']';
        return out;
    }
    
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
