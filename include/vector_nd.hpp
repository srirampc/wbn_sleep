#ifndef VECTOR_ND_HPP
#define VECTOR_ND_HPP

#include <cassert>
#include <cstdlib>
#include <vector>

template <typename T> class VectorInterface {
  protected:
    std::size_t _nelts;
    T _init;
    std::vector<T> _rdata;

  public:
    VectorInterface() : _init(0), _rdata(0) {}
    VectorInterface(std::size_t n, T init)
        : _nelts(n), _init(init), _rdata(_nelts, _init) {}
    //
    inline std::size_t size() const { return _nelts; }
    inline T* data() { return _rdata.data(); }
    inline const T* data() const { return _rdata.data(); }
    //
    inline const std::vector<T>& data_vec() { return _rdata; }
};

template <typename T> class Vector1d : public VectorInterface<T> {
    //
    std::size_t _xdim;

  public:
    Vector1d() : VectorInterface<T>(), _xdim(0) {}

    Vector1d(std::size_t xdim, T init)
        : VectorInterface<T>(xdim, init), _xdim(xdim) {}

    inline std::size_t index(std::size_t i) const { return i; }

    inline void set(std::size_t i, T value) {
        std::size_t idx = index(i);
        assert(idx < this->_nelts);
        this->_rdata[idx] = value;
    }

    inline T& operator()(std::size_t i) {
        std::size_t idx = index(i);
        assert(idx < this->_nelts);
        return this->_rdata[idx];
    }

    inline const T& operator()(std::size_t i) const {
        std::size_t idx = index(i);
        assert(idx < this->_nelts);
        return this->_rdata[idx];
    }

    inline const T& operator[](std::size_t i) const { return operator()(i); }
};

template <typename T> class Vector2d : public VectorInterface<T> {
    //
    std::size_t _xdim, _ydim;
    //
  public:
    Vector2d() : VectorInterface<T>(), _xdim(0), _ydim(0) {}

    Vector2d(std::size_t xdim, std::size_t ydim, T init)
        : VectorInterface<T>(xdim * ydim, init), _xdim(xdim), _ydim(ydim) {}

    inline std::size_t index(std::size_t i, std::size_t j) const {
        return (i * _ydim + j);  // Row major
    }

    inline void set(std::size_t i, std::size_t j, T value) {
        std::size_t idx = index(i, j);
        assert(idx < this->_nelts);
        this->_rdata[idx] = value;
    }

    inline T& operator()(std::size_t i, std::size_t j) {
        std::size_t idx = index(i, j);
        assert(idx < this->_nelts);
        return this->_rdata[idx];
    }

    inline const T& operator()(std::size_t i, std::size_t j) const {
        std::size_t idx = index(i, j);
        assert(idx < this->_nelts);
        return this->_rdata[idx];
    }

    inline int parse_record1(const char*& c_ptr) {  // NOLINT
        char* end;
        int i = std::strtol(c_ptr, &end, 10);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        int j = std::strtol(c_ptr, &end, 10);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        T value = std::strtof(c_ptr, &end);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        set(--i, --j, value);
        return 0;
    }
};

template <typename T> class Vector3d : VectorInterface<T> {
    //
    std::size_t _xdim, _ydim, _zdim;
    std::size_t _rdelta, _sdelta;
    //
  public:
    Vector3d() : VectorInterface<T>(), _xdim(0), _ydim(0), _zdim(0) {}

    Vector3d(std::size_t xdim, std::size_t ydim, std::size_t zdim, T init)
        : VectorInterface<T>(xdim * ydim * zdim, init), _xdim(xdim),
          _ydim(ydim), _zdim(zdim), _rdelta(_ydim * _zdim), _sdelta(_zdim) {}

    inline std::size_t index(std::size_t i, std::size_t j,
                             std::size_t k) const {
        return (i * _rdelta + j * _sdelta + k);  // x-major followed by y-major
    }

    inline void set(std::size_t i, std::size_t j, std::size_t k, T value) {
        std::size_t idx = index(i, j, k);
        assert(idx < this->_nelts);
        this->_rdata[idx] = value;
    }

    inline T& operator()(std::size_t i, std::size_t j, std::size_t k) {
        std::size_t idx = index(i, j, k);
        assert(idx < this->_nelts);
        return this->_rdata[idx];
    }

    inline const T& operator()(std::size_t i, std::size_t j,
                               std::size_t k) const {
        std::size_t idx = index(i, j, k);
        assert(idx < this->_nelts);
        return this->_rdata[idx];
    }

    inline int parse_record1(const char*& c_ptr) {  // NOLINT
        char* end;
        int i = std::strtol(c_ptr, &end, 10);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        int j = std::strtol(c_ptr, &end, 10);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        int k = std::strtol(c_ptr, &end, 10);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        T value = std::strtof(c_ptr, &end);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        set(--i, --j, --k, value);
        return 0;
    }
};

template <typename T> class Vector4d : public VectorInterface<T> {
    //
    std::size_t _wdim, _xdim, _ydim, _zdim;
    std::size_t _rdelta, _sdelta, _tdelta;
    //
  public:
    Vector4d() : VectorInterface<T>(), _wdim(0), _xdim(0), _ydim(0), _zdim(0) {}

    Vector4d(std::size_t wdim, std::size_t xdim, std::size_t ydim,
             std::size_t zdim, T init)
        : VectorInterface<T>(wdim * xdim * ydim * zdim, init), _wdim(wdim),
          _xdim(xdim), _ydim(ydim), _zdim(zdim), _rdelta(_xdim * _ydim * _zdim),
          _sdelta(_ydim * _zdim), _tdelta(_zdim) {}

    inline std::size_t index(std::size_t i, std::size_t j, std::size_t k,
                             std::size_t l) const {
        return (i * _rdelta + j * _sdelta + k * _tdelta + l);
    }

    inline void set(std::size_t i, std::size_t j, std::size_t k, std::size_t l,
                    T value) {
        std::size_t idx = index(i, j, k, l);
        assert(idx < this->_nelts);
        this->_rdata[idx] = value;
    }
    inline T& operator()(std::size_t i, std::size_t j, std::size_t k,
                         std::size_t l) {
        std::size_t idx = index(i, j, k, l);
        assert(idx < this->_nelts);
        return this->_rdata[idx];
    }

    inline const T& operator()(std::size_t i, std::size_t j, std::size_t k,
                               std::size_t l) const {
        std::size_t idx = index(i, j, k, l);
        assert(idx < this->_nelts);
        return this->_rdata[idx];
    }

    inline int parse_record1(const char*& c_ptr) {  // NOLINT
        char* end;
        int i = std::strtol(c_ptr, &end, 10);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        int j = std::strtol(c_ptr, &end, 10);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        int k = std::strtol(c_ptr, &end, 10);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        int l = std::strtol(c_ptr, &end, 10);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        T value = std::strtof(c_ptr, &end);
        if (c_ptr == end) {
            return -1;
        }
        c_ptr = end;
        set(--i, --j, --k, --l, value);
        return 0;
    }
};

template <typename VType>
inline std::size_t load_vector(VType* pvec_nd, const char*& c_ptr) {  // NOLINT
    std::size_t nrecords = 0;
    do {
        if (pvec_nd->parse_record1(c_ptr) != 0)
            break;
        nrecords++;
    } while (1);  // while
    return nrecords;
}

template <typename VType>
inline std::size_t load_vector_until(VType* pvec_nd,
                                     const char*& c_ptr,  // NOLINT
                                     const char* end_ptr) {
    std::size_t nrecords = 0;
    do {
        if (pvec_nd->parse_record1(c_ptr) != 0)
            break;
        nrecords++;
        if (c_ptr >= end_ptr)
            break;
    } while (1);  // while
    return nrecords;
}

#endif  // VECTOR_ND_HPP
