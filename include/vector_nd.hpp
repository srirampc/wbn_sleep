#ifndef VECTOR_ND_HPP
#define VECTOR_ND_HPP

#include <cassert>
#include <vector>

template <typename T> class Vector1d {
    //
    std::size_t _xdim, _nelts;
    T _init;
    std::vector<T> _rdata;

  public:
    Vector1d() : _xdim(0), _nelts(0), _init(0), _rdata() {}

    Vector1d(std::size_t xdim, T init)
        : _xdim(xdim), _nelts(xdim), _init(init), _rdata(_nelts, _init) {}

    inline std::size_t index(std::size_t i) const { return i; }

    inline T& operator()(std::size_t i) {
        std::size_t idx = index(i);
        assert(idx < _nelts);
        return _rdata[idx];
    }

    inline const T& operator()(std::size_t i) const {
        std::size_t idx = index(i);
        assert(idx < _nelts);
        return _rdata[idx];
    }

    inline const T& operator[](std::size_t i) const {
        std::size_t idx = index(i);
        assert(idx < _nelts);
        return _rdata[idx];
    }

    inline std::size_t size() const { return _nelts; }
    //
    inline T* data() { return _rdata.data(); }
    inline const T* data() const { return _rdata.data(); }
};

template <typename T> class Vector2d {
    //
    std::size_t _xdim, _ydim, _nelts;
    T _init;
    std::vector<T> _rdata;
    //
  public:
    Vector2d() : _xdim(0), _ydim(0), _nelts(0), _init(0), _rdata() {}

    Vector2d(std::size_t xdim, std::size_t ydim, T init)
        : _xdim(xdim), _ydim(ydim), _nelts(xdim * ydim), _init(init),
          _rdata(_nelts, _init) {}

    inline std::size_t index(std::size_t i, std::size_t j) const {
        return (i * _ydim + j);  // Row major
    }

    inline T& operator()(std::size_t i, std::size_t j) {
        std::size_t idx = index(i, j);
        assert(idx < _nelts);
        return _rdata[idx];
    }

    inline const T& operator()(std::size_t i, std::size_t j) const {
        std::size_t idx = index(i, j);
        assert(idx < _nelts);
        return _rdata[idx];
    }

    inline std::size_t size() const { return _nelts; }
    //
    inline T* data() { return _rdata.data(); }
    inline const T* data() const { return _rdata.data(); }
};

template <typename T> class Vector3d {
    //
    std::size_t _xdim, _ydim, _zdim, _nelts;
    T _init;
    std::vector<T> _rdata;
    std::size_t _rdelta, _sdelta;
    //
  public:
    Vector3d() : _xdim(0), _ydim(0), _zdim(0), _nelts(0), _init(0), _rdata() {}

    Vector3d(std::size_t xdim, std::size_t ydim, std::size_t zdim, T init)
        : _xdim(xdim), _ydim(ydim), _zdim(zdim), _nelts(xdim * ydim * zdim),
          _init(init), _rdata(_nelts, _init), _rdelta(_ydim * _zdim),
          _sdelta(_zdim) {}

    inline std::size_t index(std::size_t i, std::size_t j,
                             std::size_t k) const {
        return (i * _rdelta + j * _sdelta + k);  // x-major followed by y-major
    }

    inline T& operator()(std::size_t i, std::size_t j, std::size_t k) {
        std::size_t idx = index(i, j, k);
        assert(idx < _nelts);
        return _rdata[idx];
    }

    inline const T& operator()(std::size_t i, std::size_t j,
                               std::size_t k) const {
        std::size_t idx = index(i, j, k);
        assert(idx < _nelts);
        return _rdata[idx];
    }

    inline std::size_t size() const { return _nelts; }
    //
    inline T* data() { return _rdata.data(); }
    inline const T* data() const { return _rdata.data(); }
};

template <typename T> class Vector4d {
    //
    std::size_t _wdim, _xdim, _ydim, _zdim, _nelts;
    T _init;
    std::vector<T> _rdata;
    std::size_t _rdelta, _sdelta, _tdelta;
    //
  public:
    Vector4d()
        : _wdim(0), _xdim(0), _ydim(0), _zdim(0), _nelts(0), _init(0),
          _rdata() {}

    Vector4d(std::size_t wdim, std::size_t xdim, std::size_t ydim,
             std::size_t zdim, T init)
        : _wdim(wdim), _xdim(xdim), _ydim(ydim), _zdim(zdim),
          _nelts(wdim * xdim * ydim * zdim), _init(init), _rdata(_nelts, _init),
          _rdelta(_xdim * _ydim * _zdim), _sdelta(_ydim * _zdim),
          _tdelta(_zdim) {}

    inline std::size_t index(std::size_t i, std::size_t j, std::size_t k,
                             std::size_t l) const {
        return (i * _rdelta + j * _sdelta + k * _tdelta + l);
    }

    inline T& operator()(std::size_t i, std::size_t j, std::size_t k,
                         std::size_t l) {
        std::size_t idx = index(i, j, k, l);
        assert(idx < _nelts);
        return _rdata[idx];
    }

    inline const T& operator()(std::size_t i, std::size_t j, std::size_t k,
                               std::size_t l) const {
        std::size_t idx = index(i, j, k, l);
        assert(idx < _nelts);
        return _rdata[idx];
    }

    inline std::size_t size() const { return _nelts; }
    //
    inline T* data() { return _rdata.data(); }
    inline const T* data() const { return _rdata.data(); }
};

#endif // VECTOR_ND_HPP
