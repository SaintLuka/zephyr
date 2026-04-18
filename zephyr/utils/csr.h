#pragma once
#include <vector>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <numeric>
#include <span>

namespace zephyr {

/// @brief CSR container that mimics std::vector<std::vector<T>>:
///        - csr[i] returns a row view (iterable, indexable)
///        - stores all elements in a single packed array
///
/// Typical build pattern:
///   Csr<T> csr;
///   csr.assign_row_sizes(sizes);     // builds offsets and allocates values
///   // fill:
///   for (Index i=0; i < csr.n_rows(); ++i)
///     for (Index k=0; k<csr.row_size(i); ++k) csr[i][k] = ...
template <class T, class Index = std::size_t>
class Csr {
    static_assert(std::is_integral_v<Index>, "Csr<T,Index>: Index must be integral.");
public:
    using value_type = T;
    using index_type = Index;
    using size_type  = std::size_t;

    std::vector<Index> offsets;  ///< size = n_rows()+1, offsets[0]=0, offsets.back()=values.size()
    std::vector<T>     values;   ///< packed row data

    Csr() = default;

    bool empty() const noexcept {
        return values.empty();
    }

    size_type size() const noexcept {
        return n_rows();
    }

    void clear() {
        offsets.clear();
        values.clear();
    }

    void shrink_to_fit() {
        offsets.shrink_to_fit();
        values.shrink_to_fit();
    }

    size_type n_rows() const noexcept {
        return offsets.empty() ? 0 : (offsets.size() - 1);
    }

    size_type nnz() const noexcept {
        return values.size();
    }

    /// @brief Allocate CSR structure from per-row sizes (like vector<vector<T>> with known sizes).
    ///        This builds offsets and resizes values to sum(sizes).
    void assign_row_sizes(const std::vector<Index>& row_sizes) {
        offsets.resize(row_sizes.size() + 1);
        offsets[0] = 0;
        for (size_type i = 0; i < row_sizes.size(); ++i) {
            offsets[i + 1] = offsets[i] + row_sizes[i];
        }
        values.resize(static_cast<size_type>(offsets.back()));
    }

    /// @brief Resize to n_rows with all rows empty (values cleared).
    void resize_empty(size_type n_rows) {
        offsets.assign(n_rows + 1, Index{0});
        values.clear();
    }

    /// @brief Row access (non-const), like vector<vector<T>>::operator[].
    std::span<T> operator[](size_type row) {
        if (row + 1 >= offsets.size()) throw std::out_of_range("Csr::operator[]: row out of range");
        return {values.data() + offsets[row], static_cast<size_t>(offsets[row + 1] - offsets[row])};
    }

    /// @brief Row access (const).
    std::span<const T> operator[](size_type row) const {
        if (row + 1 >= offsets.size()) throw std::out_of_range("Csr::operator[] const: row out of range");
        return {values.data() + offsets[row], static_cast<size_t>(offsets[row + 1] - offsets[row])};
    }

    /// @brief Number of items in row.
    Index row_size(size_type row) const {
        if (row + 1 >= offsets.size()) throw std::out_of_range("Csr::row_size: row out of range");
        return offsets[row + 1] - offsets[row];
    }

    /// @brief Validate CSR invariants (throws on error).
    void validate_or_throw() const {
        if (offsets.empty()) return;
        if (offsets.size() < 2) throw std::logic_error("Csr: offsets must have size >= 2");
        if (offsets.front() != 0) throw std::logic_error("Csr: offsets[0] must be 0");
        for (size_type i = 0; i + 1 < offsets.size(); ++i) {
            if (offsets[i] > offsets[i + 1]) throw std::logic_error("Csr: offsets must be non-decreasing");
        }
        if (static_cast<size_type>(offsets.back()) != values.size()) {
            throw std::logic_error("Csr: offsets.back() must equal values.size()");
        }
    }
};

} // namespace zephyr