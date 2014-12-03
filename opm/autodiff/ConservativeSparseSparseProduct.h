// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2008-2011 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef EIGEN_CONSERVATIVESPARSESPARSEPRODUCT_H
#define EIGEN_CONSERVATIVESPARSESPARSEPRODUCT_H

#warning "Using overloaded Eigen::ConservativeSparseSparseProduct.h"

#include <algorithm>
#include <iterator>
#include <functional>
#include <limits>
#include <vector>

#include <Eigen/Core>

namespace Eigen {

// forward declaration of SparseMatrix
template<typename _Scalar, int _Options, typename _Index>
class SparseMatrix;


namespace internal {

template < unsigned int depth >
struct QuickSort
{
  template <typename T>
  static inline void sort(T begin, T end)
  {
    if (begin != end)
    {
      T middle = std::partition (begin, end,
                                 std::bind2nd(std::less<typename std::iterator_traits<T>::value_type>(), *begin)
                                );
      QuickSort< depth-1 >::sort(begin, middle);

      // std::sort (max(begin + 1, middle), end);
      T new_middle = begin;
      QuickSort< depth-1 >::sort(++new_middle, end);
    }
  }
};

template <>
struct QuickSort< 0 >
{
  template <typename T>
  static inline void sort(T begin, T end)
  {
    // fall back to standard insertion sort
    std::sort( begin, end );
  }
};


template<typename Lhs, typename Rhs, typename ResultType>
static void conservative_sparse_sparse_product_impl(const Lhs& lhs, const Rhs& rhs, ResultType& res)
{
  // if one of the matrices does not contain non zero elements
  // the result will only contain an empty matrix
  if( lhs.nonZeros() == 0 || rhs.nonZeros() == 0 )
    return ;

  typedef typename remove_all<Lhs>::type::Scalar Scalar;
  typedef typename remove_all<Lhs>::type::Index Index;

  // make sure to call innerSize/outerSize since we fake the storage order.
  Index rows = lhs.innerSize();
  Index cols = rhs.outerSize();
  eigen_assert(lhs.outerSize() == rhs.innerSize());

  std::vector<bool> mask(rows,false);
  Matrix<Scalar,Dynamic,1> values(rows);
  Matrix<Index,Dynamic,1>  indices(rows);

  // estimate the number of non zero entries
  // given a rhs column containing Y non zeros, we assume that the respective Y columns
  // of the lhs differs in average of one non zeros, thus the number of non zeros for
  // the product of a rhs column with the lhs is X+Y where X is the average number of non zero
  // per column of the lhs.
  // Therefore, we have nnz(lhs*rhs) = nnz(lhs) + nnz(rhs)
  Index estimated_nnz_prod = lhs.nonZeros() + rhs.nonZeros();

  res.setZero();
  res.reserve(Index(estimated_nnz_prod));

  //const Scalar epsilon = std::numeric_limits< Scalar >::epsilon();
  const Scalar epsilon = 1e-15 ;

  // we compute each column of the result, one after the other
  for (Index j=0; j<cols; ++j)
  {
    Index nnz = 0;
    for (typename Rhs::InnerIterator rhsIt(rhs, j); rhsIt; ++rhsIt)
    {
      const Scalar y = rhsIt.value();
      for (typename Lhs::InnerIterator lhsIt(lhs, rhsIt.index()); lhsIt; ++lhsIt)
      {
        const Index i = lhsIt.index();
        const Scalar val = lhsIt.value() * y;
        if( std::abs( val ) > epsilon )
        {
          if(!mask[i])
          {
            mask[i] = true;
            values[i] = val;
            indices[nnz] = i;
            ++nnz;
          }
          else
            values[i] += val;
        }
      }
    }

    if( nnz > 1 )
    {
      // sort indices for sorted insertion to avoid later copying
      QuickSort< 1 >::sort( indices.data(), indices.data()+nnz );
    }

    res.startVec(j);
    // ordered insertion
    // still using insertBackByOuterInnerUnordered since we know what we are doing
    for(Index k=0; k<nnz; ++k)
    {
      const Index i = indices[k];
      res.insertBackByOuterInnerUnordered(j,i) = values[i];
      mask[i] = false;
    }

#if 0
    // alternative ordered insertion code:

    Index t200 = rows/(log2(200)*1.39);
    Index t = (rows*100)/139;

    // FIXME reserve nnz non zeros
    // FIXME implement fast sort algorithms for very small nnz
    // if the result is sparse enough => use a quick sort
    // otherwise => loop through the entire vector
    // In order to avoid to perform an expensive log2 when the
    // result is clearly very sparse we use a linear bound up to 200.
    //if((nnz<200 && nnz<t200) || nnz * log2(nnz) < t)
    //res.startVec(j);
    if(true)
    {
      if(nnz>1) std::sort(indices.data(),indices.data()+nnz);
      for(Index k=0; k<nnz; ++k)
      {
        Index i = indices[k];
        res.insertBackByOuterInner(j,i) = values[i];
        mask[i] = false;
      }
    }
    else
    {
      // dense path
      for(Index i=0; i<rows; ++i)
      {
        if(mask[i])
        {
          mask[i] = false;
          res.insertBackByOuterInner(j,i) = values[i];
        }
      }
    }
#endif

  }
  res.finalize();
}

} // end namespace internal

namespace internal {

template<typename Lhs, typename Rhs, typename ResultType,
  int LhsStorageOrder = (traits<Lhs>::Flags&RowMajorBit) ? RowMajor : ColMajor,
  int RhsStorageOrder = (traits<Rhs>::Flags&RowMajorBit) ? RowMajor : ColMajor,
  int ResStorageOrder = (traits<ResultType>::Flags&RowMajorBit) ? RowMajor : ColMajor>
struct conservative_sparse_sparse_product_selector;

template<typename Lhs, typename Rhs, typename ResultType>
struct conservative_sparse_sparse_product_selector<Lhs,Rhs,ResultType,ColMajor,ColMajor,ColMajor>
{
  typedef typename remove_all<Lhs>::type LhsCleaned;
  typedef typename LhsCleaned::Scalar Scalar;

  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res)
  {
    //typedef SparseMatrix<typename ResultType::Scalar,RowMajor,typename ResultType::Index> RowMajorMatrix;
    typedef SparseMatrix<typename ResultType::Scalar,ColMajor,typename ResultType::Index> ColMajorMatrix;
    //ColMajorMatrix resCol(lhs.rows(),rhs.cols());
    res = ColMajorMatrix(lhs.rows(),rhs.cols());
    internal::conservative_sparse_sparse_product_impl<Lhs,Rhs,ColMajorMatrix>(lhs, rhs, res);
    //internal::conservative_sparse_sparse_product_impl<Lhs,Rhs,ColMajorMatrix>(lhs, rhs, resCol);
    // sort the non zeros:
    //RowMajorMatrix resRow(resCol);
    //res = resRow;
  }
};

template<typename Lhs, typename Rhs, typename ResultType>
struct conservative_sparse_sparse_product_selector<Lhs,Rhs,ResultType,RowMajor,ColMajor,ColMajor>
{
  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res)
  {
     typedef SparseMatrix<typename ResultType::Scalar,ColMajor,typename ResultType::Index> ColMajorMatrix;
     //RowMajorMatrix rhsRow = rhs;
     //RowMajorMatrix resRow(lhs.rows(), rhs.cols());
     ColMajorMatrix lhsCol = lhs;
     res = ResultType( lhs.rows(), rhs.cols() );
     internal::conservative_sparse_sparse_product_impl<ColMajorMatrix, Rhs, ResultType>( lhsCol, rhs, res );
     //internal::conservative_sparse_sparse_product_impl<RowMajorMatrix,Lhs,RowMajorMatrix>(rhsRow, lhs, resRow);
     //res = resRow;
  }
};

template<typename Lhs, typename Rhs, typename ResultType>
struct conservative_sparse_sparse_product_selector<Lhs,Rhs,ResultType,ColMajor,RowMajor,ColMajor>
{
  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res)
  {
    typedef SparseMatrix<typename ResultType::Scalar,ColMajor,typename ResultType::Index> ColMajorMatrix;
    ColMajorMatrix rhsCol = rhs;
    res = ResultType( lhs.rows(), rhs.cols() );
    internal::conservative_sparse_sparse_product_impl<Lhs, ColMajorMatrix, ResultType>( lhs, rhsCol, res);
    /*
    typedef SparseMatrix<typename ResultType::Scalar,RowMajor,typename ResultType::Index> RowMajorMatrix;
    RowMajorMatrix lhsRow = lhs;
    RowMajorMatrix resRow(lhs.rows(), rhs.cols());
    internal::conservative_sparse_sparse_product_impl<Rhs,RowMajorMatrix,RowMajorMatrix>(rhs, lhsRow, resRow);
    res = resRow;
    */
  }
};

template<typename Lhs, typename Rhs, typename ResultType>
struct conservative_sparse_sparse_product_selector<Lhs,Rhs,ResultType,RowMajor,RowMajor,ColMajor>
{
  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res)
  {
    typedef SparseMatrix<typename ResultType::Scalar,RowMajor,typename ResultType::Index> RowMajorMatrix;
    RowMajorMatrix resRow(lhs.rows(), rhs.cols());
    internal::conservative_sparse_sparse_product_impl<Rhs,Lhs,RowMajorMatrix>(rhs, lhs, resRow);
    res = resRow;
  }
};


template<typename Lhs, typename Rhs, typename ResultType>
struct conservative_sparse_sparse_product_selector<Lhs,Rhs,ResultType,ColMajor,ColMajor,RowMajor>
{
  typedef typename traits<typename remove_all<Lhs>::type>::Scalar Scalar;

  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res)
  {
    typedef SparseMatrix<typename ResultType::Scalar,ColMajor,typename ResultType::Index> ColMajorMatrix;
    ColMajorMatrix resCol(lhs.rows(), rhs.cols());
    internal::conservative_sparse_sparse_product_impl<Lhs,Rhs,ColMajorMatrix>(lhs, rhs, resCol);
    res = resCol;
  }
};

template<typename Lhs, typename Rhs, typename ResultType>
struct conservative_sparse_sparse_product_selector<Lhs,Rhs,ResultType,RowMajor,ColMajor,RowMajor>
{
  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res)
  {
    typedef SparseMatrix<typename ResultType::Scalar,RowMajor,typename ResultType::Index> RowMajorMatrix;
    RowMajorMatrix rhsRow = rhs;
    res = ResultType( lhs.rows(), rhs.cols() );
    internal::conservative_sparse_sparse_product_impl<Lhs, RowMajorMatrix, ResultType>(rhsRow, lhs, res);
    /*
    typedef SparseMatrix<typename ResultType::Scalar,ColMajor,typename ResultType::Index> ColMajorMatrix;
    ColMajorMatrix lhsCol = lhs;
    ColMajorMatrix resCol(lhs.rows(), rhs.cols());
    internal::conservative_sparse_sparse_product_impl<ColMajorMatrix,Rhs,ColMajorMatrix>(lhsCol, rhs, resCol);
    res = resCol;
    */
  }
};

template<typename Lhs, typename Rhs, typename ResultType>
struct conservative_sparse_sparse_product_selector<Lhs,Rhs,ResultType,ColMajor,RowMajor,RowMajor>
{
  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res)
  {
    typedef SparseMatrix<typename ResultType::Scalar,RowMajor,typename ResultType::Index> RowMajorMatrix;
    RowMajorMatrix lhsRow = lhs;
    res = RowMajorMatrix( lhs.rows(), rhs.cols() );
    internal::conservative_sparse_sparse_product_impl<Rhs, RowMajorMatrix, ResultType>(rhs, lhsRow, res);
    /*
    typedef SparseMatrix<typename ResultType::Scalar,ColMajor,typename ResultType::Index> ColMajorMatrix;
    ColMajorMatrix rhsCol = rhs;
    ColMajorMatrix resCol(lhs.rows(), rhs.cols());
    internal::conservative_sparse_sparse_product_impl<Lhs,ColMajorMatrix,ColMajorMatrix>(lhs, rhsCol, resCol);
    res = resCol;
    */
  }
};

template<typename Lhs, typename Rhs, typename ResultType>
struct conservative_sparse_sparse_product_selector<Lhs,Rhs,ResultType,RowMajor,RowMajor,RowMajor>
{
  static void run(const Lhs& lhs, const Rhs& rhs, ResultType& res)
  {
    typedef SparseMatrix<typename ResultType::Scalar,RowMajor,typename ResultType::Index> RowMajorMatrix;
    //typedef SparseMatrix<typename ResultType::Scalar,ColMajor,typename ResultType::Index> ColMajorMatrix;
    res = RowMajorMatrix( lhs.rows(),rhs.cols() );
    //RowMajorMatrix resRow(lhs.rows(),rhs.cols());
    internal::conservative_sparse_sparse_product_impl<Rhs,Lhs,RowMajorMatrix>(rhs, lhs, res);
    // sort the non zeros:
    //ColMajorMatrix resCol(resRow);
    //res = resCol;
  }
};

} // end namespace internal

} // end namespace Eigen

#endif // EIGEN_CONSERVATIVESPARSESPARSEPRODUCT_H
